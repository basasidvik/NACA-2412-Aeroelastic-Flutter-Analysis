function [U_vec, freq_Hz, damping, A_hist, divergence_indicator, ...
          max_growth, dominant_eigenvalue, real_eigenvalues] = ...
          flutter_sweep(Ms, Ks, free_dofs, node_y, dof_map, c_dist, ...
                        ea_loc_frac, a_ea, AR, e_s, aero, rho, U_vec, ...
                        n_modes_keep, Cs)
% FLUTTER_SWEEP  Aeroelastic speed sweep with MAC-based modal tracking.
%
% Oscillatory modes are represented by the member of each conjugate pair
% with positive imaginary part. Branches are continued using eigenvector MAC
% correlation, not by sorting frequency at every speed. If a branch becomes
% real/overdamped, its oscillatory output is NaN at that speed rather than
% being replaced by the next mode.
%
% Outputs
%   freq_Hz, damping : tracked oscillatory branches. damping is Re(s) [1/s].
%   A_hist           : first-order state matrices at each airspeed.
%   divergence_indicator : det(Kr+Sar)/det(Kr).
%   max_growth       : max Re(s) over the COMPLETE state-space spectrum.
%   dominant_eigenvalue : eigenvalue attaining max_growth at each speed.
%   real_eigenvalues : real/overdamped eigenvalues retained separately,
%                      sorted by real value at each speed; unused entries NaN.
%
% If n_modes_keep is [] or Inf, all available oscillatory branches are
% allocated. The MAC threshold below prevents an unmatched branch from
% silently inheriting another mode's plotted row.

    Mr = Ms(free_dofs, free_dofs);
    Kr = Ks(free_dofs, free_dofs);
    if isempty(Cs)
        Cr = zeros(size(Mr));
    else
        Cr = Cs(free_dofs, free_dofs);
    end

    n = numel(free_dofs);
    nU = numel(U_vec);
    if isempty(n_modes_keep) || isinf(n_modes_keep)
        n_modes_keep = n;
    end

    freq_Hz = nan(n_modes_keep, nU);
    damping = nan(n_modes_keep, nU);
    A_hist = cell(1,nU);
    divergence_indicator = nan(1,nU);
    max_growth = nan(1,nU);
    dominant_eigenvalue = nan(1,nU);
    real_eigenvalues = nan(2*n, nU);

    det_Kr = det(Kr);
    if det_Kr == 0
        error('Reduced structural stiffness matrix is singular.');
    end

    I = eye(n);
    Z = zeros(n);
    imag_tol = 1e-6;
    mac_threshold = 0.25;
    previous_vectors = [];

    for k = 1:nU
        U = U_vec(k);
        [Da, Sa, ~] = aero_strip_matrices(node_y, dof_map, c_dist, ...
                                           ea_loc_frac, a_ea, AR, e_s, ...
                                           aero, U, rho);
        Dar = Da(free_dofs, free_dofs);
        Sar = Sa(free_dofs, free_dofs);
        Keff = Kr + Sar;
        Ceff = Cr + Dar;

        divergence_indicator(k) = real(det(Keff) / det_Kr);

        Minv = Mr \ I;
        A = [Z, I; -Minv*Keff, -Minv*Ceff];
        [V, s_matrix] = eig(A);
        s = diag(s_matrix);

        [max_growth(k), dominant_idx] = max(real(s));
        dominant_eigenvalue(k) = s(dominant_idx);

        real_idx = abs(imag(s)) <= imag_tol;
        real_s = sort(real(s(real_idx)), 'ascend');
        real_eigenvalues(1:numel(real_s),k) = real_s;

        osc_idx = imag(s) > imag_tol;
        osc_s = s(osc_idx);
        osc_V = V(:,osc_idx);
        [osc_s, sort_idx] = sort_by_frequency(osc_s);
        osc_V = osc_V(:,sort_idx);

        current_vectors = nan(2*n, n_modes_keep);
        current_s = nan(n_modes_keep,1);
        assigned_rows = false(n_modes_keep,1);
        assigned_cols = false(numel(osc_s),1);

        if k == 1
            n_initial = min(n_modes_keep, numel(osc_s));
            current_s(1:n_initial) = osc_s(1:n_initial);
            current_vectors(:,1:n_initial) = osc_V(:,1:n_initial);
            assigned_rows(1:n_initial) = true;
            assigned_cols(1:n_initial) = true;
        else
            % Greedy maximum-MAC assignment. Each branch can receive at most
            % one candidate, and each candidate can be used once.
            % Store explicit row indices. This avoids MATLAB implicit
            % expansion between row and column logical masks.
            previous_valid = find(~isnan(previous_vectors(1,:)));
            previous_valid = previous_valid(previous_valid >= 1 & ...
                                            previous_valid <= n_modes_keep);
            if ~isempty(previous_valid) && ~isempty(osc_V)
                mac = -Inf(n_modes_keep, numel(osc_s));
                old_V = previous_vectors(:,previous_valid);
                mac_valid = abs(old_V' * osc_V).^2;
                mac_valid = mac_valid ./ max(sum(abs(old_V).^2,1).' * ...
                                             sum(abs(osc_V).^2,1), eps);
                mac(previous_valid,:) = mac_valid;
                while true
                    available = intersect(previous_valid, find(~assigned_rows));
                    available_cols = find(~assigned_cols);
                    available = available(available >= 1 & ...
                                          available <= size(mac,1));
                    available_cols = available_cols(available_cols >= 1 & ...
                                                    available_cols <= size(mac,2));
                    if isempty(available) || isempty(available_cols)
                        break;
                    end
                    mac_work = mac(available,available_cols);
                    [best_mac, linear_idx] = max(mac_work(:));
                    if best_mac < mac_threshold
                        break;
                    end
                    [ir, ic] = ind2sub(size(mac_work), linear_idx);
                    row = available(ir);
                    col = available_cols(ic);
                    current_s(row) = osc_s(col);
                    current_vectors(:,row) = osc_V(:,col);
                    assigned_rows(row) = true;
                    assigned_cols(col) = true;
                    mac(row,:) = -Inf;
                    mac(:,col) = -Inf;
                end
            end

            % New branches occupy genuinely unused rows. Existing unmatched
            % rows remain NaN at this airspeed instead of being compressed.
            unused_rows = setdiff(find(~assigned_rows), previous_valid);
            unused_cols = find(~assigned_cols);
            n_new = min(numel(unused_rows), numel(unused_cols));
            for j = 1:n_new
                row = unused_rows(j);
                col = unused_cols(j);
                current_s(row) = osc_s(col);
                current_vectors(:,row) = osc_V(:,col);
            end
        end

        valid = ~isnan(current_s);
        freq_Hz(valid,k) = imag(current_s(valid))/(2*pi);
        damping(valid,k) = real(current_s(valid));
        previous_vectors = current_vectors;
        A_hist{k} = A;
    end
end

function [s_sorted, idx] = sort_by_frequency(s)
    [~, idx] = sort(imag(s), 'ascend');
    s_sorted = s(idx);
end
