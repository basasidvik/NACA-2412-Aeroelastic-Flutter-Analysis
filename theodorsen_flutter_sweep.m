function [U_vec, freq_Hz, growth_rate, converged, eig_hist, ...
          max_growth, dominant_eigenvalue, real_eigenvalues, ...
          spectrum_real_parts] = ...
    theodorsen_flutter_sweep(Ms, Ks, free_dofs, node_y, dof_map, ...
                             c_dist, ea_loc_frac, a_ea, AR, e_s, aero, ...
                             rho, U_vec, n_modes_keep, Cs, model)
% THEODORSEN_FLUTTER_SWEEP  Frequency-domain Theodorsen p-k-style sweep.
%
% At every airspeed, each mode is solved by fixed-point iteration:
%   1. evaluate C(k) using the current trial eigenvalue;
%   2. assemble complex frequency-dependent aerodynamic matrices;
%   3. solve the resulting state-space eigenproblem;
%   4. select the eigenvalue using displacement-eigenvector MAC, with a
%      nearest-eigenvalue fallback.
%
% This is a practical Theodorsen extension with a low-order added-mass
% approximation, not a full unsteady-aerodynamics production solver. See
% the continuation LaTeX file for the included assumptions.

    Mr = Ms(free_dofs, free_dofs);
    Kr = Ks(free_dofs, free_dofs);
    if isempty(Cs)
        Cr = zeros(size(Mr));
    else
        Cr = Cs(free_dofs, free_dofs);
    end
    if nargin < 16 || isempty(model)
        model = 'jones';
    end

    n = numel(free_dofs);
    nU = numel(U_vec);
    if isempty(n_modes_keep) || isinf(n_modes_keep)
        n_modes_keep = n;
    end

    freq_Hz = nan(n_modes_keep,nU);
    growth_rate = nan(n_modes_keep,nU);
    converged = false(n_modes_keep,nU);
    eig_hist = cell(1,nU);
    max_growth = nan(1,nU);
    dominant_eigenvalue = nan(1,nU);
    real_eigenvalues = nan(2*n,nU);
    % Store the 2n candidate real parts nearest the stability boundary for
    % a readable diagnostic plot; the complete maximum is still retained.
    spectrum_real_parts = nan(2*n,nU);
    previous_modes = nan(n_modes_keep,1);
    I = eye(n);
    Z = zeros(n);
    imag_tol = 1e-6;
    max_iter = 40;
    relaxation = 0.6;
    tolerance = 1e-7;
    mac_threshold = 0.25;
    previous_vectors = nan(n,n_modes_keep);

    for ku = 1:nU
        U = U_vec(ku);

        % Quasi-steady eigenvalues provide robust seeds at the first speed
        % and a fallback for any branch that has become unavailable.
        [Da0, Sa0] = aero_strip_matrices(node_y, dof_map, c_dist, ...
                                         ea_loc_frac, a_ea, AR, e_s, ...
                                         aero, U, rho);
        A0 = state_matrix(Mr, Kr + Sa0(free_dofs,free_dofs), ...
                          Cr + Da0(free_dofs,free_dofs), I, Z);
        seed_all = eig(A0);
        seed_all = seed_all(imag(seed_all) > imag_tol);
        seed_all = seed_all(sort_by_imag(seed_all));

        if ku == 1
            seeds = nan(n_modes_keep,1);
            nseed = min(n_modes_keep,numel(seed_all));
            seeds(1:nseed) = seed_all(1:nseed);
        else
            seeds = previous_modes;
            missing = isnan(seeds);
            nmissing = sum(missing);
            nfill = min(nmissing,numel(seed_all));
            if nfill > 0
                missing_idx = find(missing, nfill, 'first');
                seeds(missing_idx) = seed_all(1:nfill);
            end
        end

        solved = nan(n_modes_keep,1);
        solved_vectors = nan(n,n_modes_keep);
        spectrum_candidates = [];
        for im = 1:n_modes_keep
            if isnan(seeds(im))
                continue;
            end
            s_trial = seeds(im);
            selected_vector = nan(n,1);
            for it = 1:max_iter
                [Da, Sa, Ma] = unsteady_aero_strip_matrices( ...
                    node_y, dof_map, c_dist, ea_loc_frac, a_ea, ...
                    AR, e_s, aero, U, rho, s_trial, model);
                A = state_matrix(Mr + Ma(free_dofs,free_dofs), ...
                                 Kr + Sa(free_dofs,free_dofs), ...
                                 Cr + Da(free_dofs,free_dofs), I, Z);
                [V,eig_matrix] = eig(A);
                candidates = diag(eig_matrix);
                candidate_idx = find(imag(candidates) > imag_tol);
                if isempty(candidate_idx)
                    candidate_idx = (1:numel(candidates)).';
                end
                candidate_values = candidates(candidate_idx);
                candidate_vectors = V(1:n,candidate_idx);
                if ~isnan(previous_vectors(1,im))
                    mac = abs(previous_vectors(:,im)'*candidate_vectors).^2;
                    mac = mac ./ max(sum(abs(previous_vectors(:,im)).^2) * ...
                                     sum(abs(candidate_vectors).^2,1), eps);
                    [best_mac, idx] = max(mac);
                    if best_mac < mac_threshold
                        [~, idx] = min(abs(candidate_values-s_trial));
                    end
                else
                    [~, idx] = min(abs(candidate_values-s_trial));
                end
                s_new = candidate_values(idx);
                selected_vector = candidate_vectors(:,idx);
                if abs(s_new-s_trial) <= tolerance*(1+abs(s_trial))
                    s_trial = s_new;
                    converged(im,ku) = true;
                    break;
                end
                s_trial = relaxation*s_new + (1-relaxation)*s_trial;
            end
            solved(im) = s_trial;
            solved_vectors(:,im) = selected_vector;
            spectrum_candidates = [spectrum_candidates; eig(A)]; %#ok<AGROW>
        end

        imag_solved = imag(solved);
        imag_solved(isnan(imag_solved)) = Inf;
        [~, order] = sort(imag_solved, 'ascend');
        solved = solved(order);
        solved_vectors = solved_vectors(:,order);
        converged(:,ku) = converged(order,ku);
        eig_hist{ku} = solved;
        freq_Hz(:,ku) = imag(solved)/(2*pi);
        growth_rate(:,ku) = real(solved);
        if ~isempty(spectrum_candidates)
            candidate_growth = real(spectrum_candidates);
            max_growth(ku) = max(candidate_growth);
            positive_idx = find(imag(spectrum_candidates) >= 0);
            if ~isempty(positive_idx)
                [positive_growth, positive_local] = ...
                    max(candidate_growth(positive_idx));
            else
                positive_growth = -Inf;
                positive_local = [];
            end
            if ~isempty(positive_idx) && positive_growth >= max_growth(ku)-1e-8
                dominant_idx = positive_idx(positive_local);
            else
                [~,dominant_idx] = max(candidate_growth);
            end
            dominant_eigenvalue(ku) = spectrum_candidates(dominant_idx);
            real_s = real(spectrum_candidates(abs(imag(spectrum_candidates)) <= imag_tol));
            real_s = sort(real_s,'ascend');
            real_eigenvalues(1:min(numel(real_s),2*n),ku) = ...
                real_s(1:min(numel(real_s),2*n));
            [~,near_zero_order] = sort(abs(candidate_growth),'ascend');
            n_store = min(numel(near_zero_order),size(spectrum_real_parts,1));
            spectrum_real_parts(1:n_store,ku) = ...
                candidate_growth(near_zero_order(1:n_store));
        end
        previous_modes = solved;
        previous_vectors = solved_vectors;
    end
end

function A = state_matrix(M, K, C, I, Z)
    A = [Z, I; -(M\K), -(M\C)];
end

function idx = sort_by_imag(s)
    [~,idx] = sort(imag(s),'ascend');
end
