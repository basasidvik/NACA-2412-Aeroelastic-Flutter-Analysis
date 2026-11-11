function [q_hat, F_hat, tip_amplitude, k] = sears_gust_response( ...
    Ms, Ks, Cs, free_dofs, node_y, dof_map, c_dist, ea_loc_frac, ...
    a_ea, AR, e_s, aero, U, rho, gust_amplitude, gust_frequency, model)
% SEARS_GUST_RESPONSE  Harmonic gust response using Sears gust transfer.
%
% The input gust is a harmonic vertical velocity with amplitude
% gust_amplitude [m/s] and frequency gust_frequency [Hz]. The output q_hat
% is the complex displacement amplitude of the free structural DOFs.

    if nargin < 17 || isempty(model)
        model = 'jones';
    end
    if U <= 0 || gust_frequency <= 0
        error('U and gust_frequency must be positive.');
    end

    Mr = Ms(free_dofs,free_dofs);
    Kr = Ks(free_dofs,free_dofs);
    if isempty(Cs)
        Cr = zeros(size(Mr));
    else
        Cr = Cs(free_dofs,free_dofs);
    end

    omega = 2*pi*gust_frequency;
    s = 1i*omega;
    [Da,Sa,Ma] = unsteady_aero_strip_matrices( ...
        node_y,dof_map,c_dist,ea_loc_frac,a_ea,AR,e_s,aero,U,rho,s,model);
    M = Mr + Ma(free_dofs,free_dofs);
    C = Cr + Da(free_dofs,free_dofs);
    K = Kr + Sa(free_dofs,free_dofs);

    q_inf = 0.5*rho*U^2;
    dy = tributary_lengths(node_y);
    n_nodes = numel(node_y);
    k = zeros(n_nodes,1);
    Ffull = complex(zeros(3*n_nodes,1));
    ClA3D = aero.Cl_alpha_2D / ...
            (1 + aero.Cl_alpha_2D/(pi*AR*e_s));

    for i = 1:n_nodes
        ci = local_eval(c_dist,node_y(i));
        e_ac = (ea_loc_frac-aero.ac_loc)*ci;
        k(i) = omega*ci/(2*U); %#ok<AGROW>
        Sk = sears_function(k(i));
        Lg = q_inf*ci*ClA3D*(gust_amplitude/U)*Sk*dy(i);
        g = dof_map(i,:);
        Ffull(g(1)) = Ffull(g(1)) - Lg;
        Ffull(g(3)) = Ffull(g(3)) + e_ac*Lg;
    end

    F_hat = Ffull(free_dofs);
    q_hat = (-omega^2*M + 1i*omega*C + K) \ F_hat;
    tip_w = find(free_dofs == dof_map(end,1),1,'first');
    tip_amplitude = abs(q_hat(tip_w));
end

function dy = tributary_lengths(node_y)
    n = numel(node_y);
    dy = zeros(n,1);
    for i = 1:n
        if i == 1
            dy(i) = 0.5*(node_y(2)-node_y(1));
        elseif i == n
            dy(i) = 0.5*(node_y(n)-node_y(n-1));
        else
            dy(i) = 0.5*(node_y(i+1)-node_y(i-1));
        end
    end
end

function v = local_eval(prop,y)
    if isa(prop,'function_handle')
        v = prop(y);
    else
        v = prop;
    end
end
