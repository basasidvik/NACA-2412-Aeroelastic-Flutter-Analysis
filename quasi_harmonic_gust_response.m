function [q_hat, F_hat, tip_amplitude] = quasi_harmonic_gust_response( ...
    Ms,Ks,Cs,free_dofs,node_y,dof_map,c_dist,ea_loc_frac,AR,e_s, ...
    aero,U,rho,gust_amplitude,gust_frequency)
% QUASI_HARMONIC_GUST_RESPONSE  Frequency-domain quasi-steady gust response.

    Mr = Ms(free_dofs,free_dofs);
    Kr = Ks(free_dofs,free_dofs);
    if isempty(Cs)
        Cr = zeros(size(Mr));
    else
        Cr = Cs(free_dofs,free_dofs);
    end
    [Da,Sa] = aero_strip_matrices(node_y,dof_map,c_dist, ...
                                  ea_loc_frac,0,AR,e_s,aero,U,rho);
    omega = 2*pi*gust_frequency;
    q_inf = 0.5*rho*U^2;
    ClA3D = aero.Cl_alpha_2D / ...
            (1 + aero.Cl_alpha_2D/(pi*AR*e_s));
    dy = tributary_lengths(node_y);
    Ffull = zeros(3*numel(node_y),1);
    for i = 1:numel(node_y)
        ci = local_eval(c_dist,node_y(i));
        e_ac = (ea_loc_frac-aero.ac_loc)*ci;
        Lg = q_inf*ci*ClA3D*(gust_amplitude/U)*dy(i);
        g = dof_map(i,:);
        Ffull(g(1)) = Ffull(g(1))-Lg;
        Ffull(g(3)) = Ffull(g(3))+e_ac*Lg;
    end
    F_hat = Ffull(free_dofs);
    q_hat = (-omega^2*Mr + 1i*omega*(Cr+Da(free_dofs,free_dofs)) + ...
             Kr+Sa(free_dofs,free_dofs)) \ F_hat;
    tip_w = find(free_dofs == dof_map(end,1),1,'first');
    tip_amplitude = abs(q_hat(tip_w));
end

function dy = tributary_lengths(node_y)
    n = numel(node_y); dy = zeros(n,1);
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
    if isa(prop,'function_handle'); v = prop(y); else; v = prop; end
end
