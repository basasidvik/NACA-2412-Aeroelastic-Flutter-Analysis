function [t, Z, F_hist] = gust_response(Ms, Ks, Cs, free_dofs, node_y, dof_map, ...
                                        c_dist, ea_loc_frac, AR, e_s, aero, ...
                                        U, rho, t_span, gust_velocity, z0)
% GUST_RESPONSE  Forced time response to a prescribed vertical gust.
%
% The gust is represented as a small additional angle of attack,
% alpha_g = w_g/U. The resulting strip lift is assembled as a generalized
% force with downward-positive plunge and nose-up-positive twist:
%   Q_w     = -L_g,       Q_theta = e*L_g.
%
% gust_velocity is a function handle @(t) returning vertical gust velocity
% [m/s]. z0 is optional and has length 2*numel(free_dofs).

    Mr = Ms(free_dofs, free_dofs);
    Kr = Ks(free_dofs, free_dofs);
    if isempty(Cs)
        Cr = zeros(size(Mr));
    else
        Cr = Cs(free_dofs, free_dofs);
    end

    [Da, Sa, ~] = aero_strip_matrices(node_y, dof_map, c_dist, ...
                                      ea_loc_frac, 0, AR, e_s, aero, U, rho);
    Kr = Kr + Sa(free_dofs, free_dofs);
    Cr = Cr + Da(free_dofs, free_dofs);

    n = numel(free_dofs);
    n_nodes = numel(node_y);
    q_inf = 0.5*rho*U^2;
    ClA3D = aero.Cl_alpha_2D / (1 + aero.Cl_alpha_2D/(pi*AR*e_s));
    dy = tributary_lengths(node_y);

    if nargin < 15 || isempty(z0)
        z0 = zeros(2*n,1);
    end
    if numel(z0) ~= 2*n
        error('z0 must have length 2*numel(free_dofs).');
    end

    force_at_time = @(tt) gust_force(tt, gust_velocity, node_y, dof_map, ...
                                     free_dofs, c_dist, ea_loc_frac, ...
                                     aero, q_inf, ClA3D, dy, n_nodes, U);
    odefun = @(tt,z) [z(n+1:end); Mr \ (force_at_time(tt) - ...
                         Cr*z(n+1:end) - Kr*z(1:n))];

    % The full FEM system contains widely separated frequencies and is
    % stiff; ode15s avoids the impractical step sizes required by ode45.
    [t, Z] = ode15s(odefun, t_span, z0);
    F_hist = zeros(numel(t), n);
    for k = 1:numel(t)
        F_hist(k,:) = force_at_time(t(k)).';
    end
end

function F = gust_force(t, gust_velocity, node_y, dof_map, free_dofs, ...
                        c_dist, ea_loc_frac, aero, q_inf, ClA3D, dy, n_nodes, U)
    Ffull = zeros(3*n_nodes,1);
    alpha_g = gust_velocity(t) / U;
    for i = 1:n_nodes
        ci = local_eval(c_dist, node_y(i));
        ei = (ea_loc_frac - aero.ac_loc) * ci;
        Lg = q_inf * ci * ClA3D * alpha_g * dy(i);
        g = dof_map(i,:);
        Ffull(g(1)) = Ffull(g(1)) - Lg;
        Ffull(g(3)) = Ffull(g(3)) + ei*Lg;
    end
    F = Ffull(free_dofs);
end

function v = local_eval(prop, y)
    if isa(prop, 'function_handle')
        v = prop(y);
    else
        v = prop;
    end
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
