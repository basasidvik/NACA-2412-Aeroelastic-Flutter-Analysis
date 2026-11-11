function [t, tip_displacement, eta, Phi, modal_frequency_Hz, metrics] = ...
    modal_gust_response(Ms, Ks, Cs, free_dofs, node_y, dof_map, ...
                        c_dist, ea_loc_frac, a_ea, AR, e_s, aero, ...
                        U, rho, t_span, gust_velocity, n_modes, z0)
% MODAL_GUST_RESPONSE  Stiff time response using a truncated FEM modal basis.
%
% The full FEM matrices are used to calculate the modes. The aerodynamic
% matrices and gust forces are then projected onto the first n_modes
% structural modes, and the reduced time-domain system is integrated with
% ode15s. This is intended for convergence studies such as 4, 8, and 12
% retained modes.

    if nargin < 18 || isempty(z0)
        z0 = zeros(2*n_modes,1);
    end
    if n_modes < 1 || n_modes > numel(free_dofs)
        error('n_modes must be between 1 and the number of free DOFs.');
    end

    Mr = Ms(free_dofs,free_dofs);
    Kr = Ks(free_dofs,free_dofs);
    if isempty(Cs)
        Cr = zeros(size(Mr));
    else
        Cr = Cs(free_dofs,free_dofs);
    end

    % Generalized structural eigenproblem. The full FEM model supplies the
    % basis; modes are ordered from lowest to highest frequency.
    [V,D] = eig(Kr,Mr);
    lambda = real(diag(D));
    [lambda,order] = sort(lambda,'ascend');
    Phi = V(:,order(1:n_modes));
    lambda = lambda(1:n_modes);
    modal_frequency_Hz = sqrt(max(lambda,0))/(2*pi);

    [Da,Sa,~] = aero_strip_matrices(node_y,dof_map,c_dist, ...
                                     ea_loc_frac,a_ea,AR,e_s,aero,U,rho);
    Mmod = Phi'*Mr*Phi;
    Cmod = Phi'*(Cr + Da(free_dofs,free_dofs))*Phi;
    Kmod = Phi'*(Kr + Sa(free_dofs,free_dofs))*Phi;

    n_nodes = numel(node_y);
    q_inf = 0.5*rho*U^2;
    ClA3D = aero.Cl_alpha_2D / ...
            (1 + aero.Cl_alpha_2D/(pi*AR*e_s));
    dy = tributary_lengths(node_y);
    tip_w = find(free_dofs == dof_map(end,1),1,'first');

    force_full = @(tt) gust_force(tt,gust_velocity,node_y,dof_map, ...
                                  free_dofs,c_dist,ea_loc_frac,aero, ...
                                  q_inf,ClA3D,dy,n_nodes,U);
    force_modal = @(tt) Phi'*force_full(tt);

    n = n_modes;
    odefun = @(tt,z) [z(n+1:end); Mmod \ ...
        (force_modal(tt) - Cmod*z(n+1:end) - Kmod*z(1:n))];

    if numel(z0) ~= 2*n
        error('z0 must have length 2*n_modes.');
    end
    [t,Z] = ode15s(odefun,t_span,z0);
    eta = Z(:,1:n);
    tip_displacement = (Phi(tip_w,:)*eta.').';

    metrics.max_abs_tip = max(abs(tip_displacement));
    metrics.rms_tip = sqrt(mean(tip_displacement.^2));
    metrics.final_tip = tip_displacement(end);
    metrics.solver = 'ode15s';
end

function F = gust_force(t,gust_velocity,node_y,dof_map,free_dofs, ...
                        c_dist,ea_loc_frac,aero,q_inf,ClA3D,dy,n_nodes,U)
    Ffull = zeros(3*n_nodes,1);
    alpha_g = gust_velocity(t)/U;
    for i = 1:n_nodes
        ci = local_eval(c_dist,node_y(i));
        e_ac = (ea_loc_frac-aero.ac_loc)*ci;
        Lg = q_inf*ci*ClA3D*alpha_g*dy(i);
        g = dof_map(i,:);
        Ffull(g(1)) = Ffull(g(1))-Lg;
        Ffull(g(3)) = Ffull(g(3))+e_ac*Lg;
    end
    F = Ffull(free_dofs);
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
