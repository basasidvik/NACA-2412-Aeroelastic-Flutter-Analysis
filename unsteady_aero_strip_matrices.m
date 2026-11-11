function [Da, Sa, Ma, Ck, k] = unsteady_aero_strip_matrices(node_y, dof_map, ...
                                      c_dist, ea_loc_frac, a_ea, AR, e_s, ...
                                      aero, U, rho, s, model)
% UNSTEADY_AERO_STRIP_MATRICES  Frequency-dependent Theodorsen strip model.
%
% The matrices are evaluated for a trial harmonic eigenvalue s. The
% circulatory lift uses
%
% alpha_c = theta + s*w/U + (1/2-a_ea)*c*s*theta/U,
%
% multiplied by Theodorsen's C(k), where k=abs(imag(s))*c/(2U). The result
% is returned as complex, frequency-dependent damping/stiffness matrices.
%
% model is 'exact' (Hankel functions) or 'jones' (low-order rational
% approximation). Ma is a low-order non-circulatory added-mass approximation.

    if nargin < 12 || isempty(model)
        model = 'exact';
    end

    if U <= 0
        error('U must be positive.');
    end
    if ~isscalar(s)
        error('s must be a scalar trial eigenvalue.');
    end

    n_nodes = size(dof_map,1);
    n_dof = 3*n_nodes;
    Da = complex(zeros(n_dof));
    Sa = complex(zeros(n_dof));
    Ma = complex(zeros(n_dof));
    Ck = complex(zeros(n_nodes,1));
    k = zeros(n_nodes,1);

    q_inf = 0.5*rho*U^2;
    ClA3D = aero.Cl_alpha_2D / ...
            (1 + aero.Cl_alpha_2D/(pi*AR*e_s));
    dy = tributary_lengths(node_y);

    for i = 1:n_nodes
        ci = local_eval(c_dist, node_y(i));
        e_ac = (ea_loc_frac - aero.ac_loc) * ci;
        pitch_rate_arm = (0.5 - a_ea) * ci;
        k(i) = abs(imag(s))*ci/(2*U);
        if strcmpi(model,'jones')
            Ck(i) = theodorsen_rational_function(k(i));
        else
            Ck(i) = theodorsen_function(k(i));
        end

        Ltheta = q_inf*ci*ClA3D*Ck(i)*dy(i);
        Lwdot = q_inf*ci*ClA3D*Ck(i)/U*dy(i);
        Lthetadot = q_inf*ci*ClA3D*Ck(i)*pitch_rate_arm/U*dy(i);

        % The circulatory resultant is applied at the aerodynamic center.
        % Thus the same AC-to-EA arm is used for the lift-generated moment.
        Mtheta = e_ac*Ltheta;
        Mwdot = e_ac*Lwdot;
        Mthetadot = e_ac*Lthetadot;

        Sa_local = complex(zeros(3));
        Da_local = complex(zeros(3));
        Sa_local(1,3) = Ltheta;
        Sa_local(3,3) = -Mtheta;
        Da_local(1,1) = Lwdot;
        Da_local(1,3) = Lthetadot;
        Da_local(3,1) = -Mwdot;
        Da_local(3,3) = -Mthetadot;

        % Low-order non-circulatory added-mass approximation for a flat
        % plate strip. It is retained as an effective mass matrix.
        mu = pi*rho*ci^2/4*dy(i);
        Ma_local = complex(zeros(3));
        Ma_local(1,1) = mu;
        Ma_local(3,1) = -e_ac*mu;
        Ma_local(3,3) = -e_ac*mu*ci/2;
        Da_local(1,3) = Da_local(1,3) + mu*U;
        Da_local(3,3) = Da_local(3,3) - e_ac*mu*U;

        g = dof_map(i,:);
        Sa(g,g) = Sa(g,g) + Sa_local;
        Da(g,g) = Da(g,g) + Da_local;
        Ma(g,g) = Ma(g,g) + Ma_local;
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

function v = local_eval(prop, y)
    if isa(prop, 'function_handle')
        v = prop(y);
    else
        v = prop;
    end
end
