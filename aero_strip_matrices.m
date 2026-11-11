function [Da, Sa, Qstatic] = aero_strip_matrices(node_y, dof_map, c_dist, ea_loc_frac, ...
                                                   a_ea, AR, e_s, aero, U, rho)
% AERO_STRIP_MATRICES  Quasi-steady strip-theory aerodynamic damping/stiffness
% matrices, lumped at FEM nodes ("divide the wing into sections").
%
%   [Da, Sa, Qstatic] = aero_strip_matrices(node_y, dof_map, c_dist, ea_loc_frac, ...
%                                            a_ea, AR, e_s, aero, U, rho)
%
% Inputs
%   node_y      : [n_nodes x 1] spanwise station of each node (from beam_fem_matrices)
%   dof_map     : [n_nodes x 3] global DOF index per node, [w, w', theta]
%   c_dist      : chord length [m], scalar or function handle c(y)
%   ea_loc_frac : elastic axis location, fraction of chord from leading edge (e.g. 0.33)
%   a_ea        : elastic axis location, fraction from MID-chord (a_ea = 2*ea_loc_frac - 1).
%                 Accepted for forward compatibility (e.g. an unsteady/Theodorsen
%                 extension) but NOT used by the baseline model below -- see note.
%   AR          : wing aspect ratio (b^2/S), for the Prandtl finite-span correction
%   e_s         : Oswald span efficiency factor (~0.8-0.95 typical)
%   aero        : struct from naca2412_data() (Cl_alpha_2D, Cm_ac, ac_loc, ...)
%   U           : freestream airspeed [m/s]
%   rho         : air density [kg/m^3]
%
% Outputs
%   Da, Sa    : aerodynamic damping and stiffness matrices, FULL size 3*n_nodes,
%               to be added to the structural matrices as
%                   Ms*qdd + (Cs + Da)*qd + (Ks + Sa)*q = 0
%               (reduce with free_dofs exactly as Ms, Ks are reduced before use)
%   Qstatic   : generalized static/trim aerodynamic force vector (from alpha0
%               and Cm_ac), FULL size 3*n_nodes -- not part of the linear
%               stability matrices, used for static trim/deflection checks.
%
% Modeling notes
%   - Each node is treated as one aerodynamic "strip" of width equal to the
%     trapezoidal tributary length around that node (half of each adjacent
%     element). This directly implements "divide the span into sections and
%     sum the effective lift force", per node, rather than a consistent
%     shape-function-weighted integral.
%   - The slope DOF w' carries no direct aerodynamic generalized force in this
%     lumped formulation (a consistent-FEM treatment would add one); this is a
%     deliberate simplification, noted in the accompanying LaTeX document.
%   - Strips are NOT cross-coupled aerodynamically beyond the single global
%     finite-span correction to Cl_alpha (Eq. for C_Lalpha^3D in the document).
%     A spanwise-varying induced-angle correction (vortex-lattice) is the
%     natural next refinement.
%   - BASELINE MODEL SCOPE: only the angle-of-attack / plunge-rate circulatory
%     lift (alpha_eff = theta + wdot/U) is implemented. A full quasi-steady
%     strip model also includes a pitch-RATE-dependent circulatory lift term
%     (the classic "(1/2 - a) c/U * thetadot" contribution of Theodorsen's
%     quasi-steady limit) and its own, DISTINCT moment arm about the elastic
%     axis -- that moment arm is NOT simply e(y) (the AC-to-EA offset used
%     below for the angle-of-attack term), and getting it wrong is worse than
%     leaving it out. That term is intentionally omitted here rather than
%     included on unverified footing; add it back using a checked reference
%     (e.g. Theodorsen 1935, or Hodges & Pierce "Introduction to Structural
%     Dynamics and Aeroelasticity") if you need the extra pitch damping.
%   - PRACTICAL CONSEQUENCE: because this baseline has no built-in torsional
%     aerodynamic damping, the torsion-dominated mode can show very weak
%     positive damping (mild instability) even at low airspeed when
%     structural damping (Cs) is left at zero -- this is an expected result
%     of the omission above, not an assembly bug (confirmed by an independent
%     re-implementation during development). Supplying even small structural
%     (Rayleigh) damping removes this artifact; run_flutter_analysis.m does
%     this by default.

    n_nodes = size(dof_map,1);
    n_dof   = 3*n_nodes;

    Da = zeros(n_dof);
    Sa = zeros(n_dof);
    Qstatic = zeros(n_dof,1);

    q_inf = 0.5*rho*U^2;

    % --- Finite-span corrected lift-curve slope (Prandtl lifting-line) ---
    ClA2D = aero.Cl_alpha_2D;
    ClA3D = ClA2D / (1 + ClA2D/(pi*AR*e_s));

    % --- Trapezoidal tributary width (strip length) at each node ---
    dy = tributary_lengths(node_y);

    for i = 1:n_nodes
        yi  = node_y(i);
        ci  = local_eval(c_dist, yi);
        ei  = (ea_loc_frac - aero.ac_loc) * ci;   % AC-to-EA offset (positive: EA aft of AC)

        % Per-unit-span sensitivities of L'(y,t), scaled by strip width dy(i)
        % to give nodal (lumped) force/moment sensitivities. (No pitch-rate
        % circulatory term -- see "BASELINE MODEL SCOPE" above.)
        Ltheta     = q_inf * ci * ClA3D          * dy(i);  % dL/dtheta
        Lwdot      = q_inf * ci * ClA3D / U      * dy(i);  % dL/d(wdot)

        Mtheta     = ei * Ltheta;
        Mwdot      = ei * Lwdot;

        % Local 3x3 blocks, DOF order [w, w', theta].
        %
        % Sign convention: w is positive DOWNWARD (document Sec. 1), while the
        % physical lift L is positive UPWARD (standard aero convention). The
        % generalized force conjugate to w is therefore -L (upward lift resists
        % increasing downward deflection: classic typical-section convention,
        % m*hddot + K*h = -L). The moment M_ea is already expressed about the
        % elastic axis consistent with theta positive nose-up, so it enters the
        % theta-row with no sign flip: I_theta*thetaddot + K*theta = +M_ea.
        Sa_local = zeros(3);
        Sa_local(1,3) = +Ltheta;     % -d(-L)/dtheta      into row w
        Sa_local(3,3) = -Mtheta;     % -dM/dtheta         into row theta

        Da_local = zeros(3);
        Da_local(1,1) = +Lwdot;      % -d(-L)/d(wdot)      into row w
        Da_local(3,1) = -Mwdot;      % -dM/d(wdot)         into row theta

        g = dof_map(i,:);
        Sa(g,g) = Sa(g,g) + Sa_local;
        Da(g,g) = Da(g,g) + Da_local;

        % Static camber-moment contribution (not part of Da/Sa)
        if isfield(aero,'alpha0')
            alpha0 = aero.alpha0;
        else
            alpha0 = 0;
        end
        Lstatic = q_inf * ci * ClA3D * (-alpha0) * dy(i);
        Mstatic = q_inf * ci^2 * aero.Cm_ac * dy(i);
        Qstatic(g(1)) = Qstatic(g(1)) - Lstatic;
        Qstatic(g(3)) = Qstatic(g(3)) + ei*Lstatic;
        Qstatic(g(3)) = Qstatic(g(3)) + Mstatic;
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
