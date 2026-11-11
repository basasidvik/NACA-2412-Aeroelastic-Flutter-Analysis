function [Ms, Ks, node_y, dof_map, free_dofs] = beam_fem_matrices(L, n_el, EI, GJ, m, I_theta)
% BEAM_FEM_MATRICES  Coupled bending-torsion cantilever beam FEM model.
%
%   [Ms, Ks, node_y, dof_map, free_dofs] = beam_fem_matrices(L, n_el, EI, GJ, m, I_theta)
%
% Inputs
%   L       : span (root to tip) [m]
%   n_el    : number of beam elements (n_el+1 nodes)
%   EI      : bending stiffness [N*m^2]       (scalar = uniform, or function handle EI(y))
%   GJ      : torsional stiffness [N*m^2]     (scalar or function handle GJ(y))
%   m       : mass per unit length [kg/m]     (scalar or function handle m(y))
%   I_theta : mass moment of inertia per unit length about elastic axis [kg*m]
%             (scalar or function handle I_theta(y))
%
% Each node has 3 DOFs ordered [w, w', theta] (transverse deflection, slope, twist).
% Global DOF numbering: node i occupies DOFs 3*(i-1)+1 : 3*(i-1)+3.
% The root (node 1) is clamped: w=w'=theta=0 there are removed from the free set.
%
% Outputs
%   Ms, Ks     : global mass and stiffness matrices, FULL size 3*(n_el+1)
%                (unconstrained; use free_dofs to extract the reduced cantilever system)
%   node_y     : [n_el+1 x 1] spanwise coordinate of each node
%   dof_map    : [n_el+1 x 3] DOF index for each node's [w, w', theta]
%   free_dofs  : index vector of DOFs remaining after removing the clamped root DOFs

    n_nodes = n_el + 1;
    Le = L / n_el;
    node_y = linspace(0, L, n_nodes)';

    n_dof = 3 * n_nodes;
    Ms = zeros(n_dof);
    Ks = zeros(n_dof);

    dof_map = reshape(1:n_dof, 3, n_nodes)';   % node i -> dof_map(i,:) = [w_i, wp_i, th_i]

    % Evaluate (possibly y-dependent) properties at element midpoints
    getval = @(prop, yq) local_eval(prop, yq);

    for e = 1:n_el
        y1 = node_y(e);
        y2 = node_y(e+1);
        ymid = 0.5*(y1+y2);

        EIe = getval(EI, ymid);
        GJe = getval(GJ, ymid);
        me  = getval(m, ymid);
        Ite = getval(I_theta, ymid);

        % --- Bending (Euler-Bernoulli) element matrices, DOFs [w1 w1' w2 w2'] ---
        Kb = (EIe/Le^3) * [ 12,      6*Le,   -12,      6*Le;
                             6*Le,  4*Le^2,  -6*Le,   2*Le^2;
                           -12,     -6*Le,    12,     -6*Le;
                             6*Le,  2*Le^2,  -6*Le,   4*Le^2];

        Mb = (me*Le/420) * [156,    22*Le,   54,    -13*Le;
                              22*Le, 4*Le^2,  13*Le, -3*Le^2;
                              54,    13*Le,   156,   -22*Le;
                            -13*Le, -3*Le^2, -22*Le,  4*Le^2];

        % --- Torsion (rod-type) element matrices, DOFs [theta1 theta2] ---
        Kt = (GJe/Le) * [ 1, -1; -1, 1];
        Mt = (Ite*Le/6) * [ 2, 1; 1, 2];

        % --- Assemble into 6x6 element matrix, order [w1 w1' th1 w2 w2' th2] ---
        Ke = zeros(6);
        Me = zeros(6);

        bend_idx = [1 2 4 5];   % w1, w1', w2, w2'
        tors_idx = [3 6];       % th1, th2

        Ke(bend_idx, bend_idx) = Kb;
        Me(bend_idx, bend_idx) = Mb;
        Ke(tors_idx, tors_idx) = Kt;
        Me(tors_idx, tors_idx) = Mt;

        % --- Scatter into global matrices ---
        g = [dof_map(e,:), dof_map(e+1,:)];   % 6 global dof indices for this element
        Ks(g,g) = Ks(g,g) + Ke;
        Ms(g,g) = Ms(g,g) + Me;
    end

    % Remove clamped root DOFs (node 1: w, w', theta)
    root_dofs = dof_map(1,:);
    free_dofs = setdiff(1:n_dof, root_dofs);

end

function v = local_eval(prop, y)
    if isa(prop, 'function_handle')
        v = prop(y);
    else
        v = prop;
    end
end
