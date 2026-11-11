function aero = naca2412_data()
% NACA2412_DATA  Section aerodynamic properties for the NACA 2412 airfoil.
%
%   aero = naca2412_data() returns a struct with:
%     Cl_alpha_2D : 2D lift-curve slope [1/rad]
%     alpha0      : zero-lift angle of attack [rad]
%     Cm_ac       : moment coefficient about the aerodynamic center (quarter-chord)
%     ac_loc      : aerodynamic center location, fraction of chord from LE
%
% Values are representative low-Mach, pre-stall figures. Replace with
% experimental/panel-method data for a specific Reynolds number if available.

    aero.Cl_alpha_2D = 6.0;              % [1/rad], experimental (thin-airfoil theory gives 2*pi = 6.283)
    aero.alpha0      = -2.1*pi/180;      % [rad]
    aero.Cm_ac       = -0.047;           % about quarter-chord, ~constant pre-stall
    aero.ac_loc      = 0.25;             % fraction of chord from leading edge

end
