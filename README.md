# NACA-2412-airfoil-vibration-response
# NACA 2412 Aeroelastic Flutter Analysis

MATLAB implementation of a finite-span NACA 2412 wing aeroelastic model using Euler–Bernoulli beam finite elements and quasi-steady strip theory.

The model represents a cantilever wing with bending, slope, and torsional degrees of freedom at each spanwise node. It computes structural mass and stiffness matrices, aerodynamic stiffness and damping matrices, modal frequencies, eigenvalue growth rates, static-divergence indicators, and forced gust responses.

## Main capabilities

- Euler–Bernoulli bending and Saint-Venant torsion FEM model.
- Finite-span lift-curve correction using Prandtl lifting-line theory.
- Quasi-steady aerodynamic strip loads.
- Complete state-space eigenvalue analysis over an airspeed sweep.
- MAC-based oscillatory mode tracking.
- Separate retention and plotting of real/overdamped eigenvalues.
- Complete-spectrum instability detection.
- Static-divergence determinant monitoring.
- V–f, V–g, real-eigenvalue, and divergence plots.
- Time-domain response to an initial disturbance.
- Forced response to a prescribed vertical gust.

## Files

### `run_flutter_analysis.m`

Top-level analysis script. Defines the wing geometry, material and structural properties, aerodynamic parameters, damping model, discretization, and airspeed sweep. It runs the aeroelastic analysis, reports flutter and divergence results, and generates plots.

### `beam_fem_matrices.m`

Builds the structural finite-element model. It assembles the bending and torsional element mass and stiffness matrices and removes the clamped-root degrees of freedom.

### `aero_strip_matrices.m`

Constructs the aerodynamic damping and stiffness matrices using node-based strip theory. It includes the finite-span lift-curve correction and the prescribed downward-positive plunge/nose-up-positive twist sign convention.

### `flutter_sweep.m`

Performs the airspeed sweep and state-space eigenvalue analysis. It:

- computes the complete eigenvalue spectrum;
- tracks oscillatory modes using eigenvector MAC correlation;
- leaves unavailable oscillatory branches as `NaN`;
- retains real and overdamped eigenvalues separately;
- reports the dominant eigenvalue and maximum growth rate;
- evaluates the static-divergence determinant indicator.

### `gust_response.m`

Computes the time-domain response to a prescribed vertical gust velocity. The gust is converted into additional aerodynamic lift and pitching-moment forces and integrated using `ode45`.

### `naca2412_data.m`

Provides representative low-Mach NACA 2412 aerodynamic properties, including lift-curve slope, zero-lift angle, aerodynamic-center location, and pitching-moment coefficient.

## Model limitations

This is a first-principles teaching and research template rather than a validated production aeroelasticity tool. The baseline model does not currently include:

- pitch-rate-dependent circulatory lift;
- Theodorsen unsteady aerodynamic effects;
- spanwise-varying induced-angle calculations;
- a vortex-lattice or CFD aerodynamic model;
- consistent shape-function integration of aerodynamic loads;
- experimentally calibrated structural damping;
- gust alleviation or atmospheric turbulence spectra.

Consequently, calculated flutter and divergence speeds should be treated as approximate and used for qualitative analysis unless the model is further calibrated against experimental or higher-fidelity aerodynamic data.

## Typical result

For the default parameters, the analysis reports a dynamic instability near 23.5 m/s and no static-divergence crossing within the swept range. Results depend strongly on the structural stiffness, mass distribution, damping, elastic-axis location, and aerodynamic assumptions
