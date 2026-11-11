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

# NACA 2412 Aeroelastic Flutter Analysis — Version 2

MATLAB framework for qualitative aeroelastic analysis of a finite-span NACA 2412 wing using beam finite elements, quasi-steady strip theory, and an exploratory unsteady-aerodynamics extension.

## Features

- Euler–Bernoulli bending and Saint-Venant torsion FEM model.
- Cantilever wing boundary conditions.
- Finite-span lift-curve correction using Prandtl lifting-line theory.
- Quasi-steady aerodynamic stiffness and damping matrices.
- MAC-based oscillatory-mode tracking.
- Complete-spectrum growth-rate diagnostics.
- Real and overdamped eigenvalue handling.
- Static trim from NACA 2412 camber effects.
- Determinant and singular-value divergence diagnostics.
- Initial-condition and harmonic gust response.
- Exploratory Theodorsen circulatory extension.
- Low-order added-mass approximation.
- Jones rational approximation of Theodorsen’s function.
- Sears-based harmonic gust response.
- V–f, V–g, real-part, divergence, and displacement-comparison plots.

## Main analysis branches

### Quasi-steady branch

Run:
run_flutter_analysis

This branch provides structural modal analysis, quasi-steady aeroelastic stability, static trim, divergence diagnostics, and forced gust response.

### Theodorsen extension

Run:

run_theodorsen_analysis

This branch uses frequency-dependent circulatory aerodynamics and a fixed-point p–k-style iteration. It also compares quasi-steady and Theodorsen/Sears harmonic displacement responses.

## Documentation

- aeroelastic_naca2412.tex — Version 2 main formulation and MATLAB code guide.
- aeroelastic_theodorsen_continuation.tex — Theodorsen extension, added-mass approximation, rational approximation, Sears gust model, and validation status.

## Validity and limitations

The structural FEM and quasi-steady strip-theory branches are intended for qualitative analysis, engineering sensitivity studies, and education. They are not validated production flutter-prediction tools.

The Theodorsen branch is exploratory. It includes frequency-dependent circulatory aerodynamics, a low-order added-mass approximation, a rational Theodorsen option, and a Sears-based gust response. It does not yet provide a fully validated nonlinear global aeroelastic spectrum or a complete finite-wing unsteady aerodynamic model.

Absolute flutter speeds and displacement amplitudes should not be treated as experimental predictions until the model is benchmarked against a classical typical-section solution, higher-fidelity aerodynamics, or measured NACA 2412 data
