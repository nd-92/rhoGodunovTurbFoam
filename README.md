# buiTurbFoam-v2212
Density-based compressible transient solver using the flux difference splitting scheme of Bui and Runge Kutta 4-stage time integration.  Primarily designed for LES, and with optional absorption-based acoustic damping.

The HEqn version of the solver takes enthalpy-based thermo packages.  The default version uses internal energy.  This solver originally started life as dbnsTurbFoam from foam-extend but has been modified to use a different flux scheme.
