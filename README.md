# NCAR Turbulence with Lagrangian Particles (NTLP)
NTLP is a MPI-based incompressible Navier-Stokes solver suitable for simulating
marine fog.  It supports three dimensional grids, periodic in the horizontal
with a second-order finite difference scheme in the vertical, along with two-way
coupled Lagrangian particles.

Intel ifort is required for compilation, as is either serial NetCDF or TecPlot.
Intel MPI is assumed to be available though NTLP should be compatible with any
MPI distribution installed with ifort.

Documentation is available in the [`documentation/`](documentation/README.md)
subdirectory.  New users should start with the following:

- [An NTLP Overview](documentation/explanation/ntlp-overview.md)
- [How to Setup and Run NTLP](documentation/tutorials/setup-and-run.md)
- [Post-processing Simulation Outputs](documentation/tutorials/post-processing-simulations.md)
