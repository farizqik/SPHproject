EWCSPH Python Viewers — Fast v5: Conditional Analytical Potential Flow

Run in Spyder using F5 on ewcsph_single_viewer.py, or in Windows CMD:
"C:\Users\Farizqi\AppData\Local\spyder-6\envs\spyder-runtime\python.exe" "C:\Users\Farizqi\SPHproject\EWCSPH_Python_Viewers\ewcsph_single_viewer.py"

NEW IN SINGLE VIEWER:
- Checks SimulationConstants.csv: boundaryCondition=free-slip (also accepts freeslip/slip),
  or testScenario=potential for legacy metadata missing boundaryCondition.
- For free-slip data ONLY, shows two matched field panels: SPH and analytical
  potential flow around the circular cylinder (identical nondimensional colour scales).
- Adds Potential-Flow Validation tab with Cp around first fluid ring, tangential
  velocities at 90/270 degrees, RMS velocity error vs radius, and existing C++ L2.
- Overlays analytical reference on tangential and arbitrary-angle radial profiles.
- No-slip cases retain one-panel field and NO analytical potential-flow tab/overlays.
- Keeps all Fast v4 features: change folder, open another viewer, MP4 export, etc.
- Axis formatting disables additive +1 style offsets while retaining x 10^n scaling.

ANALYTICAL MODEL:
- u_r/Uinf = (1-a^2/r^2) cos(theta-flowAngle),
  u_theta/Uinf = -(1+a^2/r^2) sin(theta-flowAngle),
  Cp = 1 - (|u|/Uinf)^2. The steady, inviscid model has zero vorticity.
- Reference density uses the C++ barotropic EOS if valid B and gammaEOS
  constants are available; otherwise returns rho/rho0=1.
- It is a baseline ideal-potential-flow comparison, NOT an exact solution of
  viscous vortex shedding or a weakly-compressible transient flow.
- Cp validation is calculated at the first fluid ring (not incorrectly assumed
  to be the wall), and the exact curve is evaluated at that same radius.

TESTING:
- Python syntax compilation passed.
- Automated graphical smoke tests (synthetic C++-format CSV data) passed for
  both free-slip and no-slip scenarios on a virtual display.
- MATLAB/Windows and your full real datasets were not exercised.
