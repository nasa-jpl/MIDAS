# MIDAS Architecture

MIDAS uses a **domain-decomposition method** for solving electromagnetic scattering from arbitrarily-shaped objects. The project is structured as follows:

## Code Components
- **FORTRAN Core**: Full-wave 3D solver, MPI parallelized
- **MATLAB Tools**:
  - Input preparation (shape, material properties, simulation setup)
  - Output processing (scattered fields, matrices, efficiency factors)

## Workflow
1. Prepare simulation input in MATLAB
2. Generate input files for the FORTRAN solver
3. Run FORTRAN solver using MPI
4. Read and visualize outputs with MATLAB tools

## Notes
- The FORTRAN code focuses on efficiency for large-scale 3D problems.
- MATLAB scripts are optional but simplify workflow and visualization.
