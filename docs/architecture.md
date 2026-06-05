# MIDAS Architecture

MIDAS uses a **domain-decomposition method** for solving electromagnetic scattering from arbitrarily-shaped objects. The project is structured as follows:

## Code Components
- **FORTRAN Core**: Full-wave 3D solver, MPI parallelized
- **MATLAB Tools**:
  - Input preparation (shape, material properties, simulation setup)
  - Output processing (scattered fields, matrices, efficiency factors)

## Key Source Modules

| File | Role |
|---|---|
| `Main_Scattering.f90` | Program entry point; orchestrates the full simulation loop over shapes and frequencies |
| `InOut_ReadWrite_Functions.f90` | All I/O: reads `Simulation_data.dat`, shape/cells files; writes geometry, S-matrix, and solution files |
| `DielComposition.f90` | Assigns refractive index to each cell per frequency; reads `.m` files, dielectric tables, or composition files |
| `Initialization.f90` | Data type definitions (`Scatterer`, `Cell`, `Dipole`, `CBFM_Block`) and global variable declarations |
| `Common_variables.f90` | Shared simulation state (wavelength, frequency, cell count, method flags, etc.) |
| `Discretization.f90` | Voxelizes the particle geometry into cubic cells |
| `Division_blocks.f90` | Partitions cells into CBFM blocks for domain decomposition |
| `Extend_blocks.f90` | Extends each block with neighboring cells to reduce edge effects |
| `MPI_distribution_blocks.f90` | Distributes CBFM blocks across MPI processes |
| `Compute_Electric_Fields.f90` | Top-level solver dispatcher; calls MoM, CBFM-E, MLCBFM-E, or RGE routines |

## Simulation Workflow

```
Simulation_data.dat
        │
        ▼
  Get_InputData          ← reads all parameters
        │
        ▼
  get_diel_values_lambdas ← loads refractive index values from m-files / table
        │
        ▼
  Discretization          ← voxelizes geometry (once, at highest frequency)
        │
        ▼
  Division_blocks  →  Extend_blocks  →  MPI_distribution_blocks
        │
        ▼
  ┌─────────────────────────────────┐
  │   Loop over frequencies         │
  │                                 │
  │   DielComposition               │ ← assigns m per cell at current frequency
  │   SetCellsParams                │ ← updates cell EM parameters
  │   Compute_Electric_Fields       │ ← solves for internal/scattered fields
  └─────────────────────────────────┘
        │
        ▼
  Write S_files / Q_files / Es_files / Ei_files / Sol_files
```

## Notes
- The FORTRAN code focuses on efficiency for large-scale 3D problems.
- MATLAB scripts are optional but simplify workflow and visualization.
- Discretization is performed once at the shortest internal wavelength (highest frequency / highest refractive index) so the mesh is valid across all requested frequencies.
