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
| `Initialization.f90` | Data type definitions (`Scatterer`, `Cell`, `Dipole`, `CBFM_Block`) and global variable declarations |
| `Discretization.f90` | Voxelizes the particle geometry into cubic cells |
| `Division_blocks.f90` | Partitions cells into CBFM blocks for domain decomposition |
| `Compute_EFields_CBFME.f90` | 1-Level CBFM-E solver: generates the characteristic basis functions per block, builds and solves the dense compressed system of linear equations (using ScaLAPACK), and reconstructs the internal total fields |
| `Compute_EFields_MoM.f90` | Full Method-of-Moments solver: assembles the distributed dense impedance matrix and solves it with ScaLAPACK `PZGESV` |
| `get_trans_Receiv.f90` | Builds the transmitter (incident) and receiver (scattering) direction sets for every integration type (`un/sd/lb/rf`) |
| `Incident_Field.f90` | Incident plane-wave field at the cells and at the receivers |
| `Compute_Scattered_Fields.f90` | Scattered electric fields at the receivers and `Es_files/` output (near-field path, `FFA=0`) |
| `Compute_Scattering_Matrices.f90` | Amplitude scattering (S) matrices saved to `S_files/` if `FFA=1` |
| `Compute_Scattering_Quantities.f90` | Orientation-averaged efficiency factors (Qext, Qsca, Qabs, Qbks) and asymmetry parameter `g`; `Q_files/`/`qtable` output |

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
