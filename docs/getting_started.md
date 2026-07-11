# Getting Started with MIDAS

This section helps the user install and set up the MIDAS code and run his first simulation. The user can also start with the file [MIDAS_Overview_5-2026.pdf](../userguides/MIDAS_Overview_5-2026.pdf) under [`userguides/`](../userguides/) for a general overview of MIDAS, its installation and start-up examples.  

## Prerequisites
- FORTRAN compiler (Intel `ifort`/`ifx` recommended)
- MPI library (e.g., Intel MPI, OpenMPI)
- MATLAB (optional — needed only for input preparation and output visualization)

## Installation

1. Clone the repository:
   ```bash
   git clone https://github.jpl.nasa.gov/ifenni/MIDAS.git
   ```

2. Select (and update if needed) the appropriate makefile from the [`makefiles/`](../makefiles/) folder and place it in the [`src/`](../src/) folder.

3. Load the required modules. For example, on NASA Discover:
   ```bash
   module load comp/intel/2021.7.0
   module load mpi/impi/2021.7.0
   ```

4. Compile the FORTRAN code and build the `midas` executable, as shown below with the Makefile_discover example:
   ```bash
   cd MIDAS/src
   make -f Makefile_discover modules
   make -f Makefile_discover
   ```

5. Copy the executable to your run directory and submit the job:
   ```bash
   cp midas ../examples/Run_Sims_<name>/
   cd ../examples/Run_Sims_<name>
   ```

## Running the Examples

A set of `Run_Sims_*` folders is provided in the [`examples/`](../examples/) directory to help first-time users get started. Each example folder contains a ready-to-use input set.

The submitted job will use the input files in the following structure:

```
Run_Sims_<name>/
├── midas                    ← compiled executable
├── job_MIDAS_compute        ← job submission script
└── inputs/
    ├── Simulation_data.dat  ← main configuration file
    ├── shape.dat            ← particle geometry (if applicable)
    └── diel/
        └── <m_files>        ← refractive index files
```

More job files examples, adapted to the JPL Gattaca2, TACC LS6 and NCCS Discover servers, are available under [`tools/`](../tools/).
See the [Usage](usage.md) page for a detailed description of each input file.
