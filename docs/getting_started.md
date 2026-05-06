# Getting Started with MIDAS

This section helps you set up the MIDAS code and run your first simulation.

## Prerequisites
- FORTRAN compiler (e.g., gfortran)
- MPI library (e.g., OpenMPI)
- MATLAB (for input/output tools)

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.jpl.nasa.gov/ifenni/MIDAS.git
   
2. Select (and update if needed) adequate makefile from the “makefiles” folder and place it in the src folder

3. Load required modules. For example, for discover:
   ```bash
   module load comp/intel/2021.7.0
   module load mpi/impi/2021.7.0 

3. Compile the FORTRAN code and build the executable midas in the src folder
    ```bash
   cd MIDAS/src
   make -f Makefile_discover modules
   make -f Makefile_discover 
  
4. Create a 'Run_sims' folder and Copy the executable to the 'Run_Sims_*' folder and submit job 
   ```bash
   cp midas ../Run_Sims/
   cd ../Run_Sims
A set of Run_Sims_* folders is provided in the Examples directory to help first-time users get started and become familiar with the code. 
The submitted job will use the input files in the 'Run_Sims_*\inputs' and ''Run_Sims_*\inputs\diel' subfolders, therefore 'Run_Sims_*' (may be renamed as desired) and inputs folders must have the structure shown in the figure below. 
