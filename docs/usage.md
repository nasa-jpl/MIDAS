# **Using MIDAS**

This section provides an overview of how to use MIDAS to compute the scattering of electromagnetic waves by complex-shaped scatterers.

## **Input Files**

MIDAS requires a specific directory structure, primarily expecting an `inputs/` folder containing the following configuration files:

### **1. Main Simulation Configuration**
*   **`Simulation_data.dat`**: The primary input file. It defines:
    *   **EM Wave Parameters**: Wavelength or frequency range and units (MHz, GHz, THz).
    *   **Scatterer Geometry**: The type of scatterer (e.g., Simple Sphere, Cylinder, or Complex Shape read from a shape file).
    *   **Numerical Methods**: Selection of the solver (e.g., `CBFM-E`, `MLCBFM-E`, `MoM`, or `RGE`).
    *   **Transmitter/Receiver Setup**: Configuration for incident and scattering directions (Uniform, Spherical T-design, or Lebedev quadrature).
    *   **Discretization**: The cell size (`Sc`) used to voxelize the scatterer.

### **2. Geometry and Shape Files**
Depending on the `type_s` defined in the main configuration, MIDAS may require:
*   **`SimShapes.dat`**: A list of available shape simulations and their paths.
*   **Shape Files (`.dat`)**: Detailed lattice positions for complex particles.
*   **`Cells.dat`**: An alternative input for predefined cell coordinates and dielectric indices.

### **3. Dielectric Properties**
Depending on the `dielcomp_option`, the code reads refractive index data from:
*   **`dielcomposition.dat`**: A file mapping refractive indices to specific cells.
*   **`dielectric_table.txt`**: A table of dielectric constants, useful for frequency-dependent materials.
*   **Individual `.m` files**: Specific files for different dielectric materials used in the scatterer.

---

## **Running Simulations**

MIDAS is designed for high-performance computing and utilizes **MPI (Message Passing Interface)** for parallel processing.

### **Execution Command**
To run the simulation, use `mpirun` to distribute the workload across multiple CPU cores:

```bash
mpirun -np <num_cores> ./midas_executable < inputs/Simulation_data.dat
```
<num_cores>: The number of processors you wish to allocate.
Note: If using the CBFM method, the number of processors should generally be less than or equal to the total number of blocks (Nblocks) to ensure optimal performance Main_Scattering.f90.

### **Output Generation** 
Once executed, the code creates a simulation-specific output folder containing:

*   **S_files/:** Scattering matrices (S-matrices) InOut_ReadWrite_Functions.f90 +1.
*   **Ei_files/ and Es_files/:** Incident and scattered electric field data Main_Scattering.f90.
*   **Sol_files/:** Internal fields and system matrices (if save_Eint or save_Zc are enabled) Main_Scattering.f90.
*   **Analysis/:** Log files and memory allocation tracking (TrackAllocate_j<rank>.dat) InOut_ReadWrite_Functions.f90.

### **Matlab Tools**