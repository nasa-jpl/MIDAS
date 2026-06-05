# **Using MIDAS**

This section provides an overview of how to use MIDAS to compute the scattering of electromagnetic waves by complex-shaped scatterers.

## **Input Files**

MIDAS requires a specific directory structure, primarily expecting an `inputs/` folder containing the following configuration files:

### **1. Main Simulation Configuration**
*   **`Simulation_data.dat`**: The primary input file. It defines:
    *   **EM Wave Parameters**: Wavelength or frequency range and units (MHz, GHz, THz) [4].
    *   **Scatterer Geometry**: The type of scatterer (e.g., Simple Sphere, Cylinder, or Complex Shape read from a shape file) [4].
    *   **Numerical Methods**: Selection of the solver (e.g., `CBFM-E`, `MLCBFM-E`, `MoM`, or `RGE`) [4].
    *   **Transmitter/Receiver Setup**: Configuration for incident and scattering directions (Uniform, Spherical T-design, or Lebedev quadrature) [4].
    *   **Discretization**: The cell size (`Sc`) used to voxelize the scatterer [4].

### **2. Geometry and Shape Files**
Depending on the `type_s` defined in the main configuration, MIDAS may require:
*   **`SimShapes.dat`**: A list of available shape simulations and their paths [5].
*   **Shape Files (`.dat`)**: Detailed lattice positions for complex particles, typically read via the `Read_ShapeFile` routine [4].
*   **`Cells.dat`**: An alternative input for predefined cell coordinates and dielectric indices [4].

### **3. Dielectric Properties**
Depending on the `dielcomp_option`, the code reads refractive index data from:
*   **`dielcomposition.dat`**: A file mapping refractive indices to specific cells [1].
*   **`dielectric_table.txt`**: A table of dielectric constants, useful for frequency-dependent materials [1].
*   **Individual `.m` files**: Specific files for different dielectric materials used in the scatterer [1].

---

## **Running Simulations**

MIDAS is designed for high-performance computing and utilizes **MPI (Message Passing Interface)** for parallel processing.

### **Execution Command**
To run the simulation, use `mpirun` to distribute the workload across multiple CPU cores:

```bash
mpirun -np <num_cores> ./midas_executable < inputs/Simulation_data.dat