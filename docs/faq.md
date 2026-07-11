
# Frequently Asked Questions (FAQ)

**Q1: Can I run MIDAS without MATLAB?**  
A: Yes. The FORTRAN solver reads plain text input files and writes plain text output files independently of MATLAB. The MATLAB tools are optional and only needed for input preparation and output visualization. They can be replaced with equivalent Python scripts. 

**Q2: How many MPI processes should I use?**  
A: When using `CBFM-E`, the number of MPI processes must not exceed the number of CBFM blocks (`Nblocks`). If it does, MIDAS will print an error and exit. Start with a number larger to `Nblocks`. For `MoM`, there is no such constraint.

**Q3: Is the code compatible with Windows?**  
A: Yes. MIDAS detects the operating system at runtime and uses the correct path separator (`\` on Windows, `/` on Linux). However, the code is primarily developed and tested on Linux HPC systems. Windows use is supported but less tested.

**Q4: What is the difference between `FFA=0` and `FFA=1`?**  
A: `FFA=1` enables the far-field approximation and computes amplitude scattering matrices (S matrices) and efficiency factors (Q). `FFA=0` computes and stores the full scattered (`Es_files/`) and incident (`Ei_files/`) electric fields at a distance R. Use `FFA=1` for precipitation particles and bulk optical properties; use `FFA=0` when scattered fields are needed (e.g., asteroid simulations).

**Q5: What `dielcomp_option` should I use for a homogeneous particle?**  
A: Use `fromonlymfile` with a single `.m` file. This assigns the same refractive index to all cells at each frequency.

*(This section will be expanded as more questions arise.)*
