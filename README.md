# MIDAS (MoM Integral-equation Domain-decomposition for Arbitrarily-shaped Scatterers)

MIDAS is an implementation of a numerically efficient 3D full-wave model for electromagnetic (EM) scattering from complex-shaped scatterers.
The code is written in FORTRAN and is MPI (Message Passing Interface) parallelized. 
Several tools, written in MATLAB, are also provided to prepare and write the FORTRAN code inputs (shape, dielectric composition, simulation data files, etc.), or read and
display its outputs (scattered fields, scattering matrices and efficiency factors, etc.).

---

## Repository layout

```
midas/
+- src/                    # source code
+- examples/               # Examples of scattering calculations
+- scripts/                # Matlab and Python scripts to pre/post process data
+- README.md
```

---

## Docs
- [Documentation Home](./docs/)
- [Getting Started](./docs/getting_started.md)

## License information
----------------------

See the file ``LICENSE.txt`` for terms & conditions for usage, and a DISCLAIMER OF ALL
WARRANTIES.

DISCLAIMER: This version of MIDAS is under peer review. Please use this software with caution, ask for assistance if needed, and let us know any feedback you may have.

Copyright (c) 2026 Ines Fenni, Helene Roussel.


