# MIDAS

MIDAS (for **M**oM **I**ntegral-equation **D**omain-decomposition for **A**rbitrarily-shaped **S**catterers) is an implementation of a numerically efficient 3D full-wave model for electromagnetic (EM) scattering from complex-shaped scatterers.
The code is written in FORTRAN and is MPI (Message Passing Interface) parallelized. 
Several tools, written in MATLAB, are also provided to prepare and write the FORTRAN code inputs (shape, dielectric composition, simulation data files, etc.), or read and
display its outputs (scattered fields, scattering matrices and efficiency factors, etc.).

## Repository layout

```
midas/
+ src/                  # source code
+ docs/                 # Markdown-based project documentation
+ userguides/           # End-user materials and presentations
+ examples/             # Examples of scattering calculations
+ tools/                # includes various user tools (make & job files, analysis scripts)
+ README.md
+ license
```

## Docs
Use the links below to navigate the documentation:

- [Getting Started](./docs/getting_started.md) � Quickstart instructions and setup
- [Architecture](./docs/architecture.md) � Overview of code structure and workflow
- [Usage](./docs/usage.md) � How to run simulations and process outputs
- [FAQ](./docs/faq.md) � Frequently asked questions
- [References](./docs/references.md) � External resources and related literature



## Contact
Ines Fenni,
Jet Propulsion Laboratory, California Institute of Technology,
ines.fenni@jpl.nasa.gov

H�l�ne Roussel,
Sorbonne Universit�,
helene.roussel@sorbonne-universite.fr

## License information

See the file ``LICENSE.txt`` for terms & conditions for usage, and a DISCLAIMER OF ALL WARRANTIES.

DISCLAIMER: MIDAS is under active development. Features may change, and unexpected issues may occur. Please use with caution, and do not hesitate to reach out for assistance or to share feedback. Users are welcome to open issues to report problems they encounter. We will consider them carefully as time allows.



