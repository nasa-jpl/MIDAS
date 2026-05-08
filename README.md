# MIDAS

**MIDAS** (**M**oM **I**ntegral-equation **D**omain-decomposition for **A**rbitrarily-shaped **S**catterers) is a numerically efficient 3D full-wave electromagnetic (EM) scattering solver for complex-shaped scatterers.

The codebase is written in **FORTRAN** and parallelized using **MPI** (Message Passing Interface).

Several companion tools written in **MATLAB** are also provided to:

- Prepare simulation inputs  
  (geometry, dielectric composition, simulation data files, etc.)
- Read and visualize outputs  
  (scattered fields, scattering matrices, efficiency factors, etc.)

---

## Repository Layout

```text
midas/
├── src/          # Source code
├── docs/         # Markdown-based project documentation
├── userguides/   # End-user materials and presentations
├── examples/     # Example scattering calculations
├── tools/        # Utility scripts, make/job files, analysis tools
├── README.md
└── LICENSE.txt
```

---

## Documentation

Use the links below to navigate the documentation:

- [Getting Started](./docs/getting_started.md) — Quickstart instructions and setup
- [Architecture](./docs/architecture.md) — Overview of code structure and workflow
- [Usage](./docs/usage.md) — Running simulations and processing outputs
- [FAQ](./docs/faq.md) — Frequently asked questions
- [References](./docs/references.md) — External resources and related literature

---

## Contact
**Ines Fenni**
Jet Propulsion Laboratory, California Institute of Technology
ines.fenni@jpl.nasa.gov

**Hélène Roussel**  
Sorbonne Université  
helene.roussel@sorbonne-universite.fr

---

## License Information

See the `LICENSE.txt` file for terms and conditions of use, including the disclaimer of warranties.

### Disclaimer

MIDAS is under active development. Features and interfaces may evolve over time, and unexpected issues may occur.

Please use the software with care and feel free to reach out with questions or feedback. Users are encouraged to open issues to report bugs or problems encountered during use. Contributions and feedback are appreciated and will be reviewed as time permits.