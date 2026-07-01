# References

This page collects external resources and background literature for the methods used in MIDAS.

## Primary user material

- **MIDAS Overview** — see [`userguides/MIDAS_Overview_5-2026.pdf`](../userguides/MIDAS_Overview_5-2026.pdf) for the most complete description of the method, inputs, and outputs.

## Methods and algorithms

MIDAS combines several established numerical-EM techniques. The categories below indicate where to look for background; specific citations are maintained in the userguide.

- **Volume-integral-equation Method of Moments (MoM)** for dielectric scatterers discretized into cubic cells.
- **Characteristic Basis Function Method (CBFM-E)** — domain-decomposition / macro-basis reduction of the MoM system, with the reduced matrix solved via ScaLAPACK.
- **Sparsification / acceleration** used during CBF generation (sparsified near-interaction matrices solved with PARDISO).
- **Angular integration / quadrature schemes** for averaging over incident and scattering directions:
  - Gauss–Legendre quadrature (`gl`)
  - Trapezoid and Simpson rules (`tr`, `sm`)
  - Spherical *t*-designs (`sd`) — e.g. the Hardin–Sloane and the "Efficient Spherical Designs" collections.
  - Lebedev quadrature on the sphere (`lb`).

## External tools and conventions

- Shape/geometry files follow the **DDSCAT** dipole-lattice convention (`shape.dat`).
- Scattering-matrix and angle conventions follow the forward-direction-at-θ=0 convention used by DDSCAT and Mie codes.

---

*This page is a living index. Detailed bibliographic citations are kept in the MIDAS userguide; please open an issue if a specific reference is needed here.*