# References

This page collects external resources and background literature for the methods used in MIDAS.

## Primary user material

- **MIDAS Overview** — see [`userguides/MIDAS_Overview_5-2026.pdf`](../userguides/MIDAS_Overview_5-2026.pdf) for the most complete description of the method, inputs, and outputs.

## Methods and algorithms

MIDAS combines several established numerical-EM techniques. The categories below indicate where to look for background; specific citations are maintained in the userguide.

- **Volume-integral-equation Method of Moments (MoM)** for dielectric scatterers discretized into cubic cells.
- **Characteristic Basis Function Method (CBFM-E)** — domain-decomposition / macro-basis reduction of the MoM system, with the reduced matrix solved via ScaLAPACK.
- **Angular integration / quadrature schemes** for averaging over incident and scattering directions:
  - Spherical *t*-designs (`sd`) — e.g. the Hardin–Sloane and the "Efficient Spherical Designs" collections.
  - Lebedev quadrature on the sphere (`lb`).

## Foundational method papers

Papers that describe the original numerical-method concepts underlying MIDAS.

- Nguyen H, Roussel H, Tabbara W. A coherent model of forest scattering and SAR imaging in the VHF and UHF-band. IEEE transactions on geoscience and remote sensing. 2006 Mar 27;44(4):838-48.

- Fenni, I., Roussel, H., Darces, M. and Mittra, R., 2014. Fast analysis of large 3-D dielectric scattering problems arising in remote sensing of forest areas using the CBFM. IEEE Transactions on Antennas and Propagation, 62(8), pp.4282-4291.

## Applications of MIDAS

Papers in which MIDAS was used to compute the scattering properties of arbitrarily shaped objects.

- Fenni, I., Haddad, Z.S., Roussel, H., Kuo, K.S. and Mittra, R., 2017. A computationally efficient 3-D full-wave model for coherent EM scattering from complex-geometry hydrometeors based on MoM/CBFM-enhanced algorithm. IEEE Transactions on Geoscience and Remote Sensing, 56(5), pp.2674-2688.

- Fenni, I., Kuo, K.S., Haynes, M.S., Haddad, Z.S. and Roussel, H., 2021. Evaluation of higher‐order quadrature schemes in improving computational efficiency for orientation‐averaged single‐scattering properties of nonspherical ice particles. *Journal of Geophysical Research: Atmospheres*, 126(11), p.e2020JD034172.

- Haynes, M.S., Fenni, I. and Davidsson, B.J., 2023. Inverse scattering under the born approximation using an object T-Matrix and full bistatic spherical sampling geometry. IEEE Transactions on Antennas and Propagation, 72(1), pp.862-876.

## External tools and conventions

- Shape/geometry files follow the **DDSCAT** dipole-lattice convention (`shape.dat`).
- Scattering-matrix and angle conventions follow the forward-direction-at-θ=0 convention used by DDSCAT and Mie codes.

---

*This page is a living index. Please open an issue if a specific reference is needed here.*