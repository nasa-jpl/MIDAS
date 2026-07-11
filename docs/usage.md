# Using MIDAS

This document describes how to configure and run MIDAS to compute electromagnetic scattering by complex-shaped particles.

---

## Directory Structure

MIDAS expects all input files to reside in an `inputs/` subdirectory relative to the working directory from which the executable is launched:

```
run_directory/
├── inputs/
│   ├── Simulation_data.dat       ← required: main configuration
│   ├── shape.dat  (or Cells.dat) ← required: particle geometry
│   ├── diel/                     ← refractive index files (.m files, table, or composition file)
│   └── SimShapes.dat             ← required only for batch shape runs
```

Output is written to a subfolder whose name is derived from the scatterer identity and size (e.g., `Results/p0004-ap=0.4mm`).

---

## Input Files

### 1. `Simulation_data.dat` — Main Configuration

This is the primary input file read by `Get_InputData`. The sections below follow the exact reading order of the code.

#### EM Wave Parameters

Specify either **frequency** or **wavelength** by starting this section with `NFreq` or `NWave`.

**Frequency mode (`NFreq`):**
```
NFreq
<number of frequencies>
minmax  (or: list)
frequency unit = GHz   (options: MHz | GHz | THz)
<freq_min>,<freq_max>  (if minmax) — or one value per line (if list)
```

**Wavelength mode (`NWave`):**
```
NWave
<number of wavelengths>
minmax  (or: list)
wavelength unit = mm   (options: m | mm | um)
<wave_min>
<wave_max>   (if minmax) — or one value per line (if list)
```

Unit mapping: `MHz ↔ m`, `GHz ↔ mm`, `THz ↔ um`.

---

#### Scatterer Geometry

```
<shape_list>      ! 0 = single shape; 1 = batch from SimShapes.dat; 2 = batch with suffix
<Outfld_name>     ! name of the root output folder

<type_s> <info_s> [cells]
! type_s values:
!   1 = sphere 
!   2 = complex shape read from shape.dat
!   3 = cylinder 

<ap>              ! effective radius (in mm if GHz, um if THz, m if MHz)
                  ! for cylinders: <ac> <lc>  (radius and length)
```

---

#### Dielectric Properties

```
<dielcomp_option>  <Ndiel>
```

| `dielcomp_option` | Description | Additional lines |
|---|---|---|
| `fromonlymfile` | Single homogeneous material; one `.m` file | 1 filename |
| `fromshapefile` | Each cell's material index is read from shape.dat | `Ndiel` filenames |
| `fromdielcompositionfile` | Per-cell refractive index from `inputs/dielcomposition.dat` | none |
| `fromdieltable` | Dielectric table from `inputs/dielectric_table.txt` | none |

> **Note:** `fromdielcompositionfile`, `fromshapefile`, and `fromdieltable` require `type_s = 2`.

**`.m` file format** (3 header lines, then data):
```
<header line 1>
<header line 2>
<header line 3>
<lambda>  <m_real>  <m_imag>
...
```
One row per simulation frequency, or a single row applied to all frequencies.

**`dielectric_table.txt` format:**
- First line: header containing `eps` (permittivity) or `m` (refractive index) label.
- Subsequent lines: one row per material, columns are `index  mr1 mi1  mr2 mi2 ...` for each frequency.
- If the header contains `eps`, the code automatically converts to refractive index.

**`dielcomposition.dat` format:**
```
<header>
Ndiel = <Ndiel> ; Nfreq = <Nfreq_dielfile>
<header>
<cell_id>  mr1 mi1  mr2 mi2  ...   (one row per cell)
```

---

#### Discretization

```
S
<Sc>    ! cell size (genrellay noted d) in um (THz/GHz) or m (MHz)
```

This cell size is used only when `type_s` ≠ 2. When `type_s = 2`, the cell size `d` (`Sc`) is instead derived from the effective radius `ap` by matching the total volume of the `N` discretization cells to that of the volume-equivalent sphere:

$$N\,d^3 = \frac{4}{3}\pi\,a_p^3 \quad\Longrightarrow\quad d = a_p\left(\frac{4\pi}{3N}\right)^{1/3}$$

The geometry is discretized once at the highest requested frequency.

---

#### Numerical Methods

```
<Nber_methods>
<method_1>
<method_2>
...
```

Available method names:

| Name | Description |
|---|---|
| `CBFM-E` | Characteristic Basis Function Method|
| `MoM` | Method of Moments (full matrix) |

---

#### Transmitters and Receivers

```
<NumIntType_t> / <NumIntType_r>   ! incident / scattering direction distribution
<Ninc_sugg>                       ! suggested number of incident directions
<Nscat_sugg>                      ! suggested number of scattering directions

! Transmitters (incident directions):
<theta_min>  <theta_max>  <NTrTheta>
<phi_min>    <phi_max>    <NTrPhi>

! Receivers (scattering directions):
<theta_min>  <theta_max>  <NRxTheta>
<phi_min>    <phi_max>    <NRxPhi>
<NPolBeta>

<wr_Sij>    ! 1 = write scattering matrices to S_files/
<wr_Qij>    ! 1 = write Q (efficiency) files to Q_files/
```

Direction distribution types (`NumIntType_t` for transmitters, `NumIntType_r` for receivers — both use the same codes):

| Code | Type | Angular averaging used |
|---|---|---|
| `un` | Uniform step in theta and phi | Adaptive quadrature over the uniform grid |
| `sd` | Spherical T-design | Equal-weight sum over design points |
| `lb` | Lebedev quadrature | Lebedev weights |
| `rf` | Custom directions read from `inputs/IncScattDirs.dat` | — |

> If `NTr = 1`, both `wr_Sij` and `wr_Qij` are automatically set to 1.

---

#### CBFM Parameters

```
<div_type>       ! block division strategy: SPH | CSH
<Navg_cells>     ! target average number of cells per block

<NipwsType>      ! incident PW distribution for CBF generation: un | sd | lb
<Nipws>          ! number of incident plane waves for CBF generation
```

Block division strategies (`div_type`):

| Code | Strategy |
|---|---|
| `SPH` | Hierarchical octree subdivision — intended for sphere/Chebyshev-like compact shapes |
| `CSH` | Geometry-driven adaptive division along the principal axes — intended for complex/arbitrary shapes |

Incident plane-wave distribution for CBF generation (`NipwsType`) maps to the same families as the receiver/transmitter codes: `un` (uniform in cos θ and φ), `sd` (spherical T-design), `lb` (Lebedev). `Nipws` is a *suggested* count; MIDAS selects the closest available design/quadrature order, and (when `set_Nipws ≠ 0`) may auto-tune it from the block size and dielectric contrast.

---

#### Far-Field / Near-Field Mode

```
<CextIntFields>  ! 0 or 1 — store extinction cross section from internal field
<FFA>            ! 0 = near-field (writes Ei/Es); 1 = far-field (writes S and Q)
<Rso>            ! observation distance
```

---

#### Save Options

```
<verbose>          ! verbosity level
<save_Zc>          ! 1 = save impedance matrix to Sol_files/
<save_Eint>        ! 1 = save internal electric field to Sol_files/
<save_Eint_Nmax>   ! maximum Nbc for which Eint is saved (memory guard)
<save_Einc>        ! 1 = save incident field (useful for inversion algorithms)
```

---

### 2. Geometry Files

#### `shape.dat` (for `type_s = 2`)

`shape.dat` follows the **DDSCAT** dipole-lattice convention and is read by `Read_ShapeFile`. Both the number of header lines and the number of data columns are auto-detected by `Inspect_ShapeFile`, so no format flag is needed:

- Any line that contains a letter, or is blank, is treated as a **header** line.
- The **column count** is inferred from the first purely numeric line.

**Reading order:**

1. **Line 1** — free-text comment (ignored).
2. **Line 2** — the total number of cells; a trailing label after the count (e.g. `= NAT`) is allowed and ignored.
3. **Any remaining header lines** are skipped.
4. **One row per cell** follows, parsed according to the detected column count.

**Accepted column formats (auto-detected):**

| Columns | Format | Notes |
|---|---|---|
| 7 | `num  ix  iy  iz  mx  my  mz` | cell id, integer lattice indices, per-axis material indices (anisotropic) |
| 5 | `num  ix  iy  iz  m_index` | cell id, integer lattice indices, single material index |
| 4 | `ix  iy  iz  m_index` | integer lattice indices, single material index (no cell id) |

The `ix iy iz` are **integer lattice coordinates**; the physical cell size `d` is applied separately (see [Discretization](#discretization)). The material index (`m_index`, or `mx my mz`) selects the corresponding `.m` file when `dielcomp_option = fromshapefile`.

#### `Cells.dat` (alternative geometry, activated with `cells` keyword)

One row per cell:
```
x(m)  y(m)  z(m)  Sc(m)  m_index  block_num
```

#### `SimShapes.dat` (for `shape_list ≥ 1`)

One entry per line:
```
path/to/shape.dat : info_label  effective_radius
```
Example:
```
inputs/shapes/p0004/shape.dat : p0004  0.400000000
```

---

## Running MIDAS

MIDAS uses MPI for parallelism. Launch with:

```bash
mpirun -np <Nprocs> ./midas
```

All processes read the same `inputs/Simulation_data.dat` automatically.

**Important constraint for CBFM/MLCBFM:**
The number of MPI processes must not exceed the number of blocks (`Nblocks`). If `Nprocs > Nblocks`, the code prints a performance error and exits:
```
Performance Error : Nprocs = X > Nblocks = Y ! Please restart with fewer processors !
```

The geometry is discretized once (at the highest frequency) and then the solver loops over all requested frequencies.

---

## Output Files

All output is written to a subfolder of `<Outfld_name>`, named automatically from the scatterer type and size:

| Folder / File | Contents | Condition |
|---|---|---|
| `Simulation_data.dat` | Copy of the input configuration | always |
| `IncScattDirs.dat` | List of all incident and scattering directions | always |
| `Cells.dat` | Voxel positions, sizes, material indices, block indices | always |
| `Blocks.dat` | Block composition and extensions | CBFM/MLCBFM only |
| `S_files/` | Scattering matrix files per incident direction and frequency | `FFA=1` and `wr_Sij=1` |
| `Q_files/` | Scattering efficiency files | `FFA=1` and `wr_Qij=1` |
| `Es_files/` | Scattered electric field | `FFA=0` |
| `Ei_files/` | Incident electric field | `FFA=0` |
| `Sol_files/` | Internal fields and/or impedance matrices | `save_Eint=1` or `save_Zc=1` |
| `Analysis/` | Debug logs and memory tracking (`TrackAllocate_j<rank>.dat`) | `debug_mode=1` |

### Scattering Matrix File Naming

Text files written to `S_files/`:
```
Smtable_<freq><unit>_kt<kkt>_<method>.dat
```
Example: `Smtable_9.40GHz_kt0001_CBFM-E.dat`

Each file contains one header line with the incident direction (theta, phi), then one row per scattering direction:
```
theta  phi  Re(Svv)  Im(Svv)  Re(Svh)  Im(Svh)  Re(Shv)  Im(Shv)  Re(Shh)  Im(Shh)
```

For multi-frequency runs, files are prefixed with `Sim<n>_` (e.g., `Sim2_Smtable_...`).

### Scattered Field File Naming (`FFA = 0`)

When the scattered field path is selected, the scattered and incident fields at the receivers are written per incident direction:

```
Es_files/Esca_<freq><unit>_kt<kkt>_<method>.dat     ← scattered field
Ei_files/Einc_<freq><unit>_kt<kkt>_<method>.dat     ← incident field at receivers
```

Each file has a header row, then one row per scattering direction with the V/H polarization components:
```
theta  phi  Re(Evv)  Im(Evv)  Re(Evh)  Im(Evh)  Re(Ehv)  Im(Ehv)  Re(Ehh)  Im(Ehh)
```
As with the S-matrix files, multi-frequency runs prefix the name with `Sim<n>_`.

---

## Batch Runs Over Multiple Shapes

Set `shape_list = 1` (or `2`) in `Simulation_data.dat` and provide `inputs/SimShapes.dat`.

- `shape_list = 1`: output folder name uses `info_label` and effective radius.
- `shape_list = 2`: output folder name also appends any suffix found in the shape file path.

MIDAS will loop over all entries in `SimShapes.dat` and produce a separate output subfolder for each shape.

---

## Dielectric Composition Summary

The `dielcomp_option` controls how a refractive index is assigned to each cell at each frequency:

```
DielComposition(m_lambdas, Cells)
```

- `fromonlymfile`: all cells share the same `m` value at each frequency.
- `fromshapefile` / `fromdielcompositionfile`: each cell has an index `n_diel` pointing to a row in `m_lambdas`.
- `fromdieltable`: `m_lambdas` is populated from a structured table file; supports both `m` and `eps` input formats.

The subroutine is called inside the frequency loop, so frequency-dependent dielectric values are correctly applied at each frequency step.
