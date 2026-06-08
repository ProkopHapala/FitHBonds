# HBondPlay — Fragment Interaction Energy Viewer

A WebGL-based interactive tool for exploring intermolecular interaction energies between two molecular fragments. It visualises how the total energy changes as a mobile fragment is moved around a fixed fragment, and supports polar (angle–distance) and cartesian (x–y) scanning modes.

---

## Quick Start

Open this folder via a local HTTP server (WebGL shaders must be fetched via HTTP):

```bash
cd web/HBondPlay
python3 -m http.server 8000
```

Then open `http://localhost:8000` in a modern browser with WebGL2 support.

---

## What It Shows

- **2D energy map** (top): interaction energy as a function of fragment-2 position
  - **Polar mode**: x-axis = approach angle, y-axis = distance from pivot
  - **Cartesian mode**: x-y plane scan around the pivot atom
- **1D line profile** (bottom): energy along the vertical axis through the pivot
- **Atom overlay**: fragment-1 atoms drawn on top of the 2D view with proper coordinate transforms

---

## Controls

### Fragment Selection
- **Fragment 1 (fixed)**: the stationary molecule. Default: `CH2NH`
- **Pivot atom**: the atom in fragment 1 used as the origin for scanning. Default: `N2`
- **Fragment 2 (mobile)**: the molecule that is moved around fragment 1

### Scan Modes
- **Polar**: scan angle (-180° … +180°) vs distance (0.5 … 10.0 Å)
- **Cartesian**: scan x-y plane with half-range (default 6 Å)

### Projection
Determines which 2D plane of fragment 1 is shown in the overlay:
- XY, XZ, YZ, ZY

### Energy Terms
Toggle which contributions are included:
- **Pauli**: repulsive term `Eij * exp(-2b·dr)`
- **London**: attractive term `-2·Eij * exp(-b·dr)`
- **Coulomb**: electrostatic `332.06 * Qi * Qj / r`
- **HB corr**: placeholder for hydrogen-bond correction

### Color Scale
- **Auto range**: automatically sets `vmax = |vmin|` so the scale is symmetric around zero
- **Manual range**: use the slider to fix `±vmax` (kcal/mol)

---

## Chemistry Features

### Electron Pairs (Lone Pairs)
Check **"Add e-pairs"** to automatically generate lone-pair positions for fragment 1 based on bond topology and the octet rule. Halogens get a sigma-hole lone pair opposite the bond.

### QEq Charges
Click **"Compute QEq Charges"** above either fragment table to run the Rappé–Goddard charge-equilibration algorithm. It uses atomic electron affinities and hardness parameters from `ElementTypes.dat` and updates the `Q` column in the parameter table. The Coulomb energy term then uses these computed charges.

---

## File Structure

| File | Purpose |
|------|---------|
| `index.html` | Main page layout |
| `app.js` | WebGL setup, rendering, UI event handling |
| `molecule.js` | `Molecule` class, QEq solver, electron-pair generation, XYZ parser |
| `shader.frag` | GLSL fragment shader that evaluates the energy per pixel |
| `style.css` | Dark-theme styling |
| `data/xyz/*.xyz` | Molecular geometries |

---

## Defaults on Load

- Fragment 1: `CH2NH`, pivot = `N2`
- Fragment 2: first available frag2 molecule
- Mode: **Cartesian** with half-range 6 Å
- Polar distance range: 0.5–10.0 Å
- Energy terms: Pauli + London + Coulomb enabled
- Colormap: seismic
