# FitHBonds — Status & Usage Guide

## Overview

FitHBonds is a C++ library with Python bindings for fitting hydrogen bond interaction models (Lennard-Jones, Morse) to quantum chemistry reference data. It was extracted from the larger FireCore repository with minimal dependencies.

---

## What Works

### Core C++ Library (`libFitREQ_lib.so`)

- **Compilation**: Builds successfully as a shared library using `g++ -std=c++20`
- **Atom type loading**: Reads `ElementTypes.dat` and `AtomTypes.dat` parameter files
- **XYZ loading**: Parses multi-sample XYZ files with energy annotations (`Etot` headers)
- **Electron pair generation**: Automatically adds epairs and sigma holes via `MMFFBuilderBase`
- **Energy evaluation**: LJ (`imodel=1`) and Morse (`imodel=2`) models
- **DOF optimization**: Gradient descent and FIRE minimization algorithms
- **OpenMP parallelization**: Multi-threaded sample evaluation

### Python Wrapper (`pyBall.FitREQ`)

- **Library loading**: Automatic `ctypes` loading from `cpp/Build/libs/Molecular/`
- **Type initialization**: `fit.loadTypes()` finds data files relative to module location
- **DOF selection**: `fit.loadDOFSelection()` loads parameter bounds
- **XYZ loading**: `fit.loadXYZ()` with optional epair insertion
- **Energy evaluation**: `fit.getEs()` returns per-sample energies
- **Optimization**: `fit.run()` with GD/MD/FIRE algorithms
- **Utilities**: `atomicUtils` for molecular manipulation

### OpenCL/GPU Path (`pyBall.OCL`)

- **FittingDriver**: Compiles templated OpenCL kernels for energy/derivative evaluation
- **Model macros**: `ENERGY_MorseQ_PAIR`, `ENERGY_LJQH2_PAIR`, `MODEL_MorseQ_PAIR`, etc.
- **Derivative tests**: Compares CPU vs GPU analytical derivatives

### Test Scripts

| Script | Status | Notes |
|--------|--------|-------|
| `opt_mini.py` | Works | CPU-only optimization; blocks on `plt.show()` |
| `check_fitREQ_ocl_cpp_derivs_.py` | Partial | CPU path works; GPU path needs matching XYZ/DOF data |
| `check_fitREQ_ocl_cpp.py` | Partial | Needs matching atom types in XYZ and DOF file |
| `check_fitREQ_derivs_ocl.py` | Partial | OCL derivative comparison |
| `plot_ocl_energy.py` | Partial | Needs appropriate XYZ/DOF files |
| `opt_2D_2.py` | Partial | Needs reference data directory |

---

## How to Use

### Quick Start — CPU Energy Evaluation

```python
import sys
sys.path.insert(0, '/path/to/FitHBonds')
from pyBall import FitREQ as fit
import numpy as np

# 1. Load atom type parameters
fit.setVerbosity(0)
fit.loadTypes()   # Finds data/ElementTypes.dat and data/AtomTypes.dat

# 2. Load DOF selection (parameter bounds)
fit.loadDOFSelection('dofSelection_H2O_Morse.dat')

# 3. Load reference geometries
nbatch = fit.loadXYZ('dimer_scan.xyz', bAddEpairs=True, bOutXYZ=False)

# 4. Read reference energies
Erefs, x0s = fit.read_xyz_data('dimer_scan.xyz')

# 5. Setup model
fit.setup(imodel=1, EvalJ=1, WriteJ=1, Regularize=1)

# 6. Set sample weights and buffers
weights = np.ones(len(Erefs)) * 0.5
fit.setWeights(weights)
fit.getBuffs()
fit.setFilter(EmodelCut=0.5, PrintOverRepulsive=-1, DiscardOverRepulsive=-1)

# 7. Evaluate model
Eerr, Es, Fs = fit.getEs(bOmp=False, bDOFtoTypes=False, bEs=True, bFs=False)
print(f"RMSE: {np.sqrt(np.mean((Es - Erefs)**2)):.4f} eV")

# 8. Optimize parameters
fit.run(iparallel=0, ialg=1, nstep=1000, Fmax=1e-4, dt=0.1, damping=0.1)

# 9. Re-evaluate
Eerr, Es_opt, _ = fit.getEs()
print(f"Optimized RMSE: {np.sqrt(np.mean((Es_opt - Erefs)**2)):.4f} eV")
```

### Building the Shared Library

```bash
cd /path/to/FitHBonds/cpp

# Debug build (default)
cmake -B build -S .
cmake --build build -j4

# Release build
 cmake -B build -S . -DRELEASE=ON
 cmake --build build -j4

# With AddressSanitizer
 cmake -B build -S . -DWITH_ASAN=ON
 cmake --build build -j4

# With OpenMP
 cmake -B build -S . -DWITH_OMP=ON
 cmake --build build -j4
```

The shared library is written to `Build/libs/Molecular/libFitREQ_lib.so`.

### Using the OpenCL Path

```python
from pyBall.OCL.NonBondFitting import setup_driver
from pyBall.OCL.FittingDriver import FittingDriver

# Create driver and compile kernel
drv = setup_driver(model_name='ENERGY_MorseQ_PAIR')

# Load data
drv.load_data('dimer_scan.xyz')
drv.load_dofs('dofSelection.dat')
drv.init_and_upload()

# Evaluate energy
Emols = drv.evaluate_objective(initial_dofs)
```

---

## File Structure

```
FitHBonds/
├── cl/                          # OpenCL kernel files
│   ├── FitREQ.cl
│   └── Forces.cl
├── cpp/
│   ├── Build/libs/Molecular/    # Compiled shared library
│   ├── common_resources/          # Symlinks to data/ and cl/
│   ├── libs/FitREQ_lib.cpp      # C++ library entry point
│   ├── math/                    # Math headers (Vec3, Mat3, etc.)
│   └── molecular/               # Molecular headers (FitREQ, MMFF, Atoms)
├── data/
│   ├── AtomTypes.dat            # MMFF atom type parameters
│   └── ElementTypes.dat         # Element parameters
├── doc/                         # This documentation
├── pyBall/                      # Python package
│   ├── FitREQ.py                # Main Python wrapper
│   ├── atomicUtils.py           # Molecular utilities
│   ├── cpp_utils_.py            # C++ compilation/loading utilities
│   ├── elements.py              # Periodic table data
│   └── OCL/                     # OpenCL Python wrappers
│       ├── FittingDriver.py
│       ├── NonBondFitting.py
│       ├── OpenCLBase.py
│       └── clUtils.py
└── tests/tFitREQ/               # Test scripts and data
    ├── wb97m_input/             # Reference XYZ scans (203 files)
    ├── dofSelection_*.dat       # DOF parameter bounds
    └── *.py                     # Test scripts
```

---

## What Still Needs to Be Done

### 1. Data File Compatibility

**Problem**: The DOF selection files (`dofSelection_H2O_Morse.dat`, etc.) specify atom types for a specific system (e.g., H₂O), but the available XYZ files are mixed complexes (`C4H3NO2-D1_H2O-A1.xyz`, etc.) with additional atom types.

**Impact**: `fit.setup()` works, but OpenCL `FittingDriver` fails with:
```
prepare_host_data() ERROR: len(dof_to_atom_list)==0 => Exiting.
```

**Solution needed**:
- Create DOF selection files matching the mixed complexes, **or**
- Filter XYZ files to only include atom types present in the DOF file, **or**
- Use pure H₂O dimer XYZ files with the H₂O DOF files

### 2. Missing Reference Data

Some test scripts reference data not present in this repo:
- `opt_2D_2.py`: `ref_path = "/home/prokop/Desktop/CARBSIS/PEOPLE/Paolo/FitREQ/DFT_2D/"`
- `opt_mini.py`: `input_2CH2NH.xyz`, `dofSelection_CH2NH_LJ.dat`
- `check_fitREQ_energy_ocl.py`: `HHalogens/porcessed/HBr-D1_HBr-A1.xyz`

**Solution**: Copy from FireCore or generate matching test data.

### 3. Hardcoded GUI Dependencies

Scripts like `opt_mini.py` call `plt.show()` which blocks execution in headless environments.

**Workaround**: Set `fit.plt = None` or use `MPLBACKEND=Agg`.

### 4. Format Warnings in C++

`sscanf`/`printf` format specifier warnings (`%i` for `uint8_t*`, `%p` for `std::vector`). These are cosmetic but should be cleaned up for C++20 compliance.

### 5. CMake Build System

The original FireCore uses `make`; a dedicated `CMakeLists.txt` for FitHBonds alone would be cleaner.

### 6. Additional Test Scripts Needed

| Script | Source | Status |
|--------|--------|--------|
| `opt_2D.py` | FireCore | Not copied |
| `check_fitREQ_derivs.py` | FireCore | Not copied |
| `check_fitREQ_ocl_cpp copy.py` | FireCore | Not copied |
| `test_export.py` | FireCore | Not copied |
| `fit_manual*.py` series | FireCore | Not copied |

---

## Known Limitations

- **C++20 required**: `std::unordered_set::contains()` is used in `MMFFBuilderBase.h`
- **OpenCL optional**: CPU path works without pyopencl; GPU path requires it
- **Single-threaded library loading**: `ctypes.CDLL` with `RTLD_LOCAL` prevents multiple loads
- **No automatic recompilation**: Set `CPP_RECOMPILE=1` or delete `.so` to force rebuild

---

## Quick Verification

```bash
cd /path/to/FitHBonds

# Test C++ library compilation
python3 -c "
import sys; sys.path.insert(0, '.')
from pyBall import FitREQ as fit
fit.loadTypes()
print('Library OK')
"

# Test with actual data
python3 -c "
import sys; sys.path.insert(0, '.')
from pyBall import FitREQ as fit
fit.loadTypes()
fit.loadDOFSelection('tests/tFitREQ/dofSelection_H2O_Morse.dat')
fit.loadXYZ('tests/tFitREQ/wb97m_input/C4H3NO2-D1_H2O-A1.xyz', bAddEpairs=True)
Erefs, x0s = fit.read_xyz_data('tests/tFitREQ/wb97m_input/C4H3NO2-D1_H2O-A1.xyz')
fit.setup(imodel=1, EvalJ=1, WriteJ=1, Regularize=1)
Eerr, Es, Fs = fit.getEs()
print(f'Loaded {len(Erefs)} samples, RMSE = {((Es - Erefs)**2).mean()**0.5:.3f} eV')
"
```
