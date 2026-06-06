# FitHBonds

Hydrogen bond interaction model fitting — extracted from the [FireCore](https://github.com/ProkopHapala/FireCore) molecular simulation framework with minimal dependencies.

This repository provides tools for fitting classical force-field parameters (Lennard-Jones and Morse potentials) to quantum-chemistry reference data of hydrogen-bonded dimers. It supports both CPU evaluation via a C++ shared library and GPU evaluation via OpenCL kernels.

## What this repo does

Given a set of dimer interaction energy curves (e.g., from DFT / wB97M-V), FitHBonds:
1. **Builds molecular topology** — adds bonds, electron pairs, and sigma holes automatically
2. **Evaluates interaction energies** using LJ or Morse models with hydrogen-bond corrections
3. **Optimizes force-field parameters** (REQH: Radius, Energy, charge, H-bond strength) to minimize the RMS error vs reference
4. **Compares CPU and GPU results** — analytical derivatives checked against numerical finite differences

The code is header-only C++ (except the single `FitREQ_lib.cpp` wrapper) with Python bindings via `ctypes`.

## Two executable entry points

### 1. C++ shared library (`libFitREQ_lib.so`) — CPU path

This is the primary entry point. Python loads the compiled shared library and calls C++ functions directly.

```python
from pyBall import FitREQ as fit

fit.loadTypes()                                    # load ElementTypes.dat / AtomTypes.dat
fit.loadDOFSelection('dofSelection_H2O_Morse.dat')  # parameter bounds
fit.loadXYZ('dimer_scan.xyz', bAddEpairs=True)     # load reference geometries
fit.setup(imodel=1, EvalJ=1, WriteJ=1)             # choose LJ (1) or Morse (2)
Eerr, Es, Fs = fit.getEs()                          # evaluate all samples
fit.run(nstep=1000, Fmax=1e-4)                     # optimize DOFs
```

Key files:
- `cpp/libs/FitREQ_lib.cpp` — C++ library source
- `pyBall/FitREQ.py` — Python `ctypes` wrapper
- `cpp/Build/libs/Molecular/libFitREQ_lib.so` — compiled shared library

### 2. OpenCL GPU kernels (`pyBall.OCL`) — GPU path

For much faster energy/derivative evaluation on GPU. The OpenCL driver compiles templated kernels from `cl/Forces.cl` and `cl/FitREQ.cl`.

```python
from pyBall.OCL.NonBondFitting import setup_driver
from pyBall.OCL.FittingDriver import FittingDriver

drv = setup_driver(model_name='ENERGY_MorseQ_PAIR')
drv.load_data('dimer_scan.xyz')
drv.load_dofs('dofSelection.dat')
drv.init_and_upload()

# Evaluate energy for given DOF values
Emols = drv.evaluate_objective(initial_dofs)
```

Key files:
- `cl/FitREQ.cl` / `cl/Forces.cl` — OpenCL kernel templates
- `pyBall/OCL/FittingDriver.py` — OpenCL driver
- `pyBall/OCL/NonBondFitting.py` — high-level fitting wrappers

## Build instructions

### Prerequisites
- `g++` with C++20 support (tested on GCC 11+)
- `cmake` >= 3.10
- Python 3 with `numpy` and `matplotlib`
- Optional: `pyopencl` for GPU acceleration

### Build the C++ library

```bash
cd cpp

# Debug build (default)
cmake -B build -S .
cmake --build build -j4

# Release build (optimized)
cmake -B build -S . -DRELEASE=ON
cmake --build build -j4

# With AddressSanitizer (debug memory issues)
cmake -B build -S . -DWITH_ASAN=ON
cmake --build build -j4

# With OpenMP (multi-threaded evaluation)
cmake -B build -S . -DWITH_OMP=ON
cmake --build build -j4
```

The shared library is written to `cpp/Build/libs/Molecular/libFitREQ_lib.so`.

## Running tests

### Quick verification (CPU only)

```bash
cd tests/tFitREQ
python3 -c "
import sys; sys.path.insert(0, '../..')
from pyBall import FitREQ as fit
fit.loadTypes()
fit.loadDOFSelection('dofSelection_H2O_Morse.dat')
fit.loadXYZ('wb97m_input/C4H3NO2-D1_H2O-A1.xyz', bAddEpairs=True)
Erefs, x0s = fit.read_xyz_data('wb97m_input/C4H3NO2-D1_H2O-A1.xyz')
fit.setup(imodel=1, EvalJ=1, WriteJ=1, Regularize=1)
Eerr, Es, Fs = fit.getEs()
import numpy as np
print(f'Loaded {len(Erefs)} samples, RMSE = {np.sqrt(np.mean((Es - Erefs)**2)):.3f} eV')
"
```

### Optimization test (`opt_mini.py`)

This is the simplest end-to-end test. It loads a dimer scan, evaluates the model, and optionally optimizes parameters.

```bash
cd tests/tFitREQ
# The script uses matplotlib; for headless runs set:
#   MPLBACKEND=Agg  or  fit.plt = None
python3 opt_mini.py
```

### CPU vs GPU derivative consistency (`check_fitREQ_ocl_cpp_derivs_.py`)

Compares CPU and OpenCL analytical derivatives. Needs matching DOF/XYZ files (same atom types in both).

```bash
cd tests/tFitREQ
python3 check_fitREQ_ocl_cpp_derivs_.py \
    --xyz wb97m_input/C4H3NO2-D1_H2O-A1.xyz \
    --dof dofSelection_H2O_Morse.dat \
    --npts 10 --morse 1 --verbose 0
```

> **Note**: This script's GPU path currently fails if the XYZ contains atom types not present in the DOF file. Use a matching pure-dimer XYZ (e.g. H₂O–H₂O) for the full CPU+GPU comparison.

## File layout

```
FitHBonds/
├── cl/                          # OpenCL kernel files
│   ├── FitREQ.cl
│   └── Forces.cl
├── cpp/
│   ├── libs/FitREQ_lib.cpp      # C++ library entry point
│   ├── math/                    # Vec3, Mat3, fastmath, etc.
│   ├── molecular/               # Atoms, MMFF, Builder, FitREQ headers
│   ├── CMakeLists.txt           # minimal CMake (RELEASE / ASAN / OpenMP)
│   └── Build/libs/Molecular/    # compiled libFitREQ_lib.so
├── data/
│   ├── AtomTypes.dat            # MMFF atom parameters
│   └── ElementTypes.dat         # Element parameters
├── pyBall/                      # Python package
│   ├── FitREQ.py                # ctypes wrapper for C++ lib
│   ├── atomicUtils.py           # molecular utilities
│   ├── cpp_utils_.py            # C++ compilation helpers
│   └── OCL/                     # OpenCL Python wrappers
├── tests/tFitREQ/               # test scripts + reference XYZ scans
│   ├── wb97m_input/             # 203 dimer scan files
│   ├── dofSelection_*.dat       # DOF bounds
│   └── *.py                     # test scripts
├── tests/RigidAtomHfunc/        # rigid atom angular function tests
│   ├── plot_H2O_sp3.py          # visualize sp3 hybrid orbitals of water
│   └── h2o_sp3_viewer.html      # WebGL interactive viewer for water sp3 orbitals
└── doc/status.md                # detailed status & known issues
```

## Known limitations

- **C++20 required**: `std::unordered_set::contains()` is used in `MMFFBuilderBase.h`
- **OpenCL optional**: CPU path works without `pyopencl`; GPU path requires it
- **Data compatibility**: DOF files must match atom types in XYZ files for GPU tests
- **Library reloading**: `ctypes.CDLL` with `RTLD_LOCAL` prevents multiple loads in one Python session

See [doc/status.md](doc/status.md) for the full list of what works, what doesn't, and future todos.
