# PowerGrid — CLAUDE.md

## Project Overview

PowerGrid is a CPU and GPU-accelerated C++ library for iterative MRI (Magnetic Resonance Imaging) reconstruction. Developed by the MRFIL Research Group at the University of Illinois, Urbana-Champaign. It implements gridding-based NUFFT, field-corrected DFT, SENSE, pcSENSE, time-segmentation, and low-rank reconstructions.

**License:** University of Illinois/NCSA Open Source License
**Current version:** 1.1.0

---

## Directory Structure

```
PowerGrid/
├── CMakeLists.txt                   # Main build configuration
├── CLAUDE.md                        # This file
├── README.md                        # Installation and Docker instructions
├── CHANGELOG.md                     # Version history
├── Support/                         # External headers (Catch2, ArmaExtensions)
├── cmake/Modules/                   # Custom CMake find scripts
├── docker/                          # Dockerfiles (pg, pg-hpcsdk, powergrid-dev)
└── PowerGrid/                       # Main source
    ├── Core/                        # Core data structures and includes
    │   ├── PowerGrid.h              # Public aggregate header
    │   ├── PGIncludes.h             # Internal includes, precision macros
    │   ├── PGLog.hpp                # Structured logging and progress bars
    │   ├── pgComplex.hpp            # Custom complex<T> for GPU
    │   ├── pgCol.hpp                # GPU-aware column vector (OpenACC/Metal)
    │   ├── pgMat.hpp                # GPU-aware matrix (OpenACC/Metal)
    │   ├── pgSubview_Col.hpp        # Column subview
    │   ├── griddingTypes.h          # Gridding type definitions
    │   ├── AccelerateDispatch.hpp   # Apple Accelerate/vDSP dispatch for pgCol
    │   └── Tracer.hpp               # NVTX/Instruments performance tracing
    ├── Operators/                    # Encoding operators
    │   ├── Gnufft.h/.cpp            # NUFFT gridding operator
    │   ├── Gdft.h/.cpp              # Field-corrected DFT operator
    │   ├── GdftR2.h/.cpp            # DFT for radially sampled data
    │   ├── Gfft.h/.cpp              # Uniform FFT operator
    │   ├── SENSE.h/.cpp             # Multi-coil SENSE operator
    │   ├── pcSENSE.h/.cpp           # Phase-corrected SENSE
    │   └── pcSenseTimeSeg.h/.cpp    # pcSENSE + time segmentation
    ├── Penalties/                    # Regularization penalties
    │   ├── Robject.h/.cpp           # Abstract regularization base
    │   ├── QuadPenalty.h/.cpp       # L2 quadratic regularization
    │   └── TVPenalty.h/.cpp         # Total Variation regularization
    ├── Solvers/                      # Reconstruction solvers
    │   ├── solve_pwls_pcg.hpp       # PWLS conjugate gradient solver
    │   ├── solve_grad_desc.hpp      # Gradient descent solver
    │   └── reconSolve.h/.cpp        # Reconstruction helper functions
    ├── Gridding/                     # Gridding and time segmentation
    │   ├── gridding.h/.cpp          # Core KB gridding algorithm
    │   ├── griddingSupport.h/.cpp   # Gridding helpers (DCF, etc.)
    │   └── TimeSegmentation.h/.cpp  # Time-segmented field correction
    ├── FFT/                          # FFT implementations
    │   ├── fftCPU.h/.cpp            # FFTW wrapper
    │   ├── fftGPU.h/.cpp            # cuFFT wrapper
    │   ├── fftAccelerate.h/.cpp     # Apple vDSP/Accelerate FFT wrapper
    │   ├── fftshift.hpp             # FFT shift utilities
    │   ├── ftCpu.h/.cpp             # CPU DFT implementation
    │   └── ftCpuWithGrads.h/.cpp    # DFT with gradient support
    ├── IO/                           # I/O and data processing
    │   ├── processIsmrmrd.hpp       # ISMRMRD data processing
    │   ├── processNIFTI.hpp         # NIfTI data processing
    │   ├── acqTracking.h/.cpp       # Acquisition tracking
    │   └── directRecon.h/.cpp       # Direct gridding reconstruction
    ├── Metal/                        # Apple Metal GPU backend
    │   ├── MetalNufftPipeline.h/.mm # Full GPU NUFFT pipeline
    │   ├── MetalVectorOps.h/.mm     # Metal vector algebra
    │   └── *.metal                  # Metal compute shaders
    └── Tests/                       # Unit tests (Catch2)
```

---

## Key Data Structures

| Type | File | Description |
|------|------|-------------|
| `pgComplex<T>` | Core/pgComplex.hpp | Custom complex number for GPU use |
| `pgCol<T>` | Core/pgCol.hpp | GPU-aware column vector with OpenACC/Metal dual memory |
| `pgMat<T>` | Core/pgMat.hpp | GPU-aware matrix, column-major (Armadillo-compatible) |

`pgCol` and `pgMat` track host/device state via `isOnGPU` flag and use OpenACC pragmas for data movement. Default precision is `float`; define `ENABLE_DOUBLE_PRECISION` for `double`.

---

## Core Encoding Operators

All encoding operators follow a convention:
- `operator*` — forward transform (image -> k-space)
- `operator/` — adjoint transform (k-space -> image)

| Operator | Class | Description |
|----------|-------|-------------|
| NUFFT | `Gnufft<T>` | Kaiser-Bessel gridding with LUT, CPU/GPU backends |
| Field-corrected DFT | `Gdft<T>` | Phase correction via field map |
| Radial DFT | `GdftR2<T>` | Radially sampled variant of Gdft |
| Uniform FFT | `Gfft<T>` | Cartesian FFT wrapper |
| SENSE | `SENSE<T, Tobj>` | Multi-coil sensitivity encoding; Tobj = forward operator |
| pcSENSE | `pcSENSE<T>` | Shot-by-shot phase corrected SENSE |
| Time Segmentation | `TimeSegmentation<T, Tobj>` | Field map correction via temporal interpolation |

---

## Solvers

- **`solve_pwls_pcg()`** (`Solvers/solve_pwls_pcg.hpp`): Preconditioned Weighted Least Squares via conjugate gradient. Solves `min ||W^1/2(Ax-y)||^2 + R(x)`. Uses Polak-Ribiere-Polyak beta.
- **`solve_grad_desc()`** (`Solvers/solve_grad_desc.hpp`): Gradient descent alternative.

---

## Regularization (IRT Pattern)

All regularization operators (in `Penalties/`) inherit from `Robject<T>`:
- `QuadPenalty<T>`: L2 quadratic penalty
- `TVPenalty<T>`: Smooth TV penalty

---

## Build System

- **CMake >= 3.17**, C++17
- **GPU:** OpenACC with CUDA 12.6 (NVHPC/PGI compiler), or Apple Metal (clang++)
- **CPU:** GCC or Clang

### Key CMake Options

| Flag | Effect |
|------|--------|
| `-DOPENACC_GPU=ON` | Enable GPU acceleration |
| `-DOPENACC_MP=ON` | Enable OpenACC multi-threading |
| `-DMPISupport=ON` | Build MPI-distributed variants |
| `-DMETAL_COMPUTE=ON` | Enable Apple Metal GPU backend |
| `-DENABLE_DOUBLE_PRECISION=ON` | Switch to double precision |

### Main Build Targets

- **`PGCommon`** — shared library with all algorithms
- **`PowerGridIsmrmrd`**, **`PowerGridPcSense`**, etc. — reconstruction executables
- **`metal_tests`** — Metal GPU test suite (Catch2)
- **`cpu_tests`** — CPU-only test suite (Catch2)
- **`tests`** — Docker/standalone test suite (Catch2)

### Dependencies

| Library | Purpose |
|---------|---------|
| Armadillo (~9.2) | MATLAB-like linear algebra |
| ISMRMRD (v1.4) | HDF5-based MRI data format |
| FFTW3 | CPU FFT |
| CUDA 12.6 / cuFFT | GPU FFT (OpenACC path) |
| Boost (>= 1.43) | Program options, MPI bindings |
| HDF5 | Data I/O |
| SuperLU 5 | Sparse matrix decompositions |

---

## Testing

**Framework:** Catch2 (header-only, at `Support/Catch/`)

**Test files** (`PowerGrid/Tests/`):
- `pgComplexTests.cpp` — complex number arithmetic
- `pgColTests.cpp` — vector operations and memory management
- `pgColCplxTests.cpp` — complex-valued vector tests
- `pgMatTests.cpp` — matrix operations
- `operatorTests.cpp` — Gdft, Gnufft, SENSE, PCG operator tests
- `utilityTests.cpp` — bessi0, LUT, fftshift, Robject tests
- `SyntheticReconTests.cpp` — 2D synthetic reconstruction tests
- `Spiral3DReconTests.cpp` — 3D spiral reconstruction tests
- `MetalNufftTests.cpp` — Metal NUFFT pipeline tests

---

## Active Branches

- **`feature/metal-compute`** — Apple Metal GPU backend
- **`HPCSDK_25.1`** — main development branch (PR target)
- **Main branch (PRs):** `master`
