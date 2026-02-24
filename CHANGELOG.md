# Change Log
All notable changes to this project will be documented in this file.

## [Unreleased]

## [v1.9.0] 2025-02-24
### Major changes toward v2.0
#### GPU Backends
* Apple Metal GPU compute backend for macOS/Apple Silicon (NUFFT pipeline, DFT, vector algebra)
* Metal vector dispatch integrated into pgCol operators with automatic CPU fallback
* Accelerate/vDSP FFT backend for macOS

#### Core Data Structures
* GPU-aware `pgCol<T>` and `pgMat<T>` column vector and matrix types with OpenACC/Metal dual memory
* Custom `pgComplex<T>` for GPU-compatible complex arithmetic
* Migrated SENSE, pcSENSE, TimeSegmentation, and PCG solver hot paths from Armadillo to pgCol/pgMat

#### TUI Monitoring (pgview)
* `pgview` terminal UI viewer for real-time reconstruction monitoring
* Structured JSONL output protocol for log, progress, metrics, and image preview messages
* Auto-spawn pgview from reconstruction executables with classic spdlog+indicators fallback
* Live image preview with 3-plane MPR (Axial/Coronal/Sagittal) for 3D volumes
* Quad-panel layout (2x2 grid) for 3D, side-by-side for 2D reconstructions
* Sparkline convergence plots (error norm + roughness penalty)
* Progress bars with ETA and wall clock completion time
* iTerm2 and Kitty inline terminal image protocol support

#### Build System
* Restructured source tree into logical subdirectories (Core, Operators, Penalties, Solvers, Gridding, FFT, IO, Metal)
* CMake semantic versioning with auto-generated Version.h
* SO versioning for libPGCommon shared library
* NVIDIA HPC SDK 25.1 / CUDA 12.6 support

#### Testing
* Catch2 unit test suite: pgComplex, pgCol, pgMat, operators, utilities
* Synthetic 2D/3D reconstruction integration tests
* Metal GPU pipeline tests and benchmarks

#### Bug Fixes
* Fixed `-ffast-math` breaking `arma::min()`/`arma::max()` on Apple clang
* Fixed pgCol copy assignment not resizing when sizes differ
* Fixed pgMat column-major indexing in `at(r,c)` and copy constructor
* Fixed Robject adjoint (Ctd) implementation
* Fixed pgCol SFINAE template parameter shadowing

## [v1.1.0] 2020-03-18
## Updating with direct gridding reconstruction and low rank reconstructions
* Implements calculation of density compensation functions in 2D/3D using method of Pipe et al.
* Implements direct gridding reconstruction of coil images and sum-of-squares (SoS) reconstruction for fully sampled data
* Implements phase corrected low rank reconstruction for Oscillate MRE
    G McIlvain, AM Cerjanic, AG Christodoulou, MDJ McGarry, CL Johnson, “OSCILLATE: A Low-Rank Approach for Accelerated Magnetic Resonance Elastography,” to be presented at the 28th Annual Meeting of the International Society for Magnetic Resonance in Medicine, Sydney, Australia, April 18-23, 2020.

## [v1.0.1] 2020-03-18
## Updating with tagged docker hub builds
* Implements pcSENSE reconstructions with both DFT and NUFFT based field correction
* Implements 2D/3D field corrected non-Cartesian reconstructions

## [v0.1.0] 2016-05-08
### Initial release
