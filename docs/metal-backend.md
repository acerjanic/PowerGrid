# PowerGrid Metal GPU Backend

## Architecture Overview

The Metal backend accelerates PowerGrid's MRI reconstruction pipeline on Apple Silicon by dispatching vector algebra operations to the GPU. The implementation is structured in three phases:

### Phase 1: Metal Vector Algebra + pgCol Dispatch (Complete)

Added Metal compute shaders for element-wise vector operations and integrated dispatch into `pgCol` operators via `MetalVectorOps_dispatch.hpp`. Operations are dispatched to Metal when:

1. `METAL_COMPUTE` is defined at compile time
2. The element type is `float` or `pgComplex<float>` (double always falls back to CPU)
3. The vector size exceeds `kMinMetalSize` (4096 elements)

**Metal Kernels:** `vec_add`, `vec_sub`, `vec_mul`, `vec_div`, `vec_add_scalar`, `vec_mul_scalar`, `cvec_add`, `cvec_sub`, `cvec_mul`, `cvec_div`, `cvec_mul_scalar`, `cvec_axpy`, `rvec_cmul`, `vec_sum`, `cvec_cdot`, `cvec_norm2sq`

**Key files:**
- `PowerGrid/Metal/vectorops_metal.metal` - Metal shader source
- `PowerGrid/Metal/MetalVectorOps.h/.mm` - Objective-C++ wrapper
- `PowerGrid/Metal/MetalVectorOps_dispatch.hpp` - C++ inline dispatch helpers

### Phase 2: Type Migration to pgCol/pgMat (Complete)

Migrated hot inner loops of SENSE, pcSENSE, TimeSegmentation operators and the PCG solver from Armadillo types (`arma::Col<complex<T1>>`) to pgCol/pgMat types. This enables Metal dispatch to actually be exercised during reconstruction.

**Key changes:**
- `pgCol` view mode: non-owning wrapper for zero-copy column extraction and subvec
- `pgMat::col()` returns views instead of copies
- `aligned_alloc(16384)` for page-aligned allocation (enables Phase 3 zero-copy)
- Operators use `#ifdef METAL_COMPUTE` + `if constexpr (is_same<T1, float>)` pattern
- Gnufft/Robject interfaces remain Armadillo; conversion at call sites via `getArma()`/constructors

**Migrated files:**
- `SENSE.cpp` - Per-coil forward/adjoint
- `pcSENSE.cpp` - Per-coil-per-shot forward/adjoint
- `TimeSegmentation.cpp` - Per-segment forward/adjoint
- `solve_pwls_pcg.hpp` - PCG solver inner loop

### Phase 3: GPU-Resident Zero-Copy (Future)

Will use `newBufferWithBytesNoCopy` with page-aligned pgCol memory to eliminate CPU-GPU copies. The aligned allocation from Phase 2 is a prerequisite.

## Dispatch Chain

When `METAL_COMPUTE` is defined and `T1 = float`:

```
pgCol operator (e.g., +=, %, cdot)
  -> if constexpr (is_metal_type<T>)
    -> pg_metal::try_metal_*(args, n_elem)
      -> if n_elem >= kMinMetalSize && ctx() != nullptr
        -> metal_vec_*(ctx, args)    // GPU dispatch
      -> else return false           // CPU fallback
  -> CPU scalar loop fallback
```

For double types, `is_metal_type<double>` is `false_type`, so `if constexpr` eliminates the Metal code path entirely at compile time.

## Build Instructions

### Metal build (Apple Silicon Mac)

```bash
/opt/homebrew/bin/cmake -B build_native -S . \
  -DCMAKE_CXX_COMPILER=clang++ \
  -DOPENACC_GPU=OFF -DOPENACC_MP=OFF -DMPISupport=OFF \
  -DMETAL_COMPUTE=ON -DCMAKE_BUILD_TYPE=Debug

/opt/homebrew/bin/cmake --build build_native -j$(sysctl -n hw.logicalcpu)
```

### Docker build (Linux/GPU cluster, no Metal)

```bash
docker run --rm \
  -v $(pwd):/root/PowerGrid \
  ghcr.io/acerjanic/powergrid-hpcsdk:HPCSDK_25.1 \
  bash -c "cd /root/PowerGrid/PowerGrid/Tests && rm -rf build && mkdir -p build && cd build && cmake .. -DCMAKE_CXX_COMPILER=g++ && make -j4 && ./tests"
```

### Running tests

```bash
# All tests
./build_native/metal_tests

# Specific test categories
./build_native/metal_tests "[pgCol_view]"
./build_native/metal_tests "[pgMat_col_view]"
./build_native/metal_tests "[pgCol_has_nan]"
./build_native/metal_tests "[pgCol_rvec_cmul]"

# Benchmarks
./build_native/metal_bench
```
