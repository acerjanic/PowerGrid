# pgCol / pgMat Developer Guide

## Overview

`pgCol<T>` and `pgMat<T>` are PowerGrid's portable vector and matrix types. They support:
- CPU scalar loops (all platforms)
- OpenACC GPU offloading (NVIDIA, via nvc++)
- Apple Metal GPU dispatch (float only, Apple Silicon)

## Key Conventions

### pgMat set_size

**IMPORTANT:** `pgMat::set_size(nCols, nRows)` - first argument is columns, second is rows. This is opposite to most matrix APIs. The constructor `pgMat(nRows, nCols)` takes rows first (normal order), but internally calls `set_size(nCols, nRows)`.

### Column-Major Layout

Data is stored column-major: element `(r, c)` is at index `r + n_rows * c`. This matches Armadillo/LAPACK/Fortran conventions.

### Type System

| pgCol type | Scalar type | Metal dispatch? |
|------------|-------------|-----------------|
| `pgCol<float>` | `float` | Yes (N >= 4096) |
| `pgCol<double>` | `double` | No (CPU only) |
| `pgCol<pgComplex<float>>` | `pgComplex<float>` | Yes (N >= 4096) |
| `pgCol<pgComplex<double>>` | `pgComplex<double>` | No (CPU only) |

## View Mode

A pgCol view wraps external memory without owning it. Views enable zero-copy column extraction from pgMat and subvec slicing.

### Creating Views

```cpp
// Factory method
pgCol<float> v = pgCol<float>::view(ptr, length);

// From pgMat column
pgCol<float> col = mat.col(2);      // view into matrix column 2
col.is_view();                       // true

// Subvec
pgCol<float> sv = vec.subvec(10, 19); // view of elements 10..19
```

### View Semantics

- **Write-through:** Modifications via a view write to the underlying memory.
  ```cpp
  mat.col(1) %= scale;  // modifies matrix column 1 in-place
  ```
- **Copy produces owning:** `pgCol<float> copy(view)` deep-copies; `copy.is_view() == false`.
- **Destructor skips free:** Views don't free memory on destruction.
- **Lifetime:** The view must not outlive the underlying memory.
- **set_size breaks view:** Calling `set_size()` on a view breaks the view (sets `isView_ = false`) and allocates fresh memory. The original memory is NOT freed.

### col_copy vs col

- `mat.col(ii)` returns a non-owning view (fast, no allocation)
- `mat.col_copy(ii)` returns a deep copy (safe for long-lived use)
- `mat.set_col(ii, src)` writes a pgCol into column ii (memcpy)

## Arma <-> pgCol Conversion

```cpp
// arma -> pgCol (copies data)
arma::Col<std::complex<float>> arma_vec = ...;
pgCol<pgComplex<float>> pg_vec(arma_vec);

// pgCol -> arma (copies data)
arma::Col<std::complex<float>> back = pg_vec.getArma();
```

Both conversions are memcpy. `pgComplex<T>` and `std::complex<T>` have identical memory layout.

### Pattern for Gnufft/Robject Boundaries

```cpp
// Forward: pgCol -> arma -> Gnufft -> arma -> pgCol
Col<CxT1> result = (*G_obj) * weighted_pg.getArma();
outData_pg.set_col(ii, pgCol<pgComplex<T1>>(result));

// Adjoint: same pattern
Col<CxT1> adjResult = (*G_obj) / slice_pg.getArma();
pgCol<pgComplex<T1>> adj_pg(adjResult);
```

## Available Operators

### Element-wise (pgCol)

| Operator | Description | Metal dispatch |
|----------|-------------|----------------|
| `A + B` | addition | `cvec_add` / `vec_add` |
| `A - B` | subtraction | `cvec_sub` / `vec_sub` |
| `A % B` | element-wise multiply | `cvec_mul` / `vec_mul` |
| `A / B` | element-wise divide | `cvec_div` / `vec_div` |
| `A + s` | scalar add | `vec_add_scalar` |
| `A % s` | scalar multiply | `cvec_mul_scalar` / `vec_mul_scalar` |
| `W % X` | real * complex | `rvec_cmul` |
| `A += B` | in-place add | same kernels |
| `A %= B` | in-place multiply | same kernels |

### Free Functions

| Function | Description | Metal dispatch |
|----------|-------------|----------------|
| `sum(A)` | sum of elements | `vec_sum` |
| `cdot(A, B)` | complex dot product | `cvec_cdot` |
| `norm(A)` | L2 norm | `cvec_norm2sq` + sqrt |
| `conj(A)` | element-wise conjugate | CPU only |
| `abs(A)` | element-wise magnitude | CPU only |
| `real(A)` | extract real parts | CPU only |
| `imag(A)` | extract imaginary parts | CPU only |

### pgMat Free Functions

| Function | Description |
|----------|-------------|
| `sum(M, dim)` | row/column sum (dim=0: col-wise, dim=1: row-wise) |
| `vectorise(M)` | flatten to pgCol |
| `conj(M)` | element-wise conjugate |

## Aligned Allocation (METAL_COMPUTE)

When `METAL_COMPUTE` is defined, `set_size()` uses `std::aligned_alloc(16384, ...)` for page-aligned allocation (Apple Silicon page size = 16384). This enables Phase 3 `newBufferWithBytesNoCopy` for true zero-copy GPU buffers.

The destructor uses `std::free()` instead of `delete[]` accordingly.

## has_nan()

Detects NaN values using bit-level inspection (works under `-ffast-math`):

```cpp
pgCol<float> v(100);
v.ones();
// ... some computation ...
if (v.has_nan()) { /* handle error */ }
```

Works for `float`, `double`, `pgComplex<float>`, and `pgComplex<double>`.
