# PowerGrid

PowerGrid is a C++ library for iterative MRI image reconstruction, developed by the
MRFIL Research Groups at the University of Illinois, Urbana-Champaign.

## Operator Convention

All encoding operators in PowerGrid follow a unified interface:

| Expression | Meaning |
|------------|---------|
| `G * x`    | Forward transform: image → k-space |
| `G / y`    | Adjoint transform: k-space → image |

Both operators accept and return `arma::Col<std::complex<T1>>` vectors.

## Core Classes

### Encoding Operators

| Class | Description |
|-------|-------------|
| `Gdft<T1>` | Direct (brute-force) non-uniform DFT with optional off-resonance correction |
| `GdftR2<T1>` | Gdft extended with R2\* field map gradient support |
| `Gnufft<T1>` | Fast NUFFT via Kaiser-Bessel gridding (CPU + GPU via OpenACC) |
| `Gfft<T1>` | Uniform Cartesian FFT (wraps FFTW / cuFFT) |
| `SENSE<T1, Tobj>` | Multi-coil sensitivity encoding wrapping any encoding operator |
| `TimeSegmentation<T1, Tobj>` | Off-resonance correction via time-segmented interpolation |

### Regularization

| Class | Description |
|-------|-------------|
| `Robject<T1>` | Abstract base class for penalty functions |
| `QuadPenalty<T1>` | Quadratic (Tikhonov) penalty: β · ½ ‖C·x‖² |
| `TVPenalty<T1>` | Approximate total-variation penalty |

### Solvers

| Function | Description |
|----------|-------------|
| `solve_pwls_pcg` | Preconditioned weighted least-squares via conjugate gradient |
| `reconSolve` | High-level reconstruction entry point |

## Quick Example

```cpp
#include "PowerGrid/Gdft.h"
#include "PowerGrid/QuadPenalty.h"
#include "PowerGrid/solve_pwls_pcg.hpp"
using namespace arma;

// Build a Cartesian Gdft operator (8×8 image, 64 k-space samples)
uword Nx = 8, Ny = 8, n1 = Nx*Ny, n2 = Nx*Ny;
Col<float> kx(n2), ky(n2), kz(n2, fill::zeros);
Col<float> ix(n1), iy(n1), iz(n1, fill::zeros);
// ... fill kx/ky with Cartesian grid, ix/iy with pixel positions ...
Col<float> FM(n1, fill::zeros), t(n2, fill::zeros);
Gdft<float> G(n1, n2, kx, ky, kz, ix, iy, iz, FM, t);

// Forward transform
Col<cx_float> image(n1, fill::ones);
Col<cx_float> kspace = G * image;

// Adjoint transform
Col<cx_float> back = G / kspace;

// Iterative reconstruction via PCG
Col<float> W(n2, fill::ones);
QuadPenalty<float> R(Nx, Ny, 1, 0.01f);
Col<cx_float> x0(n1, fill::zeros);
Col<cx_float> xhat = solve_pwls_pcg<float>(x0, G, W, kspace, R, 50);
```

## License

University of Illinois/NCSA Open Source License.
See [LICENSE](https://github.com/acerjanic/PowerGrid/blob/HPCSDK_25.1/LICENSE) for details.
