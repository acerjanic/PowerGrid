/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file fftAccelerate.cpp
/// @brief vDSP/Accelerate in-place FFT for the Metal compute path.
///
/// Float overloads use vDSP_fft2d_zip / vDSP_fft3d_zop with split-complex
/// representation.  Dimensions that are not powers of two fall back to the
/// existing FFTW wrappers in fftCPU.h.
///
/// Double overloads always delegate to FFTW because Metal shaders only support
/// single precision.

#ifdef METAL_COMPUTE

#include "fftAccelerate.h"
#include "../fftCPU.h"     // for fallback to FFTW

#include <Accelerate/Accelerate.h>
#include <cmath>
#include <cstdlib>
#include <cstring>
#include <map>
#include <mutex>

// ---------------------------------------------------------------------------
// FFTSetup cache (creation is expensive; reuse across calls for same log2N)
// ---------------------------------------------------------------------------
namespace {

// Key: (log2rows, log2cols, log2slices) — slices=0 for 2D
using SetupKey = std::tuple<vDSP_Length, vDSP_Length, vDSP_Length>;

struct SetupCache {
    std::map<SetupKey, FFTSetup> cache;
    std::mutex mtx;

    FFTSetup get(vDSP_Length log2r, vDSP_Length log2c, vDSP_Length log2s = 0) {
        SetupKey k{log2r, log2c, log2s};
        std::lock_guard<std::mutex> lk(mtx);
        auto it = cache.find(k);
        if (it != cache.end()) return it->second;
        vDSP_Length maxN = std::max({log2r, log2c, log2s});
        FFTSetup s = vDSP_create_fftsetup(maxN, FFT_RADIX2);
        cache[k] = s;
        return s;
    }
};

SetupCache& getCache() {
    static SetupCache c;
    return c;
}

// Returns true if n is a power of 2 and > 0
static inline bool isPow2(uword n) {
    return (n > 0) && ((n & (n - 1)) == 0);
}

static inline vDSP_Length ilog2(uword n) {
    vDSP_Length k = 0;
    while ((uword(1) << k) < n) ++k;
    return k;
}

} // anonymous namespace

// ---------------------------------------------------------------------------
// float 2D FFT / IFFT
// ---------------------------------------------------------------------------

void fft2dAccelerate(float* d_data, uword nx, uword ny) {
    if (!isPow2(nx) || !isPow2(ny)) {
        fft2dCPU(d_data, nx, ny);
        return;
    }

    const vDSP_Length log2nx = ilog2(nx);
    const vDSP_Length log2ny = ilog2(ny);
    const vDSP_Length N      = nx * ny;
    FFTSetup setup = getCache().get(log2nx, log2ny);

    // De-interleave: [re0, im0, re1, im1, ...] → split complex
    float* re = (float*)malloc(N * sizeof(float));
    float* im = (float*)malloc(N * sizeof(float));
    DSPSplitComplex split{ re, im };
    vDSP_ctoz((const DSPComplex*)d_data, 2, &split, 1, N);

    vDSP_fft2d_zip(setup, &split, 1, 0, log2nx, log2ny, kFFTDirection_Forward);

    vDSP_ztoc(&split, 1, (DSPComplex*)d_data, 2, N);
    free(re);
    free(im);
}

void ifft2dAccelerate(float* d_data, uword nx, uword ny) {
    if (!isPow2(nx) || !isPow2(ny)) {
        ifft2dCPU(d_data, nx, ny);
        return;
    }

    const vDSP_Length log2nx = ilog2(nx);
    const vDSP_Length log2ny = ilog2(ny);
    const vDSP_Length N      = nx * ny;
    FFTSetup setup = getCache().get(log2nx, log2ny);

    float* re = (float*)malloc(N * sizeof(float));
    float* im = (float*)malloc(N * sizeof(float));
    DSPSplitComplex split{ re, im };
    vDSP_ctoz((const DSPComplex*)d_data, 2, &split, 1, N);

    vDSP_fft2d_zip(setup, &split, 1, 0, log2nx, log2ny, kFFTDirection_Inverse);

    // Normalize by 1/(nx*ny)
    float scale = 1.0f / (float)(nx * ny);
    vDSP_vsmul(re, 1, &scale, re, 1, N);
    vDSP_vsmul(im, 1, &scale, im, 1, N);

    vDSP_ztoc(&split, 1, (DSPComplex*)d_data, 2, N);
    free(re);
    free(im);
}

// ---------------------------------------------------------------------------
// float 3D FFT / IFFT
// vDSP does not have a native 3D FFT routine; we decompose into 2D + 1D passes.
// ---------------------------------------------------------------------------

// Forward 3D: FFT all nz 2D planes, then FFT along z for each (x,y) column.
void fft3dAccelerate(float* d_data, uword nx, uword ny, uword nz) {
    if (!isPow2(nx) || !isPow2(ny) || !isPow2(nz)) {
        fft3dCPU(d_data, nx, ny, nz);
        return;
    }

    // 2D FFT each z-slice in-place
    const uword sliceFloats = nx * ny * 2;
    for (uword z = 0; z < nz; ++z) {
        fft2dAccelerate(d_data + z * sliceFloats, nx, ny);
    }

    // 1D FFT along z for each (x, y) pair
    // Gather a column, FFT, scatter back
    const vDSP_Length log2nz = ilog2(nz);
    FFTSetup setup = getCache().get(log2nz, log2nz);

    float* re = (float*)malloc(nz * sizeof(float));
    float* im = (float*)malloc(nz * sizeof(float));
    DSPSplitComplex split{ re, im };

    for (uword y = 0; y < ny; ++y) {
        for (uword x = 0; x < nx; ++x) {
            // Gather column: index in interleaved array = 2*(z*nx*ny + x*ny + y) or
            // matching gridding.cpp layout: idx = 2*(ny + x*gNy + z*gNx*gNy)
            // Generic column extraction:
            for (uword z = 0; z < nz; ++z) {
                uword idx = 2 * (y + x * ny + z * nx * ny);
                re[z] = d_data[idx];
                im[z] = d_data[idx + 1];
            }

            vDSP_fft_zip(setup, &split, 1, log2nz, kFFTDirection_Forward);

            for (uword z = 0; z < nz; ++z) {
                uword idx = 2 * (y + x * ny + z * nx * ny);
                d_data[idx]     = re[z];
                d_data[idx + 1] = im[z];
            }
        }
    }
    free(re);
    free(im);
}

void ifft3dAccelerate(float* d_data, uword nx, uword ny, uword nz) {
    if (!isPow2(nx) || !isPow2(ny) || !isPow2(nz)) {
        ifft3dCPU(d_data, nx, ny, nz);
        return;
    }

    // IFFT along z columns first
    const vDSP_Length log2nz = ilog2(nz);
    FFTSetup setup = getCache().get(log2nz, log2nz);
    float scaleZ = 1.0f / (float)nz;

    float* re = (float*)malloc(nz * sizeof(float));
    float* im = (float*)malloc(nz * sizeof(float));
    DSPSplitComplex split{ re, im };

    for (uword y = 0; y < ny; ++y) {
        for (uword x = 0; x < nx; ++x) {
            for (uword z = 0; z < nz; ++z) {
                uword idx = 2 * (y + x * ny + z * nx * ny);
                re[z] = d_data[idx];
                im[z] = d_data[idx + 1];
            }
            vDSP_fft_zip(setup, &split, 1, log2nz, kFFTDirection_Inverse);
            vDSP_vsmul(re, 1, &scaleZ, re, 1, nz);
            vDSP_vsmul(im, 1, &scaleZ, im, 1, nz);
            for (uword z = 0; z < nz; ++z) {
                uword idx = 2 * (y + x * ny + z * nx * ny);
                d_data[idx]     = re[z];
                d_data[idx + 1] = im[z];
            }
        }
    }
    free(re);
    free(im);

    // IFFT each 2D slice
    const uword sliceFloats = nx * ny * 2;
    for (uword z = 0; z < nz; ++z) {
        ifft2dAccelerate(d_data + z * sliceFloats, nx, ny);
    }
}

// ---------------------------------------------------------------------------
// double overloads — always delegate to FFTW
// ---------------------------------------------------------------------------

void fft2dAccelerate(double* d_data, uword nx, uword ny) {
    fft2dCPU(d_data, nx, ny);
}

void ifft2dAccelerate(double* d_data, uword nx, uword ny) {
    ifft2dCPU(d_data, nx, ny);
}

void fft3dAccelerate(double* d_data, uword nx, uword ny, uword nz) {
    fft3dCPU(d_data, nx, ny, nz);
}

void ifft3dAccelerate(double* d_data, uword nx, uword ny, uword nz) {
    ifft3dCPU(d_data, nx, ny, nz);
}

#endif // METAL_COMPUTE
