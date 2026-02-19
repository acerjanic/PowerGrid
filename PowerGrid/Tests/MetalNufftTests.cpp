/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [MetalNufftTests.cpp]

    Synopsis    [Unit tests for Apple Metal GPU compute backend: vDSP FFT
                    wrapper (fftAccelerate) and Metal gridding bridge
                    (MetalGridding).]

    Description [All tests are gated on #ifdef METAL_COMPUTE and are intended
                    for native macOS / Apple Silicon builds only.  They are
                    no-ops when compiled without METAL_COMPUTE defined.]

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [2026-02-19]

 *****************************************************************************/

#include "catch.hpp"

#ifdef METAL_COMPUTE

#include "../PGIncludes.h"
#include "../Metal/fftAccelerate.h"
#include "../Metal/MetalGridding.h"
#include "../Gnufft.h"
#include "../griddingSupport.h"

#include <cmath>
#include <vector>

using namespace arma;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

/// Fill an interleaved-complex float buffer with a real 2-D Gaussian blob.
/// Layout: buf[2*(iy + ix*Ny)] = real, buf[2*(iy + ix*Ny)+1] = imag.
static void fillGaussian2D(float *buf, int Nx, int Ny) {
    for (int ix = 0; ix < Nx; ix++) {
        for (int iy = 0; iy < Ny; iy++) {
            float fx = ((float)ix - Nx * 0.5f) / (Nx * 0.15f);
            float fy = ((float)iy - Ny * 0.5f) / (Ny * 0.15f);
            int idx = 2 * (iy + ix * Ny);
            buf[idx]     = std::exp(-(fx * fx + fy * fy));
            buf[idx + 1] = 0.0f;
        }
    }
}

// ---------------------------------------------------------------------------
// vDSP FFT tests
// ---------------------------------------------------------------------------

TEST_CASE("vDSP FFT 2D float: round-trip identity", "[fftAccelerate]") {
    // fft2dAccelerate and ifft2dAccelerate match FFTW convention:
    //   forward  = unnormalized DFT   (FFTW_FORWARD)
    //   backward = unnormalized IDFT  (FFTW_BACKWARD = N * normalized IFFT)
    // Therefore  IFFT(FFT(x)) = N * x, not x.
    // Verify by scaling the result back by 1/N before comparison.
    const int Nx = 64, Ny = 64;
    const float invN = 1.0f / (float)(Nx * Ny);
    std::vector<float> data(2 * Nx * Ny);
    fillGaussian2D(data.data(), Nx, Ny);
    const std::vector<float> orig = data;

    fft2dAccelerate(data.data(), (uword)Nx, (uword)Ny);
    ifft2dAccelerate(data.data(), (uword)Nx, (uword)Ny);

    // Scale down by N to recover x
    for (float& v : data) v *= invN;

    float maxVal = 0.0f, maxErr = 0.0f;
    for (size_t i = 0; i < orig.size(); i++)
        maxVal = std::max(maxVal, std::abs(orig[i]));
    for (size_t i = 0; i < orig.size(); i++)
        maxErr = std::max(maxErr, std::abs(data[i] - orig[i]));

    float relErr = maxErr / (maxVal + 1e-12f);
    INFO("max relative error = " << relErr);
    REQUIRE(relErr < 1e-5f);
}

TEST_CASE("vDSP FFT 2D float: linearity F(a+b) == F(a)+F(b)", "[fftAccelerate]") {
    const int Nx = 64, Ny = 64, N = 2 * Nx * Ny;
    std::vector<float> a(N), b(N), apb(N);
    fillGaussian2D(a.data(), Nx, Ny);
    // b is a spatially shifted copy (reverse ordering)
    for (int i = 0; i < N; i++) b[i] = a[N - 1 - i];
    for (int i = 0; i < N; i++) apb[i] = a[i] + b[i];

    std::vector<float> Fa = a, Fb = b, Fapb = apb;
    fft2dAccelerate(Fa.data(),   (uword)Nx, (uword)Ny);
    fft2dAccelerate(Fb.data(),   (uword)Nx, (uword)Ny);
    fft2dAccelerate(Fapb.data(), (uword)Nx, (uword)Ny);

    float maxErr = 0.0f;
    for (int i = 0; i < N; i++)
        maxErr = std::max(maxErr, std::abs(Fapb[i] - Fa[i] - Fb[i]));

    INFO("max linearity error = " << maxErr);
    REQUIRE(maxErr < 1e-4f);
}

TEST_CASE("vDSP FFT 2D float: Parseval's theorem", "[fftAccelerate]") {
    const int Nx = 64, Ny = 64, N = Nx * Ny;
    std::vector<float> data(2 * N);
    fillGaussian2D(data.data(), Nx, Ny);

    // Spatial-domain L2 norm squared
    double spatialPower = 0.0;
    for (int i = 0; i < N; i++)
        spatialPower += (double)data[2*i] * data[2*i] + (double)data[2*i+1] * data[2*i+1];

    // FFT
    std::vector<float> fdata = data;
    fft2dAccelerate(fdata.data(), (uword)Nx, (uword)Ny);

    // Spectral L2 norm squared (should equal N * spatialPower by Parseval)
    double spectralPower = 0.0;
    for (int i = 0; i < N; i++)
        spectralPower += (double)fdata[2*i] * fdata[2*i] + (double)fdata[2*i+1] * fdata[2*i+1];

    double relErr = std::abs(spectralPower - (double)N * spatialPower) /
                   ((double)N * spatialPower + 1e-12);
    INFO("Parseval relative error = " << relErr);
    REQUIRE(relErr < 1e-4);
}

// ---------------------------------------------------------------------------
// Metal gridding bridge tests
// ---------------------------------------------------------------------------

TEST_CASE("Metal gridding: context create/destroy does not crash", "[MetalGridding]") {
    float gridOS = 2.0f, kernelWidth = 4.0f;
    float beta = (float)(MRI_PI * std::sqrt(
        (double)(gridOS - 0.5f) * (gridOS - 0.5f) *
        (double)(kernelWidth * kernelWidth * 4.0f) / (double)(gridOS * gridOS) - 0.8));

    float *LUT = nullptr;
    uword sizeLUT = 0;
    calculateLUT(beta, kernelWidth, LUT, sizeLUT);
    REQUIRE(LUT    != nullptr);
    REQUIRE(sizeLUT > 0);

    const int nSamples = 4;
    float kx[] = { -16.0f, 16.0f, -16.0f, 16.0f };
    float ky[] = { -16.0f, -16.0f, 16.0f,  16.0f };
    float kz[] = {   0.0f,   0.0f,  0.0f,   0.0f };

    MetalGriddingContext *ctx = metal_gridding_create(
        128, 128, 1,   // gridNx, gridNy, gridNz
         64,  64, 1,   // imageNx, imageNy, imageNz
        gridOS, kernelWidth, LUT, (int)sizeLUT,
        kx, ky, kz, nSamples);

    // On Apple M1+ this succeeds; returns nullptr on Intel or non-Apple6 GPU.
    if (ctx) {
        metal_gridding_destroy(ctx);
    }
    free(LUT);
    SUCCEED();   // no crash → pass
}

TEST_CASE("Metal gridding: adjoint scatters energy onto grid", "[MetalGridding]") {
    float gridOS = 2.0f, kernelWidth = 4.0f;
    float beta = (float)(MRI_PI * std::sqrt(
        (double)(gridOS - 0.5f) * (gridOS - 0.5f) *
        (double)(kernelWidth * kernelWidth * 4.0f) / (double)(gridOS * gridOS) - 0.8));

    float *LUT = nullptr;
    uword sizeLUT = 0;
    calculateLUT(beta, kernelWidth, LUT, sizeLUT);

    // Single DC sample (kx=ky=0)
    const int nSamples = 1;
    float kx[] = { 0.0f }, ky[] = { 0.0f }, kz[] = { 0.0f };

    MetalGriddingContext *ctx = metal_gridding_create(
        128, 128, 1, 64, 64, 1,
        gridOS, kernelWidth, LUT, (int)sizeLUT,
        kx, ky, kz, nSamples);

    if (!ctx) { free(LUT); return; }   // non-Apple6 hardware: skip gracefully

    std::vector<float> dIn(2 * nSamples, 0.0f);
    dIn[0] = 1.0f;   // real = 1, imag = 0
    std::vector<float> gridOut(2 * 128 * 128, 0.0f);

    metal_gridding_adjoint_2D(ctx, dIn.data(), gridOut.data());

    float energy = 0.0f;
    for (float v : gridOut) energy += v * v;
    REQUIRE(energy > 0.0f);   // scatter must deposit something

    metal_gridding_destroy(ctx);
    free(LUT);
}

TEST_CASE("Metal gridding: adjoint twice does not accumulate (zero-init)", "[MetalGridding]") {
    float gridOS = 2.0f, kernelWidth = 4.0f;
    float beta = (float)(MRI_PI * std::sqrt(
        (double)(gridOS - 0.5f) * (gridOS - 0.5f) *
        (double)(kernelWidth * kernelWidth * 4.0f) / (double)(gridOS * gridOS) - 0.8));

    float *LUT = nullptr;
    uword sizeLUT = 0;
    calculateLUT(beta, kernelWidth, LUT, sizeLUT);

    const int nSamples = 1;
    float kx[] = { 0.0f }, ky[] = { 0.0f }, kz[] = { 0.0f };

    MetalGriddingContext *ctx = metal_gridding_create(
        128, 128, 1, 64, 64, 1,
        gridOS, kernelWidth, LUT, (int)sizeLUT,
        kx, ky, kz, nSamples);

    if (!ctx) { free(LUT); return; }

    std::vector<float> dIn(2, 0.0f);
    dIn[0] = 1.0f;
    std::vector<float> grid1(2 * 128 * 128, 0.0f);
    std::vector<float> grid2(2 * 128 * 128, 0.0f);

    // Two separate adjoint calls — each should produce the same result (zero-init)
    metal_gridding_adjoint_2D(ctx, dIn.data(), grid1.data());
    metal_gridding_adjoint_2D(ctx, dIn.data(), grid2.data());

    // grid1 and grid2 should be identical (not doubled)
    float maxDiff = 0.0f;
    for (size_t i = 0; i < grid1.size(); i++)
        maxDiff = std::max(maxDiff, std::abs(grid1[i] - grid2[i]));

    INFO("max diff between two identical adjoint calls = " << maxDiff);
    REQUIRE(maxDiff < 1e-6f);

    metal_gridding_destroy(ctx);
    free(LUT);
}

// ---------------------------------------------------------------------------
// Gnufft<float> adjoint-consistency test (Metal backend)
// ---------------------------------------------------------------------------

TEST_CASE("Gnufft<float> Metal: adjoint consistency <Ax,y> = <x,A'y>",
          "[Gnufft][Metal]") {
    const uword Nx = 32, Ny = 32;
    const uword nSamples = 64;
    const float gridOS = 2.0f;

    // Simple radial-like trajectory inside [-Nx/2, Nx/2]
    Col<float> kx(nSamples), ky(nSamples), kz(nSamples);
    kz.zeros();
    for (uword i = 0; i < nSamples; i++) {
        float angle = (float)i / (float)nSamples * (float)MRI_PI;
        float r     = (float)(i % (Nx / 2));
        kx(i) = r * std::cos(angle);
        ky(i) = r * std::sin(angle);
    }

    // Image-space coords (stored but not used in Metal transform path)
    Col<float> ix(Nx * Ny), iy(Nx * Ny), iz(Nx * Ny);
    iz.zeros();
    for (uword row = 0; row < Ny; row++) {
        for (uword col = 0; col < Nx; col++) {
            ix(col + row * Nx) = (float)col / (float)Nx - 0.5f;
            iy(col + row * Nx) = (float)row / (float)Ny - 0.5f;
        }
    }

    Gnufft<float> G(nSamples, gridOS, Nx, Ny, 1, kx, ky, kz, ix, iy, iz);

    // Deterministic complex test vectors (avoid arma_rng for reproducibility)
    Col<cx_float> x(Nx * Ny), y(nSamples);
    for (uword i = 0; i < Nx * Ny; i++)
        x(i) = cx_float(std::cos((float)i * 0.1f), std::sin((float)i * 0.07f));
    for (uword i = 0; i < nSamples; i++)
        y(i) = cx_float(std::cos((float)i * 0.13f), std::sin((float)i * 0.17f));

    Col<cx_float> Ax  = G * x;
    Col<cx_float> AHy = G / y;

    // cdot(a, b) = conj(a)^T * b = inner product
    cx_float lhs = cdot(Ax, y);    // <Ax, y>
    cx_float rhs = cdot(x, AHy);   // <x, A^H y>

    float relErr = std::abs(lhs - rhs) / (std::abs(lhs) + 1e-10f);
    INFO("<Ax,y> = " << lhs);
    INFO("<x,A'y> = " << rhs);
    INFO("relative error = " << relErr);
    REQUIRE(relErr < 5e-3f);   // 0.5% tolerance for float NUFFT
}

#endif // METAL_COMPUTE
