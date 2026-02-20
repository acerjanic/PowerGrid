/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file MetalBenchmarks.cpp
/// @brief Standalone throughput benchmark: Metal/vDSP vs CPU/FFTW.
///
/// Measures median latency (ms) per call for:
///   1. 2D FFT (vDSP vs FFTW, float, several grid sizes)
///   2. End-to-end Gnufft<float> Metal path vs Gnufft<double> CPU path
///      at several (image-size, nSamples) combinations.
///
/// Build with METAL_COMPUTE=ON, then run:
///   ./metal_bench
///
/// Note: Gnufft<double> uses the CPU/FFTW path even when METAL_COMPUTE is
/// defined because the Metal path is float-only (metalCtx stays nullptr).

#ifdef METAL_COMPUTE

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

#include "PGIncludes.h"
#include "Gnufft.h"
#include "Metal/fftAccelerate.h"
#include "fftCPU.h"
#include "pgCol.hpp"
#include "pgComplex.hpp"

using namespace arma;
using Clock = std::chrono::high_resolution_clock;
using Ms    = std::chrono::duration<double, std::milli>;

// Prevent compiler from optimizing away a computed result.
template<typename T>
static void doNotOptimize(const T& val) {
    asm volatile("" : : "r,m"(&val) : "memory");
}

// ---------------------------------------------------------------------------
// Timing helper: warm up, then collect nruns samples and return the median.
// ---------------------------------------------------------------------------
template <typename Fn>
static double median_ms(Fn fn, int nwarmup = 3, int nruns = 15) {
    for (int i = 0; i < nwarmup; i++) fn();
    std::vector<double> t(nruns);
    for (int i = 0; i < nruns; i++) {
        auto t0 = Clock::now();
        fn();
        t[i] = Ms(Clock::now() - t0).count();
    }
    std::sort(t.begin(), t.end());
    return t[nruns / 2];
}

// ---------------------------------------------------------------------------
// Build a simple radial-like 2D trajectory within [-Nx/2 .. Nx/2].
// ---------------------------------------------------------------------------
static void makeTrajectory(uword nS, uword Nx,
                           Col<float>& kx, Col<float>& ky, Col<float>& kz) {
    kx.set_size(nS); ky.set_size(nS); kz.zeros(nS);
    const uword nSpokes = 8;              // number of radial spokes
    const float maxR    = (float)(Nx / 2) * 0.95f;
    for (uword i = 0; i < nS; i++) {
        float angle = (float)(i % nSpokes) / (float)nSpokes * MRI_PI;
        float r     = maxR * (float)(i / nSpokes) / (float)(nS / nSpokes + 1);
        kx(i) = r * std::cos(angle);
        ky(i) = r * std::sin(angle);
    }
}

// ---------------------------------------------------------------------------
int main()
// ---------------------------------------------------------------------------
{
    printf("=== PowerGrid Metal/vDSP vs CPU/FFTW throughput (median of 15 runs) ===\n\n");

    // -----------------------------------------------------------------------
    // Part 1: 2D FFT comparison — vDSP vs FFTW, float, power-of-2 sizes.
    // -----------------------------------------------------------------------
    printf("[ 2D FFT — float, in-place, forward only ]\n");
    printf("  %-10s  %10s  %10s  %8s\n", "Grid size", "vDSP (ms)", "FFTW (ms)", "Speedup");

    for (int sz : {64, 128, 256, 512}) {
        const int N = sz * sz;
        // Use separate buffers so each run operates on realistic (non-zero) data.
        std::vector<float> vbuf(2 * N, 1.0f);
        std::vector<float> fbuf(2 * N, 1.0f);

        double tv = median_ms([&]{ fft2dAccelerate(vbuf.data(), (uword)sz, (uword)sz); });
        double tf = median_ms([&]{ fft2dCPU<float>(fbuf.data(), (uword)sz, (uword)sz); });

        printf("  %4dx%-4d    %10.3f  %10.3f  %7.2fx\n", sz, sz, tv, tf,
               (tf > 0.0) ? tf / tv : 0.0);
    }
    printf("\n");

    // -----------------------------------------------------------------------
    // Part 2: End-to-end Gnufft — Metal<float> vs CPU<double>.
    //
    // Two image sizes are tested so the FFT cost (proportional to N log N) and
    // gridding cost (proportional to nSamples * kernelWidth^2) can be separated.
    // -----------------------------------------------------------------------
    const float gridOS = 2.0f;

    for (uword Nx : {128u, 256u}) {
        const uword Ny = Nx;

        // Image-space coordinates
        Col<float> ix(Nx*Ny), iy(Nx*Ny), iz(Nx*Ny);
        iz.zeros();
        for (uword r = 0; r < Ny; r++)
            for (uword c = 0; c < Nx; c++) {
                ix(c + r*Nx) = (float)c / (float)Nx - 0.5f;
                iy(c + r*Nx) = (float)r / (float)Ny - 0.5f;
            }
        Col<double> ixd = conv_to<Col<double>>::from(ix);
        Col<double> iyd = conv_to<Col<double>>::from(iy);
        Col<double> izd(Nx*Ny, arma::fill::zeros);

        printf("[ Gnufft 2D  image %lux%lu  gridOS=%.1f ]\n",
               (unsigned long)Nx, (unsigned long)Ny, (double)gridOS);
        printf("  %-9s  %-8s  %12s  %12s  %8s\n",
               "nSamples", "op", "Metal (ms)", "CPU (ms)", "Speedup");

        for (uword nS : {1024u, 4096u, 16384u, 65536u}) {
            Col<float> kx, ky, kz;
            makeTrajectory(nS, Nx, kx, ky, kz);
            Col<double> kxd = conv_to<Col<double>>::from(kx);
            Col<double> kyd = conv_to<Col<double>>::from(ky);
            Col<double> kzd(nS, arma::fill::zeros);

            // Construct both operators (constructor timing not measured)
            Gnufft<float>  Gm(nS, gridOS, Nx, Ny, 1, kx, ky, kz, ix, iy, iz);
            Gnufft<double> Gc(nS, (double)gridOS, Nx, Ny, 1,
                               kxd, kyd, kzd, ixd, iyd, izd);

            // Fixed input vectors (const refs, not modified by the operators)
            Col<cx_float>  xf(Nx*Ny, arma::fill::randn);
            Col<cx_float>  yf(nS,    arma::fill::randn);
            Col<cx_double> xd = conv_to<Col<cx_double>>::from(xf);
            Col<cx_double> yd = conv_to<Col<cx_double>>::from(yf);

            // Fewer repetitions for large problems to keep total runtime reasonable
            int nr = (nS >= 16384) ? 8 : 15;

            // Forward: G * x  (image → k-space)
            double tFwdM = median_ms([&]{ (void)(Gm * xf); }, 2, nr);
            double tFwdC = median_ms([&]{ (void)(Gc * xd); }, 2, nr);
            printf("  %-9u  forward   %10.1f  %10.1f  %7.1fx\n",
                   (unsigned)nS, tFwdM, tFwdC,
                   (tFwdM > 0.0) ? tFwdC / tFwdM : 0.0);

            // Adjoint: G / y  (k-space → image)
            double tAdjM = median_ms([&]{ (void)(Gm / yf); }, 2, nr);
            double tAdjC = median_ms([&]{ (void)(Gc / yd); }, 2, nr);
            printf("  %-9u  adjoint   %10.1f  %10.1f  %7.1fx\n",
                   (unsigned)nS, tAdjM, tAdjC,
                   (tAdjM > 0.0) ? tAdjC / tAdjM : 0.0);
        }
        printf("\n");
    }

    // -----------------------------------------------------------------------
    // Part 3: pgCol vector algebra — Metal GPU vs CPU scalar loops.
    //
    // Tests key operations from the PCG solver and per-coil pipelines.
    // For each operation, the Metal threshold is 4096 elements, so we
    // test sizes above that where Metal dispatch is active.
    // -----------------------------------------------------------------------
    printf("[ pgCol vector algebra — Metal GPU vs CPU ]\n");
    printf("  %-10s  %-24s  %10s  %10s  %8s\n",
           "N", "Operation", "Metal (ms)", "CPU (ms)", "Speedup");

    // Helper: disable Metal temporarily by using double (which always falls
    // through to CPU) and compare against float (which uses Metal).
    // We'll time the pgCol operators directly.

    for (uword N : {8192u, 65536u, 262144u, 1048576u}) {
        // ----- Real element-wise: float add -----
        {
            pgCol<float> Af(N), Bf(N);
            for (uword i = 0; i < N; i++) {
                Af.at(i) = (float)(i % 1000) * 0.001f;
                Bf.at(i) = (float)((i + 500) % 1000) * 0.001f;
            }
            // Metal path (float, N >= 4096)
            double tMetal = median_ms([&]{ auto r = Af + Bf; doNotOptimize(r); }, 3, 15);
            // CPU path (double)
            pgCol<double> Ad(N), Bd(N);
            for (uword i = 0; i < N; i++) {
                Ad.at(i) = (double)Af.at(i);
                Bd.at(i) = (double)Bf.at(i);
            }
            double tCPU = median_ms([&]{ auto r = Ad + Bd; doNotOptimize(r); }, 3, 15);
            printf("  %-10lu  %-24s  %10.3f  %10.3f  %7.2fx\n",
                   (unsigned long)N, "real add (float vs dbl)",
                   tMetal, tCPU, (tMetal > 0.0) ? tCPU / tMetal : 0.0);
        }

        // ----- Complex element-wise multiply (Hadamard) -----
        {
            pgCol<pgComplex<float>> Af(N), Bf(N);
            for (uword i = 0; i < N; i++) {
                float re = (float)(i % 1000) * 0.001f - 0.5f;
                float im = (float)((i + 333) % 1000) * 0.001f - 0.5f;
                Af.at(i) = pgComplex<float>(re, im);
                Bf.at(i) = pgComplex<float>(im, re);
            }
            double tMetal = median_ms([&]{ auto r = Af % Bf; doNotOptimize(r); }, 3, 15);
            // CPU via double complex
            pgCol<pgComplex<double>> Ad(N), Bd(N);
            for (uword i = 0; i < N; i++) {
                Ad.at(i) = pgComplex<double>(Af.at(i).real(), Af.at(i).imag());
                Bd.at(i) = pgComplex<double>(Bf.at(i).real(), Bf.at(i).imag());
            }
            double tCPU = median_ms([&]{ auto r = Ad % Bd; doNotOptimize(r); }, 3, 15);
            printf("  %-10lu  %-24s  %10.3f  %10.3f  %7.2fx\n",
                   (unsigned long)N, "cplx mul (float vs dbl)",
                   tMetal, tCPU, (tMetal > 0.0) ? tCPU / tMetal : 0.0);
        }

        // ----- Complex sum (reduction) -----
        {
            pgCol<float> Af(N);
            for (uword i = 0; i < N; i++) {
                Af.at(i) = (float)(i % 1000) * 0.001f - 0.5f;
            }
            double tMetal = median_ms([&]{ float r = sum(Af); doNotOptimize(r); }, 3, 15);
            pgCol<double> Ad(N);
            for (uword i = 0; i < N; i++) Ad.at(i) = (double)Af.at(i);
            double tCPU = median_ms([&]{ double r = sum(Ad); doNotOptimize(r); }, 3, 15);
            printf("  %-10lu  %-24s  %10.3f  %10.3f  %7.2fx\n",
                   (unsigned long)N, "real sum (float vs dbl)",
                   tMetal, tCPU, (tMetal > 0.0) ? tCPU / tMetal : 0.0);
        }

        // ----- cdot (complex dot product) -----
        {
            pgCol<pgComplex<float>> Af(N), Bf(N);
            for (uword i = 0; i < N; i++) {
                float re = (float)(i % 1000) * 0.001f - 0.5f;
                float im = (float)((i + 333) % 1000) * 0.001f - 0.5f;
                Af.at(i) = pgComplex<float>(re, im);
                Bf.at(i) = pgComplex<float>(im, re);
            }
            double tMetal = median_ms([&]{ auto r = cdot(Af, Bf); doNotOptimize(r); }, 3, 15);
            // CPU reference: compute cdot manually for double
            pgCol<pgComplex<double>> Ad(N), Bd(N);
            for (uword i = 0; i < N; i++) {
                Ad.at(i) = pgComplex<double>(Af.at(i).real(), Af.at(i).imag());
                Bd.at(i) = pgComplex<double>(Bf.at(i).real(), Bf.at(i).imag());
            }
            double tCPU = median_ms([&]{ auto r = cdot(Ad, Bd); doNotOptimize(r); }, 3, 15);
            printf("  %-10lu  %-24s  %10.3f  %10.3f  %7.2fx\n",
                   (unsigned long)N, "cdot (float vs dbl)",
                   tMetal, tCPU, (tMetal > 0.0) ? tCPU / tMetal : 0.0);
        }

        printf("\n");
    }

    printf("Done.\n");
    return 0;
}

#else // METAL_COMPUTE not defined

#include <cstdio>
int main() {
    puts("METAL_COMPUTE is not enabled. Rebuild with -DMETAL_COMPUTE=ON.");
    return 1;
}

#endif // METAL_COMPUTE
