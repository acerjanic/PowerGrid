/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file GdftBenchmarks.cpp
/// @brief Gdft and pcSENSE benchmarks: Metal GPU DFT vs CPU, with synthetic
///        linear ramp phase maps for phase-corrected SENSE reconstructions.

#include "catch.hpp"

#include "../Core/PGIncludes.h"
#include "../Operators/Gdft.h"
#include "../Operators/GdftR2.h"
#include "../Operators/pcSENSE.h"
#include "../Solvers/QuadPenalty.h"
#include "../Solvers/solve_pwls_pcg.hpp"

#include "Helpers/SyntheticPhantom.hpp"
#include "Helpers/SyntheticCoils.hpp"
#include "Helpers/SyntheticTrajectory.hpp"
#include "Helpers/SyntheticFieldMap.hpp"

using namespace arma;

// ---------------------------------------------------------------------------
// Helper: generate linear ramp shot phase maps
//
// Models the phase accumulated from constant gradient offsets (e.g. diffusion
// encoding, motion-induced phase). Each shot gets a different linear ramp
// direction uniformly spaced around the unit circle.
//
// PMap(x,y,shot_j) = maxPhase * (cos(angle_j)*x + sin(angle_j)*y)
// where angle_j = 2*pi*j/Ns and x,y are in [-1,1].
//
// Returns vectorized (Ni*Ns) column, where Ni = Nx*Ny.
// Phase values are in RADIANS.
// ---------------------------------------------------------------------------
template<typename T1>
Col<T1> linearRampPhaseMap2D(uword Nx, uword Ny, uword Ns, T1 maxPhase)
{
    uword Ni = Nx * Ny;
    Mat<T1> PMap(Ni, Ns);

    for (uword ss = 0; ss < Ns; ss++) {
        T1 angle = (T1)2.0 * M_PI * (T1)ss / (T1)Ns;
        T1 gx = maxPhase * std::cos(angle);
        T1 gy = maxPhase * std::sin(angle);

        for (uword jj = 0; jj < Ny; jj++) {
            T1 y = (T1)2.0 * ((T1)jj - (T1)Ny / (T1)2.0) / (T1)Ny;
            for (uword ii = 0; ii < Nx; ii++) {
                T1 x = (T1)2.0 * ((T1)ii - (T1)Nx / (T1)2.0) / (T1)Nx;
                PMap(ii + jj * Nx, ss) = gx * x + gy * y;
            }
        }
    }

    return vectorise(PMap);
}

// ---------------------------------------------------------------------------
// Helper: NRMSE
// ---------------------------------------------------------------------------
template<typename T1>
T1 nrmse_gdft(const Col<std::complex<T1>>& x, const Col<std::complex<T1>>& ref) {
    T1 errNorm = norm(x - ref, 2);
    T1 refNorm = norm(ref, 2);
    if (refNorm < (T1)1e-20) return errNorm;
    return errNorm / refNorm;
}

// ---------------------------------------------------------------------------
// Helper: image coordinates matching pcSENSE's Cube(Ny,Nx,Nz) layout
// (same as gnufftImageCoords2D in SyntheticReconTests.cpp)
// ---------------------------------------------------------------------------
template<typename T1>
void pcSENSEImageCoords2D(uword Nx, uword Ny,
                           Col<T1>& ix, Col<T1>& iy, Col<T1>& iz)
{
    uword Ni = Nx * Ny;
    ix.set_size(Ni);
    iy.set_size(Ni);
    iz.zeros(Ni);

    for (uword jj = 0; jj < Nx; jj++) {       // x
        for (uword ii = 0; ii < Ny; ii++) {    // y
            uword idx = ii + jj * Ny;
            ix(idx) = ((T1)jj - (T1)Nx / (T1)2.0) / (T1)Nx;
            iy(idx) = ((T1)ii - (T1)Ny / (T1)2.0) / (T1)Ny;
        }
    }
}

// =====================================================================
// Test: Gdft forward-adjoint adjointness check
// =====================================================================
TEST_CASE("Gdft adjointness", "[Gdft_synthetic][adjointness]") {
    typedef float T1;
    typedef std::complex<T1> CxT1;

    const uword Nx = 16;
    const uword Ny = 16;
    const uword Ni = Nx * Ny;
    const uword nSpokes = 32;
    const uword nReadout = 32;
    const uword Nd = nSpokes * nReadout;

    // Trajectory (radial, same as NUFFT tests)
    Col<T1> kx, ky, kz;
    radialTrajectory2D<T1>(nSpokes, nReadout, Nx, kx, ky, kz);

    // Image coordinates
    Col<T1> ix, iy, iz;
    pcSENSEImageCoords2D<T1>(Nx, Ny, ix, iy, iz);

    // Field map (small, ~20 Hz)
    Col<T1> fmapHz = syntheticFieldMap2D<T1>(Nx, Ny, (T1)20.0);
    Col<T1> fmapRad = fmapHz * (T1)(2.0 * M_PI);

    // Timing vector
    Col<T1> tvec = linearTimingVector<T1>(Nd, (T1)0.005);

    // Gdft operator
    Gdft<T1> G(Nd, Ni, kx, ky, kz, ix, iy, iz, fmapRad, tvec);

    // Random vectors
    arma_rng::set_seed(42);
    Col<CxT1> x = randn<Col<CxT1>>(Ni);
    Col<CxT1> y = randn<Col<CxT1>>(Nd);

    // <Gx, y> == <x, G'y>?
    Col<CxT1> Gx = G * x;
    Col<CxT1> Gty = G / y;

    CxT1 lhs = cdot(Gx, y);
    CxT1 rhs = cdot(x, Gty);

    T1 relErr = std::abs(lhs - rhs) / (std::abs(lhs) + (T1)1e-10);
    std::cout << "Gdft adjointness: relErr = " << relErr << std::endl;
    REQUIRE(relErr < (T1)1e-3);
}

// =====================================================================
// Test: pcSENSE adjointness with linear ramp phase maps
// =====================================================================
TEST_CASE("pcSENSE adjointness with linear ramp phase", "[pcSENSE_synthetic][adjointness]") {
    typedef float T1;
    typedef std::complex<T1> CxT1;

    const uword Nx = 16;
    const uword Ny = 16;
    const uword Ni = Nx * Ny;
    const uword Nc = 4;
    const uword Ns = 4;
    const uword nSpokes = 16;
    const uword nReadout = 32;
    const uword Nd = nSpokes * nReadout;

    // Trajectory: Ns copies of the same radial trajectory (one per shot)
    Col<T1> kx1, ky1, kz1;
    radialTrajectory2D<T1>(nSpokes, nReadout, Nx, kx1, ky1, kz1);
    Col<T1> kx = repmat(kx1, Ns, 1);
    Col<T1> ky = repmat(ky1, Ns, 1);
    Col<T1> kz = repmat(kz1, Ns, 1);

    // Timing vector: same for all shots
    Col<T1> tvec1 = linearTimingVector<T1>(Nd, (T1)0.005);
    Col<T1> tvec = repmat(tvec1, Ns, 1);

    // Coil maps
    Col<CxT1> SENSEmap = syntheticCoils2D<T1>(Nx, Ny, Nc);

    // Field map
    Col<T1> fmapHz = syntheticFieldMap2D<T1>(Nx, Ny, (T1)20.0);
    Col<T1> fmapRad = fmapHz * (T1)(2.0 * M_PI);

    // Linear ramp shot phase: maxPhase = pi/2 radians
    Col<T1> shotPhase = linearRampPhaseMap2D<T1>(Nx, Ny, Ns, (T1)(M_PI / 2.0));

    // pcSENSE operator
    pcSENSE<T1> P(kx, ky, kz, Nx, Ny, 1, Nc, tvec, SENSEmap, fmapRad, shotPhase);

    // Random vectors
    arma_rng::set_seed(77);
    Col<CxT1> x = randn<Col<CxT1>>(Ni);
    Col<CxT1> y = randn<Col<CxT1>>(Nd * Ns * Nc);

    // <Px, y> == <x, P'y>?
    Col<CxT1> Px = P * x;
    Col<CxT1> Pty = P / y;

    REQUIRE(Px.n_elem == Nd * Ns * Nc);
    REQUIRE(Pty.n_elem == Ni);

    CxT1 lhs = cdot(Px, y);
    CxT1 rhs = cdot(x, Pty);

    T1 relErr = std::abs(lhs - rhs) / (std::abs(lhs) + (T1)1e-10);
    std::cout << "pcSENSE adjointness (linear ramp phase): relErr = " << relErr << std::endl;
    REQUIRE(relErr < (T1)1e-2);
}

// =====================================================================
// Test: pcSENSE reconstruction with linear ramp phase
// =====================================================================
TEST_CASE("pcSENSE recon with linear ramp phase", "[pcSENSE_synthetic][recon]") {
    typedef float T1;
    typedef std::complex<T1> CxT1;

    const uword Nx = 32;
    const uword Ny = 32;
    const uword Ni = Nx * Ny;
    const uword Nc = 4;
    const uword Ns = 2;
    const uword nSpokes = 32;
    const uword nReadout = 32;
    const uword Nd = nSpokes * nReadout;
    const uword niter = 10;

    // Phantom
    Col<CxT1> phantom = sheppLogan2D<T1>(Nx, Ny);

    // Trajectory: same for all shots
    Col<T1> kx1, ky1, kz1;
    radialTrajectory2D<T1>(nSpokes, nReadout, Nx, kx1, ky1, kz1);
    Col<T1> kx = repmat(kx1, Ns, 1);
    Col<T1> ky = repmat(ky1, Ns, 1);
    Col<T1> kz = repmat(kz1, Ns, 1);

    // Timing vector
    Col<T1> tvec1 = linearTimingVector<T1>(Nd, (T1)0.005);
    Col<T1> tvec = repmat(tvec1, Ns, 1);

    // Coil maps
    Col<CxT1> SENSEmap = syntheticCoils2D<T1>(Nx, Ny, Nc);

    // Field map (mild, ~10 Hz)
    Col<T1> fmapHz = syntheticFieldMap2D<T1>(Nx, Ny, (T1)10.0);
    Col<T1> fmapRad = fmapHz * (T1)(2.0 * M_PI);

    // Linear ramp shot phase: pi/4 radians max
    Col<T1> shotPhase = linearRampPhaseMap2D<T1>(Nx, Ny, Ns, (T1)(M_PI / 4.0));

    // pcSENSE operator
    pcSENSE<T1> P(kx, ky, kz, Nx, Ny, 1, Nc, tvec, SENSEmap, fmapRad, shotPhase);

    // Generate data: y = P * phantom
    Col<CxT1> y = P * phantom;
    REQUIRE(y.n_elem == Nd * Ns * Nc);
    REQUIRE(!y.has_nan());

    // Weights and regularizer
    Col<T1> W = ones<Col<T1>>(Nd * Ns * Nc);
    QuadPenalty<T1> R(Nx, Ny, 1, (T1)1e-3, 2);
    Col<CxT1> x0 = zeros<Col<CxT1>>(Ni);

    // Reconstruct
    Col<CxT1> xhat = solve_pwls_pcg<T1>(x0, P, W, y, R, niter);
    REQUIRE(xhat.n_elem == Ni);
    REQUIRE(!xhat.has_nan());

    T1 err = nrmse_gdft<T1>(xhat, phantom);
    std::cout << "pcSENSE recon NRMSE (linear ramp phase) = " << err << std::endl;
    REQUIRE(err < (T1)0.60);
}

// =====================================================================
// Benchmark: Standalone Gdft forward/adjoint timing
// =====================================================================
TEST_CASE("Gdft forward/adjoint timing", "[Gdft_bench][benchmark]") {
    typedef float T1;
    typedef std::complex<T1> CxT1;

    struct GdftSize {
        uword Nx, Ny, nSpokes, nReadout;
        const char* label;
    };

    GdftSize sizes[] = {
        {16, 16, 32, 32, "16x16 / 1024 samples"},
        {32, 32, 64, 64, "32x32 / 4096 samples"},
        {64, 64, 64, 64, "64x64 / 4096 samples"},
        {128, 128, 32, 32, "128x128 / 1024 samples"},
        {256, 256, 32, 32, "256x256 / 1024 samples"},
    };

    std::cout << std::endl;
    std::cout << "Gdft Metal GPU DFT benchmark (forward + adjoint, median of 5):" << std::endl;
    std::cout << "  " << std::left << std::setw(28) << "Config"
              << std::setw(6) << "Ni"
              << std::setw(8) << "Nd"
              << std::setw(12) << "Fwd (ms)"
              << std::setw(12) << "Adj (ms)"
              << std::endl;

    for (const auto& sz : sizes) {
        uword Ni = sz.Nx * sz.Ny;
        uword Nd = sz.nSpokes * sz.nReadout;

        Col<T1> kx, ky, kz;
        radialTrajectory2D<T1>(sz.nSpokes, sz.nReadout, sz.Nx, kx, ky, kz);

        Col<T1> ix, iy, iz;
        pcSENSEImageCoords2D<T1>(sz.Nx, sz.Ny, ix, iy, iz);

        Col<T1> fmapRad = syntheticFieldMap2D<T1>(sz.Nx, sz.Ny, (T1)20.0) * (T1)(2.0 * M_PI);
        Col<T1> tvec = linearTimingVector<T1>(Nd, (T1)0.005);

        Gdft<T1> G(Nd, Ni, kx, ky, kz, ix, iy, iz, fmapRad, tvec);

        arma_rng::set_seed(42);
        Col<CxT1> x = randn<Col<CxT1>>(Ni);
        Col<CxT1> y = randn<Col<CxT1>>(Nd);

        // Warmup
        for (int i = 0; i < 2; i++) {
            (void)(G * x);
            (void)(G / y);
        }

        // Benchmark forward
        const int nruns = 5;
        std::vector<double> fwd_times(nruns), adj_times(nruns);

        for (int i = 0; i < nruns; i++) {
            auto t0 = std::chrono::high_resolution_clock::now();
            Col<CxT1> result = G * x;
            auto t1 = std::chrono::high_resolution_clock::now();
            fwd_times[i] = std::chrono::duration<double, std::milli>(t1 - t0).count();
        }

        for (int i = 0; i < nruns; i++) {
            auto t0 = std::chrono::high_resolution_clock::now();
            Col<CxT1> result = G / y;
            auto t1 = std::chrono::high_resolution_clock::now();
            adj_times[i] = std::chrono::duration<double, std::milli>(t1 - t0).count();
        }

        std::sort(fwd_times.begin(), fwd_times.end());
        std::sort(adj_times.begin(), adj_times.end());

        std::cout << "  " << std::left << std::setw(28) << sz.label
                  << std::setw(6) << Ni
                  << std::setw(8) << Nd
                  << std::setw(12) << std::fixed << std::setprecision(2) << fwd_times[nruns/2]
                  << std::setw(12) << adj_times[nruns/2]
                  << std::endl;

        // Sanity check: results should be finite
        REQUIRE(!Col<CxT1>(G * x).has_nan());
        REQUIRE(!Col<CxT1>(G / y).has_nan());
    }
    std::cout << std::endl;
}

// =====================================================================
// Benchmark: pcSENSE reconstruction with linear ramp phase maps
// =====================================================================
TEST_CASE("pcSENSE recon benchmark", "[pcSENSE_bench][benchmark]") {
    typedef float T1;
    typedef std::complex<T1> CxT1;

    struct BenchConfig {
        uword Nx, Ny, Nc, Ns, nSpokes, nReadout, niter;
        const char* label;
    };

    BenchConfig configs[] = {
        {32, 32, 4, 2, 32, 32, 10, "32x32, Nc=4, Ns=2"},
        {32, 32, 4, 4, 32, 32, 10, "32x32, Nc=4, Ns=4"},
        {64, 64, 4, 2, 32, 32, 10, "64x64, Nc=4, Ns=2"},
        {128, 128, 4, 2, 32, 32, 5, "128x128, Nc=4, Ns=2"},
    };

    std::cout << std::endl;
    std::cout << "pcSENSE reconstruction benchmark (linear ramp phase):" << std::endl;
    std::cout << "  " << std::left << std::setw(28) << "Config"
              << std::setw(10) << "Ni"
              << std::setw(10) << "Nd*Ns*Nc"
              << std::setw(14) << "Setup (ms)"
              << std::setw(14) << "Recon (ms)"
              << std::setw(10) << "NRMSE"
              << std::endl;

    for (const auto& cfg : configs) {
        uword Ni = cfg.Nx * cfg.Ny;
        uword Nd = cfg.nSpokes * cfg.nReadout;

        // Phantom
        Col<CxT1> phantom = sheppLogan2D<T1>(cfg.Nx, cfg.Ny);

        // Trajectory: same for all shots
        Col<T1> kx1, ky1, kz1;
        radialTrajectory2D<T1>(cfg.nSpokes, cfg.nReadout, cfg.Nx, kx1, ky1, kz1);
        Col<T1> kx = repmat(kx1, cfg.Ns, 1);
        Col<T1> ky = repmat(ky1, cfg.Ns, 1);
        Col<T1> kz = repmat(kz1, cfg.Ns, 1);

        // Timing vector
        Col<T1> tvec1 = linearTimingVector<T1>(Nd, (T1)0.005);
        Col<T1> tvec = repmat(tvec1, cfg.Ns, 1);

        // Coil maps
        Col<CxT1> SENSEmap = syntheticCoils2D<T1>(cfg.Nx, cfg.Ny, cfg.Nc);

        // Field map
        Col<T1> fmapRad = syntheticFieldMap2D<T1>(cfg.Nx, cfg.Ny, (T1)10.0) * (T1)(2.0 * M_PI);

        // Linear ramp phase
        Col<T1> shotPhase = linearRampPhaseMap2D<T1>(cfg.Nx, cfg.Ny, cfg.Ns, (T1)(M_PI / 4.0));

        auto t0 = std::chrono::high_resolution_clock::now();

        // pcSENSE operator
        pcSENSE<T1> P(kx, ky, kz, cfg.Nx, cfg.Ny, 1, cfg.Nc, tvec,
                      SENSEmap, fmapRad, shotPhase);

        auto t_setup = std::chrono::high_resolution_clock::now();

        // Generate data
        Col<CxT1> y = P * phantom;

        // Weights and regularizer
        Col<T1> W = ones<Col<T1>>(Nd * cfg.Ns * cfg.Nc);
        QuadPenalty<T1> R(cfg.Nx, cfg.Ny, 1, (T1)1e-3, 2);
        Col<CxT1> x0 = zeros<Col<CxT1>>(Ni);

        auto t_recon_start = std::chrono::high_resolution_clock::now();

        // Reconstruct
        Col<CxT1> xhat = solve_pwls_pcg<T1>(x0, P, W, y, R, cfg.niter);

        auto t_recon_end = std::chrono::high_resolution_clock::now();

        T1 err = nrmse_gdft<T1>(xhat, phantom);

        double setup_ms = std::chrono::duration<double, std::milli>(t_setup - t0).count();
        double recon_ms = std::chrono::duration<double, std::milli>(t_recon_end - t_recon_start).count();

        std::cout << "  " << std::left << std::setw(28) << cfg.label
                  << std::setw(10) << Ni
                  << std::setw(10) << (Nd * cfg.Ns * cfg.Nc)
                  << std::setw(14) << std::fixed << std::setprecision(1) << setup_ms
                  << std::setw(14) << recon_ms
                  << std::setw(10) << std::setprecision(4) << err
                  << std::endl;

        REQUIRE(!xhat.has_nan());
        REQUIRE(err < (T1)0.80);
    }
    std::cout << std::endl;
}

// =====================================================================
// Benchmark: pcSENSE 256x256 reconstruction (tagged [bench256] for
// comparison with SENSE/TimeSeg benchmarks at the same image size)
// =====================================================================
TEST_CASE("pcSENSE recon 256x256 benchmark", "[bench256][pcSENSE_bench256]") {
    typedef float T1;
    typedef std::complex<T1> CxT1;

    const uword Nx = 256;
    const uword Ny = 256;
    const uword Ni = Nx * Ny;
    const uword Nc = 4;
    const uword Ns = 2;
    const uword nSpokes = 32;
    const uword nReadout = 32;
    const uword Nd = nSpokes * nReadout;
    const uword niter = 5;

    auto t0 = std::chrono::high_resolution_clock::now();

    // Phantom
    Col<CxT1> phantom = sheppLogan2D<T1>(Nx, Ny);

    // Trajectory: same for all shots
    Col<T1> kx1, ky1, kz1;
    radialTrajectory2D<T1>(nSpokes, nReadout, Nx, kx1, ky1, kz1);
    Col<T1> kx = repmat(kx1, Ns, 1);
    Col<T1> ky = repmat(ky1, Ns, 1);
    Col<T1> kz = repmat(kz1, Ns, 1);

    // Timing vector
    Col<T1> tvec1 = linearTimingVector<T1>(Nd, (T1)0.005);
    Col<T1> tvec = repmat(tvec1, Ns, 1);

    // Coil maps
    Col<CxT1> SENSEmap = syntheticCoils2D<T1>(Nx, Ny, Nc);

    // Field map
    Col<T1> fmapRad = syntheticFieldMap2D<T1>(Nx, Ny, (T1)10.0) * (T1)(2.0 * M_PI);

    // Linear ramp shot phase: pi/4 radians max
    Col<T1> shotPhase = linearRampPhaseMap2D<T1>(Nx, Ny, Ns, (T1)(M_PI / 4.0));

    // pcSENSE operator
    pcSENSE<T1> P(kx, ky, kz, Nx, Ny, 1, Nc, tvec, SENSEmap, fmapRad, shotPhase);

    auto t_setup = std::chrono::high_resolution_clock::now();

    // Generate data
    Col<CxT1> y = P * phantom;

    // Weights and regularizer
    Col<T1> W = ones<Col<T1>>(Nd * Ns * Nc);
    QuadPenalty<T1> R(Nx, Ny, 1, (T1)1e-3, 2);
    Col<CxT1> x0 = zeros<Col<CxT1>>(Ni);

    auto t_recon_start = std::chrono::high_resolution_clock::now();

    // Reconstruct
    Col<CxT1> xhat = solve_pwls_pcg<T1>(x0, P, W, y, R, niter);

    auto t_recon_end = std::chrono::high_resolution_clock::now();

    REQUIRE(xhat.n_elem == Ni);
    REQUIRE(!xhat.has_nan());

    T1 err = nrmse_gdft<T1>(xhat, phantom);

    double setup_ms = std::chrono::duration<double, std::milli>(t_setup - t0).count();
    double recon_ms = std::chrono::duration<double, std::milli>(t_recon_end - t_recon_start).count();
    double total_ms = std::chrono::duration<double, std::milli>(t_recon_end - t0).count();

    std::cout << "pcSENSE 256x256 (Nc=" << Nc << ", Ns=" << Ns
              << ", " << nSpokes << " spokes, " << niter << " iters):" << std::endl;
    std::cout << "  Setup:  " << setup_ms << " ms" << std::endl;
    std::cout << "  Recon:  " << recon_ms << " ms" << std::endl;
    std::cout << "  Total:  " << total_ms << " ms" << std::endl;
    std::cout << "  NRMSE:  " << err << std::endl;

    // Loose threshold: 1024 samples for 65536 pixels (64x undersampling)
    REQUIRE(err < (T1)0.95);
}
