/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file SyntheticReconTests.cpp
/// @brief Integration tests for SENSE, pcSENSE, and TimeSegmentation using
///        synthetic phantoms, coil maps, and trajectories.

#include "catch.hpp"

#include "../PGIncludes.h"
#include "../Gnufft.h"
#include "../SENSE.h"
#include "../TimeSegmentation.h"
#include "../QuadPenalty.h"
#include "../solve_pwls_pcg.hpp"

#include "SyntheticPhantom.hpp"
#include "SyntheticCoils.hpp"
#include "SyntheticTrajectory.hpp"
#include "SyntheticFieldMap.hpp"

using namespace arma;

// ---------------------------------------------------------------------------
// Helper: compute image-space coordinates matching Gnufft's expected layout
// (same as pcSENSE constructor, based on Cube<T1>(Ny, Nx, Nz) indexing)
// ---------------------------------------------------------------------------
template<typename T1>
void gnufftImageCoords2D(uword Nx, uword Ny,
                         Col<T1>& ix, Col<T1>& iy, Col<T1>& iz)
{
    uword Ni = Nx * Ny;
    ix.set_size(Ni);
    iy.set_size(Ni);
    iz.zeros(Ni);

    // pcSENSE uses Cube(Nx, Ny, Nz) with loops: outer y (ii), inner x (jj)
    // ix(ii, jj, kk) = (jj - Nx/2) / Nx  (x varies with inner index)
    // iy(ii, jj, kk) = (ii - Ny/2) / Ny  (y varies with outer index)
    // After vectorise, the Cube stores column-major: (row, col, slice)
    for (uword jj = 0; jj < Nx; jj++) {       // x (column of cube)
        for (uword ii = 0; ii < Ny; ii++) {    // y (row of cube)
            uword idx = ii + jj * Ny;
            ix(idx) = ((T1)jj - (T1)Nx / (T1)2.0) / (T1)Nx;
            iy(idx) = ((T1)ii - (T1)Ny / (T1)2.0) / (T1)Ny;
        }
    }
}

// ---------------------------------------------------------------------------
// Helper: Normalized Root Mean Square Error between two complex vectors
// ---------------------------------------------------------------------------
template<typename T1>
T1 nrmse(const Col<std::complex<T1>>& x, const Col<std::complex<T1>>& ref) {
    T1 errNorm = norm(x - ref, 2);
    T1 refNorm = norm(ref, 2);
    if (refNorm < (T1)1e-20) return errNorm;
    return errNorm / refNorm;
}

// ---------------------------------------------------------------------------
// Test 1: SENSE forward-adjoint consistency (adjointness / dot-product test)
// ---------------------------------------------------------------------------
TEST_CASE("SENSE adjointness with radial trajectory", "[SENSE_synthetic][adjointness]") {
    typedef float T1;
    typedef std::complex<T1> CxT1;

    const uword Nx = 32;
    const uword Ny = 32;
    const uword Ni = Nx * Ny;
    const uword Nc = 4;
    const uword nSpokes = 64;
    const uword nReadout = 64;
    const uword Nd = nSpokes * nReadout;

    // Generate trajectory
    Col<T1> kx, ky, kz;
    radialTrajectory2D<T1>(nSpokes, nReadout, Nx, kx, ky, kz);

    // Image space coordinates
    Col<T1> ix, iy, iz;
    gnufftImageCoords2D<T1>(Nx, Ny, ix, iy, iz);

    // Create Gnufft operator
    Gnufft<T1> G(Nd, (T1)2.0, Nx, Ny, 1, kx, ky, kz, ix, iy, iz);

    // Coil sensitivity maps
    Col<CxT1> SENSEmap = syntheticCoils2D<T1>(Nx, Ny, Nc);

    // Create SENSE operator
    SENSE<T1, Gnufft<T1>> S(G, SENSEmap, Nd, Ni, Nc);

    // Generate random vectors
    arma_rng::set_seed(42);
    Col<CxT1> x = randn<Col<CxT1>>(Ni);
    Col<CxT1> y = randn<Col<CxT1>>(Nd * Nc);

    // Forward: Ax
    Col<CxT1> Ax = S * x;
    REQUIRE(Ax.n_elem == Nd * Nc);

    // Adjoint: A'y
    Col<CxT1> Aty = S / y;
    REQUIRE(Aty.n_elem == Ni);

    // Dot product test: <Ax, y> should approximately equal <x, A'y>
    CxT1 lhs = cdot(Ax, y);
    CxT1 rhs = cdot(x, Aty);

    T1 relErr = std::abs(lhs - rhs) / (std::abs(lhs) + (T1)1e-10);
    std::cout << "SENSE adjointness: |<Ax,y> - <x,A'y>| / |<Ax,y>| = " << relErr << std::endl;
    REQUIRE(relErr < (T1)1e-3);
}

// ---------------------------------------------------------------------------
// Test 2: SENSE reconstruction from synthetic data
// ---------------------------------------------------------------------------
TEST_CASE("SENSE reconstruction with Shepp-Logan phantom", "[SENSE_synthetic][recon]") {
    typedef float T1;
    typedef std::complex<T1> CxT1;

    const uword Nx = 64;
    const uword Ny = 64;
    const uword Ni = Nx * Ny;
    const uword Nc = 4;
    const uword nSpokes = 128;
    const uword nReadout = 128;
    const uword Nd = nSpokes * nReadout;
    const uword niter = 10;

    // Generate phantom
    Col<CxT1> phantom = sheppLogan2D<T1>(Nx, Ny);
    REQUIRE(phantom.n_elem == Ni);
    REQUIRE(norm(phantom, 2) > (T1)0.0);

    // Generate trajectory
    Col<T1> kx, ky, kz;
    radialTrajectory2D<T1>(nSpokes, nReadout, Nx, kx, ky, kz);

    // Image space coordinates
    Col<T1> ix, iy, iz;
    gnufftImageCoords2D<T1>(Nx, Ny, ix, iy, iz);

    // Create Gnufft
    Gnufft<T1> G(Nd, (T1)2.0, Nx, Ny, 1, kx, ky, kz, ix, iy, iz);

    // Coil maps
    Col<CxT1> SENSEmap = syntheticCoils2D<T1>(Nx, Ny, Nc);

    // SENSE operator
    SENSE<T1, Gnufft<T1>> S(G, SENSEmap, Nd, Ni, Nc);

    // Generate noiseless data: y = S * phantom
    Col<CxT1> y = S * phantom;
    REQUIRE(y.n_elem == Nd * Nc);
    REQUIRE(!y.has_nan());

    // Uniform weights
    Col<T1> W = ones<Col<T1>>(Nd * Nc);

    // Regularizer
    QuadPenalty<T1> R(Nx, Ny, 1, (T1)1e-3, 2);

    // Initial estimate: zero
    Col<CxT1> x0 = zeros<Col<CxT1>>(Ni);

    // Reconstruct
    Col<CxT1> xhat = solve_pwls_pcg<T1>(x0, S, W, y, R, niter);
    REQUIRE(xhat.n_elem == Ni);
    REQUIRE(!xhat.has_nan());

    // Check NRMSE
    T1 err = nrmse<T1>(xhat, phantom);
    std::cout << "SENSE recon NRMSE = " << err << std::endl;
    // Loose tolerance: radial undersampling + limited iterations + regularization
    REQUIRE(err < (T1)0.50);
}

// ---------------------------------------------------------------------------
// Test 3: SENSE forward-adjoint with TimeSegmentation
// ---------------------------------------------------------------------------
TEST_CASE("SENSE+TimeSeg adjointness", "[TimeSeg_synthetic][adjointness]") {
    typedef float T1;
    typedef std::complex<T1> CxT1;

    const uword Nx = 32;
    const uword Ny = 32;
    const uword Ni = Nx * Ny;
    const uword nSpokes = 64;
    const uword nReadout = 64;
    const uword Nd = nSpokes * nReadout;
    const uword Nc = 4;
    const uword L = 4;   // time segments

    // Trajectory
    Col<T1> kx, ky, kz;
    radialTrajectory2D<T1>(nSpokes, nReadout, Nx, kx, ky, kz);

    // Image coordinates
    Col<T1> ix, iy, iz;
    gnufftImageCoords2D<T1>(Nx, Ny, ix, iy, iz);

    // Gnufft
    Gnufft<T1> G(Nd, (T1)2.0, Nx, Ny, 1, kx, ky, kz, ix, iy, iz);

    // Field map: ±50 Hz quadratic, convert to radians/sec
    Col<T1> fmapHz = syntheticFieldMap2D<T1>(Nx, Ny, (T1)50.0);
    Col<T1> fmapRad = fmapHz * (T1)(2.0 * M_PI);

    // Timing vector
    Col<T1> tvec = linearTimingVector<T1>(Nd, (T1)0.01);

    // TimeSegmentation operator
    TimeSegmentation<T1, Gnufft<T1>> TS(G, fmapRad, tvec, Nd, Ni, L, 1, 1);

    // Coil maps
    Col<CxT1> SENSEmap = syntheticCoils2D<T1>(Nx, Ny, Nc);

    // SENSE with TimeSegmentation
    SENSE<T1, TimeSegmentation<T1, Gnufft<T1>>> S(TS, SENSEmap, Nd, Ni, Nc);

    // Random vectors
    arma_rng::set_seed(123);
    Col<CxT1> x = randn<Col<CxT1>>(Ni);
    Col<CxT1> y = randn<Col<CxT1>>(Nd * Nc);

    // Adjointness check
    Col<CxT1> Ax = S * x;
    Col<CxT1> Aty = S / y;

    CxT1 lhs = cdot(Ax, y);
    CxT1 rhs = cdot(x, Aty);

    T1 relErr = std::abs(lhs - rhs) / (std::abs(lhs) + (T1)1e-10);
    std::cout << "SENSE+TimeSeg adjointness: relErr = " << relErr << std::endl;
    REQUIRE(relErr < (T1)1e-2);
}

// ---------------------------------------------------------------------------
// Test 4: Gnufft round-trip (forward then adjoint recovers something useful)
// ---------------------------------------------------------------------------
TEST_CASE("Gnufft forward-adjoint round-trip", "[SENSE_synthetic][Gnufft]") {
    typedef float T1;
    typedef std::complex<T1> CxT1;

    const uword Nx = 32;
    const uword Ny = 32;
    const uword Ni = Nx * Ny;
    const uword nSpokes = 64;
    const uword nReadout = 64;
    const uword Nd = nSpokes * nReadout;

    // Trajectory
    Col<T1> kx, ky, kz;
    radialTrajectory2D<T1>(nSpokes, nReadout, Nx, kx, ky, kz);

    // Image coordinates
    Col<T1> ix, iy, iz;
    gnufftImageCoords2D<T1>(Nx, Ny, ix, iy, iz);

    // Gnufft
    Gnufft<T1> G(Nd, (T1)2.0, Nx, Ny, 1, kx, ky, kz, ix, iy, iz);

    // Phantom
    Col<CxT1> phantom = sheppLogan2D<T1>(Nx, Ny);

    // Forward
    Col<CxT1> kdata = G * phantom;
    REQUIRE(kdata.n_elem == Nd);
    REQUIRE(!kdata.has_nan());

    // Adjoint
    Col<CxT1> img = G / kdata;
    REQUIRE(img.n_elem == Ni);
    REQUIRE(!img.has_nan());

    // The adjoint of forward should recover something correlated with the original
    // Compute correlation: real(dot(img, phantom)) / (norm(img) * norm(phantom))
    T1 corr = std::abs(cdot(img, phantom)) / (norm(img, 2) * norm(phantom, 2));
    std::cout << "Gnufft round-trip correlation = " << corr << std::endl;
    REQUIRE(corr > (T1)0.1);
}

// ---------------------------------------------------------------------------
// Test 5: Shepp-Logan phantom sanity checks
// ---------------------------------------------------------------------------
TEST_CASE("Shepp-Logan phantom properties", "[SENSE_synthetic][phantom]") {
    typedef float T1;
    typedef std::complex<T1> CxT1;

    const uword Nx = 128;
    const uword Ny = 128;

    Col<CxT1> phantom = sheppLogan2D<T1>(Nx, Ny);
    REQUIRE(phantom.n_elem == Nx * Ny);

    // Should have nonzero values (the phantom isn't all black)
    REQUIRE(norm(phantom, 2) > (T1)0.0);

    // Real parts should be in approximately [0, 1] (modified Shepp-Logan)
    Col<T1> realPart = real(phantom);
    REQUIRE(realPart.min() >= (T1)-0.01);
    REQUIRE(realPart.max() <= (T1)1.01);

    // Imaginary parts should be zero
    Col<T1> imagPart = imag(phantom);
    REQUIRE(norm(imagPart, 2) < (T1)1e-10);

    // Center pixel should be nonzero (inside the skull)
    uword centerIdx = Nx / 2 + (Ny / 2) * Nx;
    REQUIRE(std::abs(phantom(centerIdx)) > (T1)0.01);
}

// ---------------------------------------------------------------------------
// Test 6: Coil map properties
// ---------------------------------------------------------------------------
TEST_CASE("Synthetic coil map properties", "[SENSE_synthetic][coils]") {
    typedef float T1;
    typedef std::complex<T1> CxT1;

    const uword Nx = 64;
    const uword Ny = 64;
    const uword Nc = 8;
    const uword Ni = Nx * Ny;

    Col<CxT1> maps = syntheticCoils2D<T1>(Nx, Ny, Nc);
    REQUIRE(maps.n_elem == Ni * Nc);

    // Reshape to matrix for easier checking
    Mat<CxT1> SMap = reshape(maps, Ni, Nc);

    // Sum-of-squares should be approximately 1 everywhere
    Col<T1> sos(Ni, fill::zeros);
    for (uword cc = 0; cc < Nc; cc++) {
        sos += real(SMap.col(cc) % conj(SMap.col(cc)));
    }
    // After normalization, SoS should be close to 1 at center voxels
    // (edge voxels far from all coils may have near-zero sensitivity)
    uword centerIdx = Nx / 2 + (Ny / 2) * Nx;
    REQUIRE(std::abs(sos(centerIdx) - (T1)1.0) < (T1)0.01);
    REQUIRE(std::abs(sos.max() - (T1)1.0) < (T1)0.01);

    // Each coil should have different phase patterns
    // (check that not all maps are identical)
    T1 diffNorm = norm(SMap.col(0) - SMap.col(1), 2);
    REQUIRE(diffNorm > (T1)0.1);
}

// ---------------------------------------------------------------------------
// Test 7: Radial trajectory properties
// ---------------------------------------------------------------------------
TEST_CASE("Radial trajectory properties", "[SENSE_synthetic][trajectory]") {
    typedef float T1;

    const uword nSpokes = 128;
    const uword nReadout = 256;
    const uword Nx = 256;

    Col<T1> kx, ky, kz;
    radialTrajectory2D<T1>(nSpokes, nReadout, Nx, kx, ky, kz);

    REQUIRE(kx.n_elem == nSpokes * nReadout);
    REQUIRE(ky.n_elem == nSpokes * nReadout);

    // k-space coordinates should be within [-Nx/2, Nx/2]
    T1 kmax = (T1)Nx / (T1)2.0;
    REQUIRE(kx.max() <= kmax * (T1)1.01);
    REQUIRE(kx.min() >= -kmax * (T1)1.01);
    REQUIRE(ky.max() <= kmax * (T1)1.01);
    REQUIRE(ky.min() >= -kmax * (T1)1.01);

    // kz should be all zeros for 2D
    REQUIRE(norm(kz, 2) < (T1)1e-10);
}

// ---------------------------------------------------------------------------
// Test 8: Spiral trajectory properties
// ---------------------------------------------------------------------------
TEST_CASE("Spiral trajectory properties", "[SENSE_synthetic][trajectory]") {
    typedef float T1;

    const uword nInterleaves = 16;
    const uword nSamplesPerArm = 512;
    const uword Nx = 128;

    Col<T1> kx, ky, kz;
    spiralTrajectory2D<T1>(nInterleaves, nSamplesPerArm, Nx, kx, ky, kz);

    REQUIRE(kx.n_elem == nInterleaves * nSamplesPerArm);

    // Spiral should start at center (k=0) for first interleave
    REQUIRE(std::abs(kx(0)) < (T1)1e-5);
    REQUIRE(std::abs(ky(0)) < (T1)1e-5);

    // Should reach near kmax at the end of each arm
    T1 kmax = (T1)Nx / (T1)2.0;
    T1 lastR = std::sqrt(kx(nSamplesPerArm - 1) * kx(nSamplesPerArm - 1) +
                         ky(nSamplesPerArm - 1) * ky(nSamplesPerArm - 1));
    REQUIRE(lastR > kmax * (T1)0.9);

    // kz should be all zeros
    REQUIRE(norm(kz, 2) < (T1)1e-10);
}

// ---------------------------------------------------------------------------
// Test 9: Field map properties
// ---------------------------------------------------------------------------
TEST_CASE("Synthetic field map properties", "[SENSE_synthetic][fieldmap]") {
    typedef float T1;

    const uword Nx = 64;
    const uword Ny = 64;
    T1 maxHz = (T1)100.0;

    Col<T1> fmap = syntheticFieldMap2D<T1>(Nx, Ny, maxHz);
    REQUIRE(fmap.n_elem == Nx * Ny);

    // Field map should be non-negative (quadratic)
    REQUIRE(fmap.min() >= (T1)-0.01);

    // Maximum should be close to maxHz (at corners)
    REQUIRE(fmap.max() <= maxHz * (T1)1.6);  // x^2 + 0.5*y^2 can reach ~1.5 at corners

    // Center should be near zero
    uword centerIdx = Nx / 2 + (Ny / 2) * Nx;
    REQUIRE(std::abs(fmap(centerIdx)) < maxHz * (T1)0.05);
}

// ---------------------------------------------------------------------------
// Test 10: SENSE with spiral trajectory adjointness
// ---------------------------------------------------------------------------
TEST_CASE("SENSE adjointness with spiral trajectory", "[SENSE_synthetic][adjointness]") {
    typedef float T1;
    typedef std::complex<T1> CxT1;

    const uword Nx = 32;
    const uword Ny = 32;
    const uword Ni = Nx * Ny;
    const uword Nc = 4;
    const uword nInterleaves = 8;
    const uword nSamplesPerArm = 256;
    const uword Nd = nInterleaves * nSamplesPerArm;

    // Spiral trajectory
    Col<T1> kx, ky, kz;
    spiralTrajectory2D<T1>(nInterleaves, nSamplesPerArm, Nx, kx, ky, kz);

    // Image coordinates
    Col<T1> ix, iy, iz;
    gnufftImageCoords2D<T1>(Nx, Ny, ix, iy, iz);

    // Gnufft
    Gnufft<T1> G(Nd, (T1)2.0, Nx, Ny, 1, kx, ky, kz, ix, iy, iz);

    // Coil maps
    Col<CxT1> SENSEmap = syntheticCoils2D<T1>(Nx, Ny, Nc);

    // SENSE
    SENSE<T1, Gnufft<T1>> S(G, SENSEmap, Nd, Ni, Nc);

    // Random vectors
    arma_rng::set_seed(99);
    Col<CxT1> x = randn<Col<CxT1>>(Ni);
    Col<CxT1> y = randn<Col<CxT1>>(Nd * Nc);

    Col<CxT1> Ax = S * x;
    Col<CxT1> Aty = S / y;

    CxT1 lhs = cdot(Ax, y);
    CxT1 rhs = cdot(x, Aty);

    T1 relErr = std::abs(lhs - rhs) / (std::abs(lhs) + (T1)1e-10);
    std::cout << "SENSE+spiral adjointness: relErr = " << relErr << std::endl;
    REQUIRE(relErr < (T1)1e-3);
}

// ---------------------------------------------------------------------------
// Benchmark: SENSE reconstruction at 256x256
// ---------------------------------------------------------------------------
TEST_CASE("SENSE recon 256x256 benchmark", "[bench256][SENSE_bench]") {
    typedef float T1;
    typedef std::complex<T1> CxT1;

    const uword Nx = 256;
    const uword Ny = 256;
    const uword Ni = Nx * Ny;
    const uword Nc = 8;
    const uword nSpokes = 256;
    const uword nReadout = 256;
    const uword Nd = nSpokes * nReadout;
    const uword niter = 10;

    auto t0 = std::chrono::high_resolution_clock::now();

    // Generate phantom
    Col<CxT1> phantom = sheppLogan2D<T1>(Nx, Ny);

    // Generate trajectory
    Col<T1> kx, ky, kz;
    radialTrajectory2D<T1>(nSpokes, nReadout, Nx, kx, ky, kz);

    // Image space coordinates
    Col<T1> ix, iy, iz;
    gnufftImageCoords2D<T1>(Nx, Ny, ix, iy, iz);

    // Create Gnufft
    Gnufft<T1> G(Nd, (T1)2.0, Nx, Ny, 1, kx, ky, kz, ix, iy, iz);

    // Coil maps
    Col<CxT1> SENSEmap = syntheticCoils2D<T1>(Nx, Ny, Nc);

    // SENSE operator
    SENSE<T1, Gnufft<T1>> S(G, SENSEmap, Nd, Ni, Nc);

    auto t_setup = std::chrono::high_resolution_clock::now();

    // Generate noiseless data: y = S * phantom
    Col<CxT1> y = S * phantom;

    // Uniform weights
    Col<T1> W = ones<Col<T1>>(Nd * Nc);

    // Regularizer
    QuadPenalty<T1> R(Nx, Ny, 1, (T1)1e-3, 2);

    // Initial estimate: zero
    Col<CxT1> x0 = zeros<Col<CxT1>>(Ni);

    auto t_recon_start = std::chrono::high_resolution_clock::now();

    // Reconstruct
    Col<CxT1> xhat = solve_pwls_pcg<T1>(x0, S, W, y, R, niter);

    auto t_recon_end = std::chrono::high_resolution_clock::now();

    REQUIRE(xhat.n_elem == Ni);
    REQUIRE(!xhat.has_nan());

    T1 err = nrmse<T1>(xhat, phantom);

    double setup_ms = std::chrono::duration<double, std::milli>(t_setup - t0).count();
    double recon_ms = std::chrono::duration<double, std::milli>(t_recon_end - t_recon_start).count();
    double total_ms = std::chrono::duration<double, std::milli>(t_recon_end - t0).count();

    std::cout << "SENSE 256x256 (8 coils, 256 spokes, 10 iters):" << std::endl;
    std::cout << "  Setup:  " << setup_ms << " ms" << std::endl;
    std::cout << "  Recon:  " << recon_ms << " ms" << std::endl;
    std::cout << "  Total:  " << total_ms << " ms" << std::endl;
    std::cout << "  NRMSE:  " << err << std::endl;

    REQUIRE(err < (T1)0.50);
}

// ---------------------------------------------------------------------------
// Benchmark: SENSE + TimeSegmentation reconstruction at 256x256
// ---------------------------------------------------------------------------
TEST_CASE("SENSE+TimeSeg recon 256x256 benchmark", "[bench256][TimeSeg_bench]") {
    typedef float T1;
    typedef std::complex<T1> CxT1;

    const uword Nx = 256;
    const uword Ny = 256;
    const uword Ni = Nx * Ny;
    const uword Nc = 8;
    const uword nSpokes = 256;
    const uword nReadout = 256;
    const uword Nd = nSpokes * nReadout;
    const uword L = 4;
    const uword niter = 10;

    auto t0 = std::chrono::high_resolution_clock::now();

    // Generate phantom
    Col<CxT1> phantom = sheppLogan2D<T1>(Nx, Ny);

    // Trajectory
    Col<T1> kx, ky, kz;
    radialTrajectory2D<T1>(nSpokes, nReadout, Nx, kx, ky, kz);

    // Image coordinates
    Col<T1> ix, iy, iz;
    gnufftImageCoords2D<T1>(Nx, Ny, ix, iy, iz);

    // Gnufft
    Gnufft<T1> G(Nd, (T1)2.0, Nx, Ny, 1, kx, ky, kz, ix, iy, iz);

    // Field map: ±50 Hz quadratic, convert to radians/sec
    Col<T1> fmapHz = syntheticFieldMap2D<T1>(Nx, Ny, (T1)50.0);
    Col<T1> fmapRad = fmapHz * (T1)(2.0 * M_PI);

    // Timing vector
    Col<T1> tvec = linearTimingVector<T1>(Nd, (T1)0.01);

    // TimeSegmentation operator
    TimeSegmentation<T1, Gnufft<T1>> TS(G, fmapRad, tvec, Nd, Ni, L, 1, 1);

    // Coil maps
    Col<CxT1> SENSEmap = syntheticCoils2D<T1>(Nx, Ny, Nc);

    // SENSE with TimeSegmentation
    SENSE<T1, TimeSegmentation<T1, Gnufft<T1>>> S(TS, SENSEmap, Nd, Ni, Nc);

    auto t_setup = std::chrono::high_resolution_clock::now();

    // Generate data
    Col<CxT1> y = S * phantom;

    // Weights and regularizer
    Col<T1> W = ones<Col<T1>>(Nd * Nc);
    QuadPenalty<T1> R(Nx, Ny, 1, (T1)1e-3, 2);
    Col<CxT1> x0 = zeros<Col<CxT1>>(Ni);

    auto t_recon_start = std::chrono::high_resolution_clock::now();

    // Reconstruct
    Col<CxT1> xhat = solve_pwls_pcg<T1>(x0, S, W, y, R, niter);

    auto t_recon_end = std::chrono::high_resolution_clock::now();

    REQUIRE(xhat.n_elem == Ni);
    REQUIRE(!xhat.has_nan());

    T1 err = nrmse<T1>(xhat, phantom);

    double setup_ms = std::chrono::duration<double, std::milli>(t_setup - t0).count();
    double recon_ms = std::chrono::duration<double, std::milli>(t_recon_end - t_recon_start).count();
    double total_ms = std::chrono::duration<double, std::milli>(t_recon_end - t0).count();

    std::cout << "SENSE+TimeSeg 256x256 (8 coils, L=4, 256 spokes, 10 iters):" << std::endl;
    std::cout << "  Setup:  " << setup_ms << " ms" << std::endl;
    std::cout << "  Recon:  " << recon_ms << " ms" << std::endl;
    std::cout << "  Total:  " << total_ms << " ms" << std::endl;
    std::cout << "  NRMSE:  " << err << std::endl;

    REQUIRE(err < (T1)0.50);
}
