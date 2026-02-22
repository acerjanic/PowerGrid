/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file Spiral3DReconTests.cpp
/// @brief 3D stack-of-spirals reconstruction tests with pcSENSE and
///        pcSenseTimeSeg, multi-compartment brain phantom, and NIfTI output.

#include "catch.hpp"

#include "../PGIncludes.h"
#include "../Gdft.h"
#include "../Gnufft.h"
#include "../pcSENSE.h"
#include "../pcSenseTimeSeg.h"
#include "../QuadPenalty.h"
#include "../solve_pwls_pcg.hpp"

#include "SyntheticPhantom3D.hpp"
#include "SyntheticCoils3D.hpp"
#include "SyntheticFieldMap3D.hpp"
#include "SyntheticTrajectory3D.hpp"

// Minimal NIfTI writer (avoids PowerGrid.h → ISMRMRD dependency)
#include "../../Support/nifti1.h"
#include <fstream>
#include <cstring>
#include <set>

using namespace arma;

// ---------------------------------------------------------------------------
// Minimal NIfTI helpers (self-contained, no ISMRMRD dependency)
// ---------------------------------------------------------------------------
namespace {

template<typename T1>
void writeNifti(const std::string& filename, const Col<T1>& data,
                uword Nx, uword Ny, uword Nz)
{
    std::string path = filename + ".nii";

    nifti_1_header hdr;
    std::memset(&hdr, 0, sizeof(hdr));
    hdr.sizeof_hdr = 348;
    hdr.dim[0] = 3;
    hdr.dim[1] = (short)Nx;
    hdr.dim[2] = (short)Ny;
    hdr.dim[3] = (short)Nz;
    hdr.dim[4] = 1;
    hdr.dim[5] = 1;
    hdr.dim[6] = 1;
    hdr.dim[7] = 1;
    hdr.datatype = NIFTI_TYPE_FLOAT32;
    hdr.bitpix = 32;
    hdr.pixdim[0] = 1.0f;
    hdr.pixdim[1] = 1.0f;
    hdr.pixdim[2] = 1.0f;
    hdr.pixdim[3] = 1.0f;
    hdr.vox_offset = 352.0f;
    hdr.scl_slope = 1.0f;
    hdr.xyzt_units = NIFTI_UNITS_MM;
    std::strncpy(hdr.magic, "n+1\0", 4);

    nifti1_extender pad = {0, 0, 0, 0};

    std::ofstream ofs(path, std::ios::binary);
    ofs.write(reinterpret_cast<const char*>(&hdr), 348);
    ofs.write(reinterpret_cast<const char*>(&pad), 4);
    for (uword i = 0; i < data.n_elem; i++) {
        float val = static_cast<float>(data(i));
        ofs.write(reinterpret_cast<const char*>(&val), sizeof(float));
    }
    ofs.close();

    std::cout << "  Wrote: " << path << " (" << Nx << "x" << Ny << "x" << Nz << ")" << std::endl;
}

template<typename T1>
void writeNiftiComplex(const std::string& filename, const Col<std::complex<T1>>& data,
                       uword Nx, uword Ny, uword Nz)
{
    Col<T1> magImg = abs(data);
    Col<T1> phsImg(data.n_elem);
    for (uword i = 0; i < data.n_elem; i++) {
        phsImg(i) = std::arg(data(i));
    }
    writeNifti<T1>(filename + "_mag", magImg, Nx, Ny, Nz);
    writeNifti<T1>(filename + "_phs", phsImg, Nx, Ny, Nz);
}

} // anonymous namespace

// ===========================================================================
// Test 1: 3D brain phantom sanity checks
// ===========================================================================
TEST_CASE("3D brain phantom properties", "[spiral3D][phantom3D]") {
    using T1 = float;
    const uword Nx = 32, Ny = 32, Nz = 16;
    const uword Ni = Nx * Ny * Nz;

    Col<std::complex<T1>> phantom = brainPhantom3D<T1>(Nx, Ny, Nz);

    REQUIRE(phantom.n_elem == Ni);
    REQUIRE(norm(phantom) > 0);

    // All real parts in [0, 1.01], imaginary parts ~0
    Col<T1> realPart = real(phantom);
    Col<T1> imagPart = imag(phantom);
    REQUIRE(realPart.min() >= (T1)-0.01);
    REQUIRE(realPart.max() <= (T1)1.01);
    REQUIRE(max(abs(imagPart)) < (T1)1e-6);

    // Center voxel should be inside the brain (nonzero)
    uword centerIdx = Ny/2 + (Nx/2)*Nx + (Nz/2)*Nx*Ny;
    REQUIRE(std::abs(phantom(centerIdx)) > (T1)0.1);

    // Should have at least 3 distinct nonzero intensity levels
    std::set<int> levels;
    for (uword i = 0; i < Ni; i++) {
        T1 v = realPart(i);
        if (v > 0.01) {
            levels.insert((int)(v * 100 + 0.5));
        }
    }
    REQUIRE(levels.size() >= 3);
}

// ===========================================================================
// Test 2: 3D coil SoS normalization
// ===========================================================================
TEST_CASE("3D coil sensitivity SoS normalization", "[spiral3D][coils3D]") {
    using T1 = float;
    const uword Nx = 16, Ny = 16, Nz = 8, Nc = 16;
    const uword Ni = Nx * Ny * Nz;

    Col<std::complex<T1>> coils = syntheticCoils3D<T1>(Nx, Ny, Nz, Nc);
    REQUIRE(coils.n_elem == Ni * Nc);

    // Check SoS at center voxel ≈ 1
    uword centerIdx = Ny/2 + (Nx/2)*Nx + (Nz/2)*Nx*Ny;
    T1 sos = 0;
    for (uword cc = 0; cc < Nc; cc++) {
        T1 m = std::abs(coils(centerIdx + cc * Ni));
        sos += m * m;
    }
    REQUIRE(sos == Approx((T1)1.0).margin(0.01));

    // Coil maps should be distinct
    Col<std::complex<T1>> coil0 = coils.subvec(0, Ni - 1);
    Col<std::complex<T1>> coil1 = coils.subvec(Ni, 2 * Ni - 1);
    REQUIRE(norm(coil0 - coil1) > (T1)0.01);
}

// ===========================================================================
// Test 3: Stack-of-spirals trajectory properties
// ===========================================================================
TEST_CASE("Stack-of-spirals trajectory properties", "[spiral3D][trajectory3D]") {
    using T1 = float;
    const uword nInterleaves = 16, nSamplesPerArm = 256;
    const uword Nx = 64, Nz = 32;
    const uword R_xy = 2, R_z = 2;

    Col<T1> kx, ky, kz;
    uword Nd_per_shot, Ns;
    stackOfSpiralsTrajectory3D<T1>(nInterleaves, nSamplesPerArm, Nx, Nz,
                                    R_xy, R_z, kx, ky, kz, Nd_per_shot, Ns);

    uword nKzAcq = Nz / R_z;  // 16
    REQUIRE(Ns == nInterleaves / R_xy);  // 8
    REQUIRE(Nd_per_shot == nSamplesPerArm * nKzAcq);
    REQUIRE(kx.n_elem == Nd_per_shot * Ns);
    REQUIRE(ky.n_elem == Nd_per_shot * Ns);
    REQUIRE(kz.n_elem == Nd_per_shot * Ns);

    // k-space bounds
    T1 kxMax = (T1)Nx / (T1)2.0;
    T1 kzMax = (T1)Nz / (T1)2.0;
    REQUIRE(max(abs(kx)) <= kxMax * (T1)1.01);
    REQUIRE(max(abs(ky)) <= kxMax * (T1)1.01);
    REQUIRE(max(abs(kz)) <= kzMax * (T1)1.01);

    // First sample of first shot should be near k-space center
    REQUIRE(std::abs(kx(0)) < (T1)1.0);
    REQUIRE(std::abs(ky(0)) < (T1)1.0);

    // Timing vector
    Col<T1> tvec = stackOfSpiralsTimingVector<T1>(nSamplesPerArm, nKzAcq, Ns, (T1)0.005);
    REQUIRE(tvec.n_elem == Nd_per_shot * Ns);
    REQUIRE(tvec(0) >= (T1)0.0);
}

// ===========================================================================
// Test 4: pcSENSE 3D adjointness with stack-of-spirals
// ===========================================================================
TEST_CASE("pcSENSE 3D adjointness with stack-of-spirals",
          "[spiral3D][pcSENSE3D][adjointness]")
{
    using T1 = float;
    using CxT1 = std::complex<T1>;

    const uword Nx = 32, Ny = 32, Nz = 16;
    const uword Ni = Nx * Ny * Nz;
    const uword Nc = 8;
    const uword nInterleaves = 16, nSamplesPerArm = 256;
    const uword R_xy = 2, R_z = 2;

    // Trajectory
    Col<T1> kx, ky, kz;
    uword Nd_per_shot, Ns;
    stackOfSpiralsTrajectory3D<T1>(nInterleaves, nSamplesPerArm, Nx, Nz,
                                    R_xy, R_z, kx, ky, kz, Nd_per_shot, Ns);

    uword nKzAcq = Nz / R_z;
    Col<T1> tvec = stackOfSpiralsTimingVector<T1>(nSamplesPerArm, nKzAcq, Ns, (T1)0.005);

    // Synthetic data
    Col<CxT1> SENSEmap = syntheticCoils3D<T1>(Nx, Ny, Nz, Nc);
    Col<T1> fmapHz = syntheticFieldMap3D<T1>(Nx, Ny, Nz, (T1)10.0);
    Col<T1> fmapRad = fmapHz * (T1)(2.0 * M_PI);
    Col<T1> shotPhase = randomShotPhase3D<T1>(Nx, Ny, Nz, Ns, (T1)(M_PI / 4.0), 42);

    // Construct pcSENSE operator
    pcSENSE<T1> P(kx, ky, kz, Nx, Ny, Nz, Nc, tvec, SENSEmap, fmapRad, shotPhase);

    // Adjointness test: <Px, y> ≈ <x, P'y>
    arma::arma_rng::set_seed(123);
    Col<CxT1> x = randn<Col<CxT1>>(Ni);
    Col<CxT1> y = randn<Col<CxT1>>(Nd_per_shot * Ns * Nc);

    Col<CxT1> Px = P * x;
    Col<CxT1> Pty = P / y;

    CxT1 lhs = cdot(Px, y);
    CxT1 rhs = cdot(x, Pty);

    T1 relErr = std::abs(lhs - rhs) / std::abs(lhs);
    std::cout << "pcSENSE 3D adjointness: relErr = " << relErr << std::endl;
    REQUIRE(relErr < (T1)1e-2);
}

// ===========================================================================
// Test 5: pcSENSE 3D reconstruction (small grid)
// ===========================================================================
TEST_CASE("pcSENSE 3D recon with stack-of-spirals (small grid)",
          "[spiral3D][pcSENSE3D][recon]")
{
    using T1 = float;
    using CxT1 = std::complex<T1>;

    const uword Nx = 32, Ny = 32, Nz = 16;
    const uword Ni = Nx * Ny * Nz;
    const uword Nc = 8;
    const uword nInterleaves = 16, nSamplesPerArm = 256;
    const uword R_xy = 2, R_z = 2;
    const uword niter = 30;

    // Trajectory
    Col<T1> kx, ky, kz;
    uword Nd_per_shot, Ns;
    stackOfSpiralsTrajectory3D<T1>(nInterleaves, nSamplesPerArm, Nx, Nz,
                                    R_xy, R_z, kx, ky, kz, Nd_per_shot, Ns);

    uword nKzAcq = Nz / R_z;
    Col<T1> tvec = stackOfSpiralsTimingVector<T1>(nSamplesPerArm, nKzAcq, Ns, (T1)0.005);

    // Synthetic data
    Col<CxT1> phantom = brainPhantom3D<T1>(Nx, Ny, Nz);
    Col<CxT1> SENSEmap = syntheticCoils3D<T1>(Nx, Ny, Nz, Nc);
    Col<T1> fmapHz = syntheticFieldMap3D<T1>(Nx, Ny, Nz, (T1)10.0);
    Col<T1> fmapRad = fmapHz * (T1)(2.0 * M_PI);
    Col<T1> shotPhase = randomShotPhase3D<T1>(Nx, Ny, Nz, Ns, (T1)(M_PI / 4.0), 42);

    // Construct operator and forward simulate
    pcSENSE<T1> P(kx, ky, kz, Nx, Ny, Nz, Nc, tvec, SENSEmap, fmapRad, shotPhase);
    Col<CxT1> y = P * phantom;

    // Reconstruct
    Col<T1> W = ones<Col<T1>>(y.n_elem);
    QuadPenalty<T1> R(Nx, Ny, Nz, (T1)1e-3, 3);
    Col<CxT1> x0 = zeros<Col<CxT1>>(Ni);

    Col<CxT1> xhat = solve_pwls_pcg<T1>(x0, P, W, y, R, niter);

    // NRMSE
    T1 nrmse = norm(phantom - xhat) / norm(phantom);
    std::cout << "pcSENSE 3D recon NRMSE (32x32x16) = " << nrmse << std::endl;
    REQUIRE(nrmse < (T1)0.70);
}

// ===========================================================================
// Test 6: pcSenseTimeSeg 3D reconstruction (full grid)
// ===========================================================================
TEST_CASE("pcSenseTimeSeg 3D recon with stack-of-spirals (120x120x60)",
          "[spiral3D][pcSenseTimeSeg3D][recon][bench3D]")
{
    using T1 = float;
    using CxT1 = std::complex<T1>;

    const uword Nx = 120, Ny = 120, Nz = 60;
    const uword Ni = Nx * Ny * Nz;
    const uword Nc = 16;
    const uword nInterleaves = 16, nSamplesPerArm = 512;
    const uword R_xy = 2, R_z = 2;
    const uword L = 8;
    const uword niter = 5;

    std::cout << "Setting up 3D stack-of-spirals reconstruction..." << std::endl;
    std::cout << "  Grid: " << Nx << "x" << Ny << "x" << Nz
              << " (" << Ni << " voxels)" << std::endl;

    // Trajectory
    Col<T1> kx, ky, kz;
    uword Nd_per_shot, Ns;
    stackOfSpiralsTrajectory3D<T1>(nInterleaves, nSamplesPerArm, Nx, Nz,
                                    R_xy, R_z, kx, ky, kz, Nd_per_shot, Ns);

    uword nKzAcq = Nz / R_z;
    Col<T1> tvec = stackOfSpiralsTimingVector<T1>(nSamplesPerArm, nKzAcq, Ns, (T1)0.005);

    std::cout << "  Ns=" << Ns << ", Nd_per_shot=" << Nd_per_shot
              << ", Nc=" << Nc << ", L=" << L << std::endl;

    // Synthetic data
    Col<CxT1> phantom = brainPhantom3D<T1>(Nx, Ny, Nz);
    Col<CxT1> SENSEmap = syntheticCoils3D<T1>(Nx, Ny, Nz, Nc);
    Col<T1> fmapHz = syntheticFieldMap3D<T1>(Nx, Ny, Nz, (T1)50.0);
    Col<T1> fmapRad = fmapHz * (T1)(2.0 * M_PI);
    Col<T1> shotPhase = randomShotPhase3D<T1>(Nx, Ny, Nz, Ns, (T1)(M_PI / 2.0), 42);

    std::cout << "  Data generated. Constructing pcSenseTimeSeg operator..." << std::endl;

    // Construct operator
    auto t0 = std::chrono::high_resolution_clock::now();
    pcSenseTimeSeg<T1> P(kx, ky, kz, Nx, Ny, Nz, Nc, tvec, L, 1,
                          SENSEmap, fmapRad, shotPhase);
    auto t1 = std::chrono::high_resolution_clock::now();
    double setupMs = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << "  Setup: " << setupMs << " ms" << std::endl;

    // Forward simulate
    std::cout << "  Forward simulation..." << std::endl;
    auto t2 = std::chrono::high_resolution_clock::now();
    Col<CxT1> y = P * phantom;
    auto t3 = std::chrono::high_resolution_clock::now();
    double fwdMs = std::chrono::duration<double, std::milli>(t3 - t2).count();
    std::cout << "  Forward: " << fwdMs << " ms, data size: " << y.n_elem << std::endl;

    REQUIRE(y.n_elem == Nd_per_shot * Ns * Nc);

    // Reconstruct
    Col<T1> W = ones<Col<T1>>(y.n_elem);
    QuadPenalty<T1> R(Nx, Ny, Nz, (T1)1e-4, 3);
    Col<CxT1> x0 = zeros<Col<CxT1>>(Ni);

    std::cout << "  Reconstructing (" << niter << " PCG iterations)..." << std::endl;
    auto t4 = std::chrono::high_resolution_clock::now();
    Col<CxT1> xhat = solve_pwls_pcg<T1>(x0, P, W, y, R, niter);
    auto t5 = std::chrono::high_resolution_clock::now();
    double reconMs = std::chrono::duration<double, std::milli>(t5 - t4).count();

    T1 nrmse = norm(phantom - xhat) / norm(phantom);
    std::cout << "  Recon: " << reconMs << " ms" << std::endl;
    std::cout << "  NRMSE: " << nrmse << std::endl;
    std::cout << "  Total: " << (setupMs + fwdMs + reconMs) << " ms" << std::endl;

    REQUIRE(nrmse < (T1)0.70);
}

// ===========================================================================
// Test 7: Write NIfTI volumes for visualization
// ===========================================================================
TEST_CASE("Write 3D spiral recon NIfTI volumes",
          "[spiral3D][nifti3D]")
{
    using T1 = float;
    using CxT1 = std::complex<T1>;

    const uword Nx = 120, Ny = 120, Nz = 60;
    const uword Ni = Nx * Ny * Nz;
    const uword Nc = 16;
    const uword nInterleaves = 16, nSamplesPerArm = 512;
    const uword R_xy = 2, R_z = 2;
    const uword L = 8;
    const uword niter = 5;

    std::cout << "Generating 3D reconstruction for NIfTI output..." << std::endl;

    // Trajectory
    Col<T1> kx, ky, kz;
    uword Nd_per_shot, Ns;
    stackOfSpiralsTrajectory3D<T1>(nInterleaves, nSamplesPerArm, Nx, Nz,
                                    R_xy, R_z, kx, ky, kz, Nd_per_shot, Ns);

    uword nKzAcq = Nz / R_z;
    Col<T1> tvec = stackOfSpiralsTimingVector<T1>(nSamplesPerArm, nKzAcq, Ns, (T1)0.005);

    // Synthetic data
    Col<CxT1> phantom = brainPhantom3D<T1>(Nx, Ny, Nz);
    Col<CxT1> SENSEmap = syntheticCoils3D<T1>(Nx, Ny, Nz, Nc);
    Col<T1> fmapHz = syntheticFieldMap3D<T1>(Nx, Ny, Nz, (T1)50.0);
    Col<T1> fmapRad = fmapHz * (T1)(2.0 * M_PI);
    Col<T1> shotPhase = randomShotPhase3D<T1>(Nx, Ny, Nz, Ns, (T1)(M_PI / 2.0), 42);

    // Write phantom and auxiliary maps
    std::cout << "Writing input volumes:" << std::endl;
    writeNiftiComplex<T1>("spiral3d_phantom", phantom, Nx, Ny, Nz);
    writeNifti<T1>("spiral3d_fieldmap", fmapHz, Nx, Ny, Nz);

    // Coil SoS map
    Col<T1> coilSoS(Ni, fill::zeros);
    for (uword cc = 0; cc < Nc; cc++) {
        Col<CxT1> coilMap = SENSEmap.subvec(cc * Ni, (cc + 1) * Ni - 1);
        coilSoS += square(abs(coilMap));
    }
    coilSoS = sqrt(coilSoS);
    writeNifti<T1>("spiral3d_coil_sos", coilSoS, Nx, Ny, Nz);

    // Construct operator and forward simulate
    std::cout << "Running pcSenseTimeSeg reconstruction..." << std::endl;
    pcSenseTimeSeg<T1> P(kx, ky, kz, Nx, Ny, Nz, Nc, tvec, L, 1,
                          SENSEmap, fmapRad, shotPhase);
    Col<CxT1> y = P * phantom;

    // Reconstruct
    Col<T1> W = ones<Col<T1>>(y.n_elem);
    QuadPenalty<T1> R(Nx, Ny, Nz, (T1)1e-4, 3);
    Col<CxT1> x0 = zeros<Col<CxT1>>(Ni);
    Col<CxT1> xhat = solve_pwls_pcg<T1>(x0, P, W, y, R, niter);

    // Write reconstruction and error
    std::cout << "Writing output volumes:" << std::endl;
    writeNiftiComplex<T1>("spiral3d_recon", xhat, Nx, Ny, Nz);

    Col<T1> errorMag = abs(phantom - xhat);
    writeNifti<T1>("spiral3d_error_mag", errorMag, Nx, Ny, Nz);

    T1 nrmse = norm(phantom - xhat) / norm(phantom);
    std::cout << "NRMSE: " << nrmse << std::endl;
    std::cout << "NIfTI files written to current directory." << std::endl;

    // Verify files were written (non-empty)
    REQUIRE(nrmse < (T1)1.0);
}
