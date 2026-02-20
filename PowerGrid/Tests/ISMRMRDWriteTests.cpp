/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file ISMRMRDWriteTests.cpp
/// @brief Tests that write synthetic ISMRMRD files, read them back,
///        reconstruct, and verify against ground truth.

#include "catch.hpp"

#include "../PowerGrid.h"
#include "../processIsmrmrd.hpp"
#include "../reconSolve.h"

#include "SyntheticPhantom.hpp"
#include "SyntheticCoils.hpp"
#include "SyntheticTrajectory.hpp"
#include "SyntheticFieldMap.hpp"

#include <chrono>
#include <sstream>
#include <cstdio>
#include "ismrmrd/version.h"

using namespace arma;

// ---------------------------------------------------------------------------
// Helper: compute image-space coordinates matching initImageSpaceCoords()
// Uses Cube(Nx, Ny, Nz) with loops: outer y (ii=0..Ny-1), inner x (jj=0..Nx-1)
// ---------------------------------------------------------------------------
static void imageCoords(uword Nx, uword Ny, uword Nz,
                        Col<float>& ix, Col<float>& iy, Col<float>& iz)
{
    initImageSpaceCoords<float>(ix, iy, iz, Nx, Ny, Nz);
}

// ---------------------------------------------------------------------------
// Helper: NRMSE
// ---------------------------------------------------------------------------
template<typename T1>
static T1 nrmse(const Col<std::complex<T1>>& x, const Col<std::complex<T1>>& ref) {
    T1 errNorm = norm(x - ref, 2);
    T1 refNorm = norm(ref, 2);
    if (refNorm < (T1)1e-20) return errNorm;
    return errNorm / refNorm;
}

// ---------------------------------------------------------------------------
// writeSyntheticISMRMRD — write a complete ISMRMRD HDF5 file from synthetic data
// ---------------------------------------------------------------------------
static void writeSyntheticISMRMRD(
    const std::string& filename,
    uword Nx, uword Ny, uword Nz,
    uword Nc,
    uword nSpokes, uword nReadout,
    const Col<float>& kx, const Col<float>& ky, const Col<float>& kz,
    const Col<float>& tvec,
    const Col<std::complex<float>>& kspaceData,
    const Col<std::complex<float>>& SENSEmap,
    const Col<float>& fieldMap,
    const Col<float>* phaseMaps = nullptr)
{
    uword Nd = nSpokes * nReadout;
    uword Ni = Nx * Ny * Nz;

    // Remove existing file to avoid appending (ISMRMRD doesn't truncate)
    std::remove(filename.c_str());

    // Create the HDF5 dataset
    ISMRMRD::Dataset dataset(filename.c_str(), "dataset", true);

    // Build XML header
    ISMRMRD::IsmrmrdHeader hdr;
    hdr.version = ISMRMRD_XMLHDR_VERSION;

    // Experimental conditions (required)
    hdr.experimentalConditions.H1resonanceFrequency_Hz = 128000000; // 3T

    // Acquisition system info
    ISMRMRD::AcquisitionSystemInformation sysInfo;
    sysInfo.receiverChannels = (unsigned short)Nc;
    hdr.acquisitionSystemInformation = sysInfo;

    // Encoding
    ISMRMRD::Encoding enc;
    enc.encodedSpace.matrixSize = ISMRMRD::MatrixSize((unsigned short)Nx,
                                                       (unsigned short)Ny,
                                                       (unsigned short)Nz);
    enc.encodedSpace.fieldOfView_mm.x = 256.0f;
    enc.encodedSpace.fieldOfView_mm.y = 256.0f;
    enc.encodedSpace.fieldOfView_mm.z = 5.0f;

    enc.reconSpace.matrixSize = ISMRMRD::MatrixSize((unsigned short)Nx,
                                                     (unsigned short)Ny,
                                                     (unsigned short)Nz);
    enc.reconSpace.fieldOfView_mm.x = 256.0f;
    enc.reconSpace.fieldOfView_mm.y = 256.0f;
    enc.reconSpace.fieldOfView_mm.z = 5.0f;

    // ISMRMRD >= 1.14: enum class TrajectoryType; <= 1.3: std::string
#if ISMRMRD_VERSION_MINOR >= 14
    enc.trajectory = ISMRMRD::TrajectoryType::OTHER;
#else
    enc.trajectory = "other";
#endif

    // Encoding limits
    enc.encodingLimits.kspace_encoding_step_1 =
        ISMRMRD::Limit(0, (unsigned short)(nSpokes - 1), 0);
    enc.encodingLimits.kspace_encoding_step_2 =
        ISMRMRD::Limit(0, 0, 0);
    enc.encodingLimits.slice =
        ISMRMRD::Limit(0, 0, 0);
    enc.encodingLimits.repetition =
        ISMRMRD::Limit(0, 0, 0);
    enc.encodingLimits.average =
        ISMRMRD::Limit(0, 0, 0);
    enc.encodingLimits.contrast =
        ISMRMRD::Limit(0, 0, 0);
    enc.encodingLimits.phase =
        ISMRMRD::Limit(0, 0, 0);
    enc.encodingLimits.segment =
        ISMRMRD::Limit(0, 0, 0);
    enc.encodingLimits.set =
        ISMRMRD::Limit(0, 0, 0);

    hdr.encoding.push_back(enc);

    // Serialize header to XML and write
    std::ostringstream oss;
    ISMRMRD::serialize(hdr, oss);
    dataset.writeHeader(oss.str());

    // Write acquisitions — one per spoke
    for (uword s = 0; s < nSpokes; s++) {
        ISMRMRD::Acquisition acq((uint16_t)nReadout, (uint16_t)Nc, (uint16_t)4);

        // Set encoding counters
        acq.idx().kspace_encode_step_1 = (uint16_t)s;
        acq.idx().kspace_encode_step_2 = 0;
        acq.idx().slice = 0;
        acq.idx().repetition = 0;
        acq.idx().average = 0;
        acq.idx().contrast = 0;
        acq.idx().phase = 0;
        acq.idx().segment = 0;
        acq.idx().set = 0;

        // Set trajectory (4D: kx, ky, kz, tvec)
        for (uword r = 0; r < nReadout; r++) {
            uword idx = s * nReadout + r;
            acq.traj(0, (uint16_t)r) = kx(idx);
            acq.traj(1, (uint16_t)r) = ky(idx);
            acq.traj(2, (uint16_t)r) = kz(idx);
            acq.traj(3, (uint16_t)r) = tvec(idx);
        }

        // Set data: kspaceData is vectorise(Mat(Nd, Nc)) = coil 0 first Nd, coil 1 next, etc.
        for (uword c = 0; c < Nc; c++) {
            for (uword r = 0; r < nReadout; r++) {
                uword dataIdx = s * nReadout + r + c * Nd;
                std::complex<float> val = kspaceData(dataIdx);
                acq.data((uint16_t)r, (uint16_t)c) = {val.real(), val.imag()};
            }
        }

        dataset.appendAcquisition(acq);
    }

    // Write SENSEMap NDArray (as complex<double>)
    {
        std::vector<size_t> dims;
        dims.push_back(Ni * Nc);
        ISMRMRD::NDArray<std::complex<double>> senArray(dims);
        for (uword i = 0; i < Ni * Nc; i++) {
            senArray.getDataPtr()[i] = std::complex<double>(
                (double)SENSEmap(i).real(), (double)SENSEmap(i).imag());
        }
        dataset.appendNDArray("SENSEMap", senArray);
    }

    // Write FieldMap NDArray (as double)
    {
        std::vector<size_t> dims;
        dims.push_back(Ni);
        ISMRMRD::NDArray<double> fmArray(dims);
        for (uword i = 0; i < Ni; i++) {
            fmArray.getDataPtr()[i] = (double)fieldMap(i);
        }
        dataset.appendNDArray("FieldMap", fmArray);
    }

    // Write PhaseMaps NDArray if provided (as double)
    if (phaseMaps != nullptr) {
        uword pmSize = phaseMaps->n_elem;
        std::vector<size_t> dims;
        dims.push_back(pmSize);
        ISMRMRD::NDArray<double> pmArray(dims);
        for (uword i = 0; i < pmSize; i++) {
            pmArray.getDataPtr()[i] = (double)(*phaseMaps)(i);
        }
        dataset.appendNDArray("PhaseMaps", pmArray);
    }
}

// ===========================================================================
// Test 1: Round-trip write/read verification (no reconstruction)
// ===========================================================================
TEST_CASE("ISMRMRD write and read-back round-trip", "[ISMRMRD][roundtrip]") {
    typedef float T1;
    typedef std::complex<T1> CxT1;

    const uword Nx = 32, Ny = 32, Nz = 1;
    const uword Ni = Nx * Ny * Nz;
    const uword Nc = 4;
    const uword nSpokes = 32, nReadout = 64;
    const uword Nd = nSpokes * nReadout;

    // Generate synthetic data
    Col<CxT1> phantom = sheppLogan2D<T1>(Nx, Ny);
    Col<CxT1> SENSEmap = syntheticCoils2D<T1>(Nx, Ny, Nc);
    Col<T1> kx, ky, kz;
    radialTrajectory2D<T1>(nSpokes, nReadout, Nx, kx, ky, kz);
    Col<T1> fieldMap;
    fieldMap.zeros(Ni);
    Col<T1> tvec = linearTimingVector<T1>(Nd, (T1)0.01);

    // Image coords
    Col<T1> ix, iy, iz;
    imageCoords(Nx, Ny, Nz, ix, iy, iz);

    // Generate k-space data via SENSE forward model
    Gnufft<T1> G(Nd, (T1)2.0, Nx, Ny, Nz, kx, ky, kz, ix, iy, iz);
    SENSE<T1, Gnufft<T1>> S(G, SENSEmap, Nd, Ni, Nc);
    Col<CxT1> kspaceData = S * phantom;

    // Write ISMRMRD file
    std::string filename = "/tmp/pg_test_roundtrip.h5";
    writeSyntheticISMRMRD(filename, Nx, Ny, Nz, Nc, nSpokes, nReadout,
                          kx, ky, kz, tvec, kspaceData, SENSEmap, fieldMap);

    // Read back
    ISMRMRD::Dataset* d = nullptr;
    ISMRMRD::IsmrmrdHeader hdr;
    acqTracking* acqTrack = nullptr;
    openISMRMRDData(filename, d, hdr, acqTrack);

    REQUIRE(hdr.encoding.size() == 1);
    REQUIRE(hdr.encoding[0].encodedSpace.matrixSize.x == Nx);
    REQUIRE(hdr.encoding[0].encodedSpace.matrixSize.y == Ny);
    REQUIRE(hdr.encoding[0].encodedSpace.matrixSize.z == Nz);
    REQUIRE(hdr.encoding[0].encodingLimits.kspace_encoding_step_1->maximum == nSpokes - 1);

    // Read SENSEMap
    Col<CxT1> senReadBack = getISMRMRDSenseMap<CxT1>(d);
    REQUIRE(senReadBack.n_elem == Ni * Nc);
    T1 senErr = norm(senReadBack - SENSEmap, 2) / norm(SENSEmap, 2);
    std::cout << "SENSEMap round-trip relative error: " << senErr << std::endl;
    REQUIRE(senErr < (T1)1e-5);

    // Read FieldMap
    Col<T1> fmReadBack = getISMRMRDFieldMap<T1>(d);
    REQUIRE(fmReadBack.n_elem == Ni);
    // Field map is zero, just check it's all zeros
    REQUIRE(norm(fmReadBack, 2) < (T1)1e-10);

    // Read acquisition data
    Col<CxT1> dataReadBack;
    Col<T1> kxRB, kyRB, kzRB, tvecRB;
    getCompleteISMRMRDAcqData<T1>(d, acqTrack, 0, 0, 0, 0, 0,
                                   dataReadBack, kxRB, kyRB, kzRB, tvecRB);

    REQUIRE(kxRB.n_elem == Nd);
    REQUIRE(kyRB.n_elem == Nd);
    REQUIRE(dataReadBack.n_elem == Nd * Nc);

    T1 kxErr = norm(kxRB - kx, 2) / norm(kx, 2);
    T1 kyErr = norm(kyRB - ky, 2) / norm(ky, 2);
    std::cout << "kx round-trip error: " << kxErr << std::endl;
    std::cout << "ky round-trip error: " << kyErr << std::endl;
    REQUIRE(kxErr < (T1)1e-5);
    REQUIRE(kyErr < (T1)1e-5);

    T1 dataErr = norm(dataReadBack - kspaceData, 2) / norm(kspaceData, 2);
    std::cout << "k-space data round-trip error: " << dataErr << std::endl;
    REQUIRE(dataErr < (T1)1e-5);

    closeISMRMRDData(d, hdr, acqTrack);
    std::remove(filename.c_str());
}

// ===========================================================================
// Test 2: SENSE + NUFFT reconstruction from ISMRMRD file
// ===========================================================================
TEST_CASE("ISMRMRD SENSE reconstruction", "[ISMRMRD][recon][SENSE]") {
    typedef float T1;
    typedef std::complex<T1> CxT1;

    const uword Nx = 64, Ny = 64, Nz = 1;
    const uword Ni = Nx * Ny * Nz;
    const uword Nc = 4;
    const uword nSpokes = 128, nReadout = 128;
    const uword Nd = nSpokes * nReadout;
    const uword niter = 10;

    std::cout << "\n=== ISMRMRD SENSE Reconstruction Test ===" << std::endl;
    std::cout << "Matrix: " << Nx << "x" << Ny << ", Coils: " << Nc
              << ", Spokes: " << nSpokes << ", Readout: " << nReadout << std::endl;

    // Generate synthetic data
    Col<CxT1> phantom = sheppLogan2D<T1>(Nx, Ny);
    Col<CxT1> SENSEmap = syntheticCoils2D<T1>(Nx, Ny, Nc);
    Col<T1> kx, ky, kz;
    radialTrajectory2D<T1>(nSpokes, nReadout, Nx, kx, ky, kz);
    Col<T1> fieldMap;
    fieldMap.zeros(Ni);
    Col<T1> tvec = linearTimingVector<T1>(Nd, (T1)0.01);

    Col<T1> ix, iy, iz;
    imageCoords(Nx, Ny, Nz, ix, iy, iz);

    // Forward model to generate k-space data
    Gnufft<T1> Gfwd(Nd, (T1)2.0, Nx, Ny, Nz, kx, ky, kz, ix, iy, iz);
    SENSE<T1, Gnufft<T1>> Sfwd(Gfwd, SENSEmap, Nd, Ni, Nc);
    Col<CxT1> kspaceData = Sfwd * phantom;

    // Write ISMRMRD file
    std::string filename = "/tmp/pg_test_sense_recon.h5";
    writeSyntheticISMRMRD(filename, Nx, Ny, Nz, Nc, nSpokes, nReadout,
                          kx, ky, kz, tvec, kspaceData, SENSEmap, fieldMap);

    // ---- Read back (same as PowerGridIsmrmrd.cpp) ----
    ISMRMRD::Dataset* d = nullptr;
    ISMRMRD::IsmrmrdHeader hdr;
    acqTracking* acqTrack = nullptr;
    openISMRMRDData(filename, d, hdr, acqTrack);

    Col<CxT1> senSlice = getISMRMRDSenseMap<CxT1>(d);
    Col<T1> fmSlice = getISMRMRDFieldMap<T1>(d);

    Col<CxT1> data;
    Col<T1> kxR, kyR, kzR, tvecR;
    getCompleteISMRMRDAcqData<T1>(d, acqTrack, 0, 0, 0, 0, 0,
                                   data, kxR, kyR, kzR, tvecR);

    uword nroR = d->getNumberOfAcquisitions() > 0 ? nReadout : 0;
    REQUIRE(data.n_elem == Nd * Nc);

    Col<T1> ixR, iyR, izR;
    imageCoords(Nx, Ny, Nz, ixR, iyR, izR);

    // ---- Reconstruct (same pipeline as PowerGridIsmrmrd -F NUFFT -t 1) ----
    auto t_setup_start = std::chrono::high_resolution_clock::now();

    Gnufft<T1> Grecon(kxR.n_rows, (T1)2.0, Nx, Ny, Nz, kxR, kyR, kzR, ixR, iyR, izR);
    TimeSegmentation<T1, Gnufft<T1>> Arecon(Grecon, fmSlice, tvecR,
                                             kxR.n_rows, Ni, 1, 1, 1);
    SENSE<T1, TimeSegmentation<T1, Gnufft<T1>>> Sgrecon(Arecon, senSlice,
                                                         kxR.n_rows, Ni, Nc);
    QuadPenalty<T1> R(Nx, Ny, Nz, (T1)1e-3, 2);

    auto t_setup_end = std::chrono::high_resolution_clock::now();

    auto t_recon_start = std::chrono::high_resolution_clock::now();
    Col<CxT1> xhat = reconSolve<T1>(data, Sgrecon, R,
                                      kxR, kyR, kzR, Nx, Ny, Nz, tvecR, niter);
    auto t_recon_end = std::chrono::high_resolution_clock::now();

    double setup_ms = std::chrono::duration<double, std::milli>(t_setup_end - t_setup_start).count();
    double recon_ms = std::chrono::duration<double, std::milli>(t_recon_end - t_recon_start).count();

    REQUIRE(xhat.n_elem == Ni);
    REQUIRE(!xhat.has_nan());

    T1 err = nrmse<T1>(xhat, phantom);
    std::cout << "Setup time:  " << setup_ms << " ms" << std::endl;
    std::cout << "Recon time:  " << recon_ms << " ms (" << niter << " iterations)" << std::endl;
    std::cout << "NRMSE:       " << err << std::endl;
    REQUIRE(err < (T1)0.50);

    closeISMRMRDData(d, hdr, acqTrack);
    std::remove(filename.c_str());
}

// ===========================================================================
// Test 3: pcSENSE reconstruction from ISMRMRD file
// ===========================================================================
TEST_CASE("ISMRMRD pcSENSE reconstruction", "[ISMRMRD][recon][pcSENSE]") {
    typedef float T1;
    typedef std::complex<T1> CxT1;

    const uword Nx = 64, Ny = 64, Nz = 1;
    const uword Ni = Nx * Ny * Nz;
    const uword Nc = 4;
    const uword Ns = 2; // shots
    const uword nSpokesPerShot = 64;
    const uword nReadout = 64;
    const uword nSpokesTotal = Ns * nSpokesPerShot;
    const uword NdPerShot = nSpokesPerShot * nReadout;
    const uword NdTotal = nSpokesTotal * nReadout;
    const uword niter = 10;

    std::cout << "\n=== ISMRMRD pcSENSE Reconstruction Test ===" << std::endl;
    std::cout << "Matrix: " << Nx << "x" << Ny << ", Coils: " << Nc
              << ", Shots: " << Ns << ", Spokes/shot: " << nSpokesPerShot << std::endl;

    // Generate synthetic data
    Col<CxT1> phantom = sheppLogan2D<T1>(Nx, Ny);
    Col<CxT1> SENSEmap = syntheticCoils2D<T1>(Nx, Ny, Nc);

    // Trajectory: all shots concatenated
    Col<T1> kx, ky, kz;
    radialTrajectory2D<T1>(nSpokesTotal, nReadout, Nx, kx, ky, kz);

    Col<T1> fieldMap;
    fieldMap.zeros(Ni);
    Col<T1> tvec = linearTimingVector<T1>(NdTotal, (T1)0.01);

    // Phase maps: simple linear ramp per shot (Ni * Ns elements)
    Col<T1> phaseMaps(Ni * Ns);
    for (uword s = 0; s < Ns; s++) {
        for (uword i = 0; i < Ni; i++) {
            // Linear phase ramp: shot 0 = 0, shot 1 = small gradient
            phaseMaps(i + s * Ni) = (T1)s * (T1)0.1 * ((T1)(i % Nx) / (T1)Nx);
        }
    }

    Col<T1> ix, iy, iz;
    imageCoords(Nx, Ny, Nz, ix, iy, iz);

    // Forward model via pcSENSE
    pcSENSE<T1> Sfwd(kx, ky, kz, Nx, Ny, Nz, Nc, tvec, SENSEmap, fieldMap, phaseMaps);
    Col<CxT1> kspaceData = Sfwd * phantom;

    // Write ISMRMRD file
    std::string filename = "/tmp/pg_test_pcsense_recon.h5";
    writeSyntheticISMRMRD(filename, Nx, Ny, Nz, Nc, nSpokesTotal, nReadout,
                          kx, ky, kz, tvec, kspaceData, SENSEmap, fieldMap, &phaseMaps);

    // Read back
    ISMRMRD::Dataset* d = nullptr;
    ISMRMRD::IsmrmrdHeader hdr;
    acqTracking* acqTrack = nullptr;
    openISMRMRDData(filename, d, hdr, acqTrack);

    Col<CxT1> senRB = getISMRMRDSenseMap<CxT1>(d);
    Col<T1> fmRB = getISMRMRDFieldMap<T1>(d);
    Col<T1> pmRB = getISMRMRDPhaseMaps<T1>(d);

    Col<CxT1> dataRB;
    Col<T1> kxRB, kyRB, kzRB, tvecRB;
    getCompleteISMRMRDAcqData<T1>(d, acqTrack, 0, 0, 0, 0, 0,
                                   dataRB, kxRB, kyRB, kzRB, tvecRB);

    // Reconstruct using pcSENSE (same pipeline as PowerGridPcSense)
    auto t_start = std::chrono::high_resolution_clock::now();

    pcSENSE<T1> Srecon(kxRB, kyRB, kzRB, Nx, Ny, Nz, Nc, tvecRB,
                        senRB, fmRB, (T1)0 - pmRB);
    QuadPenalty<T1> R(Nx, Ny, Nz, (T1)1e-3, 2);
    Col<CxT1> xhat = reconSolve<T1>(dataRB, Srecon, R,
                                      kxRB, kyRB, kzRB, Nx, Ny, Nz, tvecRB, niter);

    auto t_end = std::chrono::high_resolution_clock::now();
    double recon_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    REQUIRE(xhat.n_elem == Ni);
    REQUIRE(!xhat.has_nan());

    T1 err = nrmse<T1>(xhat, phantom);
    std::cout << "Recon time:  " << recon_ms << " ms (" << niter << " iterations)" << std::endl;
    std::cout << "NRMSE:       " << err << std::endl;
    REQUIRE(err < (T1)0.50);

    closeISMRMRDData(d, hdr, acqTrack);
    std::remove(filename.c_str());
}

// ===========================================================================
// Test 4: SENSE + TimeSegmentation with field map from ISMRMRD file
// ===========================================================================
TEST_CASE("ISMRMRD SENSE+TimeSeg reconstruction with field map", "[ISMRMRD][recon][TimeSeg]") {
    typedef float T1;
    typedef std::complex<T1> CxT1;

    const uword Nx = 64, Ny = 64, Nz = 1;
    const uword Ni = Nx * Ny * Nz;
    const uword Nc = 4;
    const uword nSpokes = 128, nReadout = 128;
    const uword Nd = nSpokes * nReadout;
    const uword niter = 10;
    const uword L = 4;

    std::cout << "\n=== ISMRMRD SENSE+TimeSeg Reconstruction Test ===" << std::endl;

    Col<CxT1> phantom = sheppLogan2D<T1>(Nx, Ny);
    Col<CxT1> SENSEmap = syntheticCoils2D<T1>(Nx, Ny, Nc);
    Col<T1> kx, ky, kz;
    radialTrajectory2D<T1>(nSpokes, nReadout, Nx, kx, ky, kz);
    // syntheticFieldMap2D returns Hz; TimeSegmentation expects rad/s.
    // PowerGridIsmrmrd passes FM from the file directly, so the file should store rad/s.
    Col<T1> fieldMapHz = syntheticFieldMap2D<T1>(Nx, Ny, (T1)50.0);
    Col<T1> fieldMap = fieldMapHz * (T1)(2.0 * M_PI); // rad/s
    Col<T1> tvec = linearTimingVector<T1>(Nd, (T1)0.01);

    Col<T1> ix, iy, iz;
    imageCoords(Nx, Ny, Nz, ix, iy, iz);

    // Forward model with TimeSegmentation (field map in rad/s)
    Gnufft<T1> Gfwd(Nd, (T1)2.0, Nx, Ny, Nz, kx, ky, kz, ix, iy, iz);
    TimeSegmentation<T1, Gnufft<T1>> Afwd(Gfwd, fieldMap, tvec, Nd, Ni, L, 1, 1);
    SENSE<T1, TimeSegmentation<T1, Gnufft<T1>>> Sfwd(Afwd, SENSEmap, Nd, Ni, Nc);
    Col<CxT1> kspaceData = Sfwd * phantom;

    // Write ISMRMRD file (store field map in rad/s — matches PowerGridIsmrmrd convention)
    std::string filename = "/tmp/pg_test_timeseg_recon.h5";
    writeSyntheticISMRMRD(filename, Nx, Ny, Nz, Nc, nSpokes, nReadout,
                          kx, ky, kz, tvec, kspaceData, SENSEmap, fieldMap);

    // Read back
    ISMRMRD::Dataset* d = nullptr;
    ISMRMRD::IsmrmrdHeader hdr;
    acqTracking* acqTrack = nullptr;
    openISMRMRDData(filename, d, hdr, acqTrack);

    Col<CxT1> senRB = getISMRMRDSenseMap<CxT1>(d);
    Col<T1> fmRB = getISMRMRDFieldMap<T1>(d);  // rad/s
    Col<CxT1> dataRB;
    Col<T1> kxRB, kyRB, kzRB, tvecRB;
    getCompleteISMRMRDAcqData<T1>(d, acqTrack, 0, 0, 0, 0, 0,
                                   dataRB, kxRB, kyRB, kzRB, tvecRB);

    // Verify field map round-trip
    REQUIRE(fmRB.n_elem == Ni);
    T1 fmErr = norm(fmRB - fieldMap, 2) / norm(fieldMap, 2);
    std::cout << "Field map round-trip error: " << fmErr << std::endl;
    REQUIRE(fmErr < (T1)1e-5);
    REQUIRE(dataRB.n_elem == Nd * Nc);

    // Reconstruct (same pipeline as PowerGridIsmrmrd -F NUFFT -t L)
    auto t_start = std::chrono::high_resolution_clock::now();

    Gnufft<T1> Grecon(kxRB.n_rows, (T1)2.0, Nx, Ny, Nz, kxRB, kyRB, kzRB, ix, iy, iz);
    TimeSegmentation<T1, Gnufft<T1>> Arecon(Grecon, fmRB, tvecRB,
                                             kxRB.n_rows, Ni, L, 1, 1);
    SENSE<T1, TimeSegmentation<T1, Gnufft<T1>>> Sgrecon(Arecon, senRB,
                                                         kxRB.n_rows, Ni, Nc);
    QuadPenalty<T1> R(Nx, Ny, Nz, (T1)1e-3, 2);
    Col<CxT1> xhat = reconSolve<T1>(dataRB, Sgrecon, R,
                                      kxRB, kyRB, kzRB, Nx, Ny, Nz, tvecRB, niter);

    auto t_end = std::chrono::high_resolution_clock::now();
    double recon_ms = std::chrono::duration<double, std::milli>(t_end - t_start).count();

    REQUIRE(xhat.n_elem == Ni);

    T1 err = nrmse<T1>(xhat, phantom);
    std::cout << "Recon time (L=" << L << "): " << recon_ms << " ms" << std::endl;
    std::cout << "NRMSE:       " << err << std::endl;

    REQUIRE(err < (T1)0.50);

    closeISMRMRDData(d, hdr, acqTrack);
    std::remove(filename.c_str());
}

// ===========================================================================
// Test 5: In-process benchmark across matrix sizes
// ===========================================================================
TEST_CASE("ISMRMRD SENSE reconstruction benchmark", "[ISMRMRD][benchmark]") {
    typedef float T1;
    typedef std::complex<T1> CxT1;

    struct BenchConfig {
        uword Nx, Ny, Nc, nSpokes, nReadout, niter;
    };

    std::vector<BenchConfig> configs = {
        {32,  32,  4, 64,   64,  10},
        {64,  64,  4, 128,  128, 10},
        {128, 128, 8, 256,  256, 10},
    };

    std::cout << "\n=== ISMRMRD SENSE Reconstruction Benchmark ===" << std::endl;
    printf("%-10s %-6s %-7s %-5s %-6s %-12s %-12s %-8s\n",
           "Size", "Coils", "Spokes", "Nro", "Iters", "Setup(ms)", "Recon(ms)", "NRMSE");
    printf("%-10s %-6s %-7s %-5s %-6s %-12s %-12s %-8s\n",
           "----", "-----", "------", "---", "-----", "---------", "---------", "-----");

    for (const auto& cfg : configs) {
        uword Nx = cfg.Nx, Ny = cfg.Ny, Nz = 1;
        uword Ni = Nx * Ny * Nz;
        uword Nc = cfg.Nc;
        uword nSpokes = cfg.nSpokes, nReadout = cfg.nReadout;
        uword Nd = nSpokes * nReadout;

        Col<CxT1> phantom = sheppLogan2D<T1>(Nx, Ny);
        Col<CxT1> SENSEmap = syntheticCoils2D<T1>(Nx, Ny, Nc);
        Col<T1> kx, ky, kz;
        radialTrajectory2D<T1>(nSpokes, nReadout, Nx, kx, ky, kz);
        Col<T1> fieldMap;
        fieldMap.zeros(Ni);
        Col<T1> tvec = linearTimingVector<T1>(Nd, (T1)0.01);

        Col<T1> ix, iy, iz;
        imageCoords(Nx, Ny, Nz, ix, iy, iz);

        Gnufft<T1> Gfwd(Nd, (T1)2.0, Nx, Ny, Nz, kx, ky, kz, ix, iy, iz);
        SENSE<T1, Gnufft<T1>> Sfwd(Gfwd, SENSEmap, Nd, Ni, Nc);
        Col<CxT1> kspaceData = Sfwd * phantom;

        std::string filename = "/tmp/pg_bench_" + std::to_string(Nx) + ".h5";
        writeSyntheticISMRMRD(filename, Nx, Ny, Nz, Nc, nSpokes, nReadout,
                              kx, ky, kz, tvec, kspaceData, SENSEmap, fieldMap);

        // Read back
        ISMRMRD::Dataset* d = nullptr;
        ISMRMRD::IsmrmrdHeader hdr;
        acqTracking* acqTrack = nullptr;
        openISMRMRDData(filename, d, hdr, acqTrack);

        Col<CxT1> senRB = getISMRMRDSenseMap<CxT1>(d);
        Col<T1> fmRB = getISMRMRDFieldMap<T1>(d);
        Col<CxT1> dataRB;
        Col<T1> kxRB, kyRB, kzRB, tvecRB;
        getCompleteISMRMRDAcqData<T1>(d, acqTrack, 0, 0, 0, 0, 0,
                                       dataRB, kxRB, kyRB, kzRB, tvecRB);

        // Setup
        auto t0 = std::chrono::high_resolution_clock::now();
        Gnufft<T1> Grecon(kxRB.n_rows, (T1)2.0, Nx, Ny, Nz, kxRB, kyRB, kzRB, ix, iy, iz);
        TimeSegmentation<T1, Gnufft<T1>> Arecon(Grecon, fmRB, tvecRB,
                                                 kxRB.n_rows, Ni, 1, 1, 1);
        SENSE<T1, TimeSegmentation<T1, Gnufft<T1>>> Sgrecon(Arecon, senRB,
                                                             kxRB.n_rows, Ni, Nc);
        QuadPenalty<T1> R(Nx, Ny, Nz, (T1)1e-3, 2);
        auto t1 = std::chrono::high_resolution_clock::now();

        // Reconstruct
        Col<CxT1> xhat = reconSolve<T1>(dataRB, Sgrecon, R,
                                          kxRB, kyRB, kzRB, Nx, Ny, Nz, tvecRB, cfg.niter);
        auto t2 = std::chrono::high_resolution_clock::now();

        double setup_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();
        double recon_ms = std::chrono::duration<double, std::milli>(t2 - t1).count();
        T1 err = nrmse<T1>(xhat, phantom);

        char sizeStr[32];
        snprintf(sizeStr, sizeof(sizeStr), "%lux%lu", (unsigned long)Nx, (unsigned long)Ny);
        printf("%-10s %-6lu %-7lu %-5lu %-6lu %-12.1f %-12.1f %-8.4f\n",
               sizeStr, (unsigned long)Nc, (unsigned long)nSpokes,
               (unsigned long)nReadout, (unsigned long)cfg.niter,
               setup_ms, recon_ms, (double)err);

        REQUIRE(!xhat.has_nan());
        REQUIRE(err < (T1)0.60);

        closeISMRMRDData(d, hdr, acqTrack);
        std::remove(filename.c_str());
    }
}
