/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [GenSyntheticISMRMRD.cpp]

    Synopsis    [Generates synthetic ISMRMRD files for testing and benchmarking
                 PowerGrid reconstruction.]

    Description [Creates a Shepp-Logan phantom with synthetic coil sensitivities,
                 generates k-space data via SENSE or pcSENSE forward model,
                 and writes to ISMRMRD HDF5 format compatible with
                 PowerGridIsmrmrd, PowerGridPcSense, and PowerGridPcSenseTimeSeg.]

    Date        [2024]

 *****************************************************************************/

#include "../PowerGrid/Core/PowerGrid.h"
#include "../PowerGrid/IO/processIsmrmrd.hpp"
#include <boost/program_options.hpp>
#include <sstream>
#include <chrono>

// Include synthetic data generators
#include "../PowerGrid/Tests/SyntheticPhantom.hpp"
#include "../PowerGrid/Tests/SyntheticCoils.hpp"
#include "../PowerGrid/Tests/SyntheticTrajectory.hpp"
#include "../PowerGrid/Tests/SyntheticFieldMap.hpp"

namespace po = boost::program_options;
using namespace arma;

// ---------------------------------------------------------------------------
// Write a complete ISMRMRD HDF5 file from synthetic data
// ---------------------------------------------------------------------------
static void writeSyntheticISMRMRD(
    const std::string& filename,
    uword Nx, uword Ny, uword Nz, uword Nc,
    uword nSpokes, uword nReadout,
    const Col<float>& kx, const Col<float>& ky, const Col<float>& kz,
    const Col<float>& tvec,
    const Col<std::complex<float>>& kspaceData,
    const Col<std::complex<float>>& SENSEmap,
    const Col<float>& fieldMap,
    const Col<float>* phaseMaps)
{
    uword Nd = nSpokes * nReadout;
    uword Ni = Nx * Ny * Nz;

    // Remove existing file to avoid appending (ISMRMRD doesn't truncate)
    std::remove(filename.c_str());

    ISMRMRD::Dataset dataset(filename.c_str(), "dataset", true);

    // Build header
    ISMRMRD::IsmrmrdHeader hdr;
    hdr.version = ISMRMRD_XMLHDR_VERSION;
    hdr.experimentalConditions.H1resonanceFrequency_Hz = 128000000;

    ISMRMRD::AcquisitionSystemInformation sysInfo;
    sysInfo.receiverChannels = (unsigned short)Nc;
    hdr.acquisitionSystemInformation = sysInfo;

    ISMRMRD::Encoding enc;
    enc.encodedSpace.matrixSize = ISMRMRD::MatrixSize(
        (unsigned short)Nx, (unsigned short)Ny, (unsigned short)Nz);
    enc.encodedSpace.fieldOfView_mm.x = 256.0f;
    enc.encodedSpace.fieldOfView_mm.y = 256.0f;
    enc.encodedSpace.fieldOfView_mm.z = 5.0f;
    enc.reconSpace.matrixSize = ISMRMRD::MatrixSize(
        (unsigned short)Nx, (unsigned short)Ny, (unsigned short)Nz);
    enc.reconSpace.fieldOfView_mm.x = 256.0f;
    enc.reconSpace.fieldOfView_mm.y = 256.0f;
    enc.reconSpace.fieldOfView_mm.z = 5.0f;
    // ISMRMRD >= 1.14: enum class TrajectoryType; <= 1.3: std::string
#if ISMRMRD_VERSION_MINOR >= 14
    enc.trajectory = ISMRMRD::TrajectoryType::OTHER;
#else
    enc.trajectory = "other";
#endif

    enc.encodingLimits.kspace_encoding_step_1 =
        ISMRMRD::Limit(0, (unsigned short)(nSpokes - 1), 0);
    enc.encodingLimits.kspace_encoding_step_2 = ISMRMRD::Limit(0, 0, 0);
    enc.encodingLimits.slice = ISMRMRD::Limit(0, 0, 0);
    enc.encodingLimits.repetition = ISMRMRD::Limit(0, 0, 0);
    enc.encodingLimits.average = ISMRMRD::Limit(0, 0, 0);
    enc.encodingLimits.contrast = ISMRMRD::Limit(0, 0, 0);
    enc.encodingLimits.phase = ISMRMRD::Limit(0, 0, 0);
    enc.encodingLimits.segment = ISMRMRD::Limit(0, 0, 0);
    enc.encodingLimits.set = ISMRMRD::Limit(0, 0, 0);

    hdr.encoding.push_back(enc);

    std::ostringstream oss;
    ISMRMRD::serialize(hdr, oss);
    dataset.writeHeader(oss.str());

    // Write acquisitions
    for (uword s = 0; s < nSpokes; s++) {
        ISMRMRD::Acquisition acq((uint16_t)nReadout, (uint16_t)Nc, (uint16_t)4);
        acq.idx().kspace_encode_step_1 = (uint16_t)s;
        acq.idx().kspace_encode_step_2 = 0;
        acq.idx().slice = 0;
        acq.idx().repetition = 0;
        acq.idx().average = 0;
        acq.idx().contrast = 0;
        acq.idx().phase = 0;
        acq.idx().segment = 0;
        acq.idx().set = 0;

        for (uword r = 0; r < nReadout; r++) {
            uword idx = s * nReadout + r;
            acq.traj(0, (uint16_t)r) = kx(idx);
            acq.traj(1, (uint16_t)r) = ky(idx);
            acq.traj(2, (uint16_t)r) = kz(idx);
            acq.traj(3, (uint16_t)r) = tvec(idx);
        }

        for (uword c = 0; c < Nc; c++) {
            for (uword r = 0; r < nReadout; r++) {
                uword dataIdx = s * nReadout + r + c * Nd;
                std::complex<float> val = kspaceData(dataIdx);
                acq.data((uint16_t)r, (uint16_t)c) = {val.real(), val.imag()};
            }
        }

        dataset.appendAcquisition(acq);
    }

    // Write NDArrays
    {
        std::vector<size_t> dims = {(size_t)(Ni * Nc)};
        ISMRMRD::NDArray<std::complex<double>> senArray(dims);
        for (uword i = 0; i < Ni * Nc; i++) {
            senArray.getDataPtr()[i] = std::complex<double>(
                (double)SENSEmap(i).real(), (double)SENSEmap(i).imag());
        }
        dataset.appendNDArray("SENSEMap", senArray);
    }

    {
        std::vector<size_t> dims = {(size_t)Ni};
        ISMRMRD::NDArray<double> fmArray(dims);
        for (uword i = 0; i < Ni; i++) {
            fmArray.getDataPtr()[i] = (double)fieldMap(i);
        }
        dataset.appendNDArray("FieldMap", fmArray);
    }

    if (phaseMaps != nullptr) {
        uword pmSize = phaseMaps->n_elem;
        std::vector<size_t> dims = {(size_t)pmSize};
        ISMRMRD::NDArray<double> pmArray(dims);
        for (uword i = 0; i < pmSize; i++) {
            pmArray.getDataPtr()[i] = (double)(*phaseMaps)(i);
        }
        dataset.appendNDArray("PhaseMaps", pmArray);
    }
}

// ---------------------------------------------------------------------------
// Main
// ---------------------------------------------------------------------------
int main(int argc, char** argv) {
    std::string outputFile;
    uword Nx = 64, Ny = 64, Nz = 1;
    uword Nc = 4;
    uword nSpokes = 128, nReadout = 128;
    uword Ns = 1;
    double fieldmapHz = 0.0;
    double readoutTime = 0.01;
    double snrDB = 0.0;

    po::options_description desc("GenSyntheticISMRMRD - Generate synthetic ISMRMRD test data");
    desc.add_options()
        ("help,h", "Show help message")
        ("output,o", po::value<std::string>(&outputFile)->required(), "Output ISMRMRD file path")
        ("Nx,x", po::value<uword>(&Nx), "Image width (default: 64)")
        ("Ny,y", po::value<uword>(&Ny), "Image height (default: 64)")
        ("Nz,z", po::value<uword>(&Nz), "Image depth (default: 1)")
        ("coils,c", po::value<uword>(&Nc), "Number of coils (default: 4)")
        ("spokes,k", po::value<uword>(&nSpokes), "Number of radial spokes (default: 128)")
        ("readout,r", po::value<uword>(&nReadout), "Samples per readout (default: 128)")
        ("fieldmap-hz,f", po::value<double>(&fieldmapHz), "Max off-resonance in Hz (default: 0)")
        ("shots,s", po::value<uword>(&Ns), "Shots for pcSENSE mode (default: 1 = SENSE)")
        ("readout-time,T", po::value<double>(&readoutTime), "Readout duration in seconds (default: 0.01)")
        ("snr", po::value<double>(&snrDB), "Add noise at this SNR in dB (default: 0 = no noise)");

    po::variables_map vm;
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        if (vm.count("help")) {
            std::cout << desc << std::endl;
            return 0;
        }
        po::notify(vm);
    } catch (po::error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cout << desc << std::endl;
        return 1;
    }

    uword Ni = Nx * Ny * Nz;
    uword nSpokesTotal = nSpokes;
    if (Ns > 1) {
        nSpokesTotal = nSpokes; // nSpokes is total across all shots
    }
    uword Nd = nSpokesTotal * nReadout;

    std::cout << "Generating synthetic ISMRMRD data..." << std::endl;
    std::cout << "  Matrix: " << Nx << "x" << Ny << "x" << Nz << std::endl;
    std::cout << "  Coils:  " << Nc << std::endl;
    std::cout << "  Spokes: " << nSpokesTotal << " x " << nReadout << " readout" << std::endl;
    std::cout << "  Shots:  " << Ns << (Ns > 1 ? " (pcSENSE mode)" : " (SENSE mode)") << std::endl;
    if (fieldmapHz > 0.0) {
        std::cout << "  Field map: +/-" << fieldmapHz << " Hz" << std::endl;
    }

    auto t0 = std::chrono::high_resolution_clock::now();

    // Generate phantom
    Col<std::complex<float>> phantom = sheppLogan2D<float>(Nx, Ny);

    // Generate coil sensitivities
    Col<std::complex<float>> SENSEmap = syntheticCoils2D<float>(Nx, Ny, Nc);

    // Generate trajectory
    Col<float> kx, ky, kz;
    radialTrajectory2D<float>(nSpokesTotal, nReadout, Nx, kx, ky, kz);

    // Field map
    Col<float> fieldMap;
    if (fieldmapHz > 0.0) {
        fieldMap = syntheticFieldMap2D<float>(Nx, Ny, (float)fieldmapHz);
    } else {
        fieldMap.zeros(Ni);
    }

    // Timing vector
    Col<float> tvec = linearTimingVector<float>(Nd, (float)readoutTime);

    // Image coordinates
    Col<float> ix, iy, iz;
    initImageSpaceCoords<float>(ix, iy, iz, Nx, Ny, Nz);

    // Generate k-space data
    Col<std::complex<float>> kspaceData;
    Col<float> phaseMaps;
    Col<float>* pmPtr = nullptr;

    if (Ns > 1) {
        // pcSENSE mode: generate phase maps and use pcSENSE forward model
        phaseMaps.set_size(Ni * Ns);
        for (uword s = 0; s < Ns; s++) {
            for (uword i = 0; i < Ni; i++) {
                phaseMaps(i + s * Ni) = (float)s * 0.1f * ((float)(i % Nx) / (float)Nx);
            }
        }
        pmPtr = &phaseMaps;

        pcSENSE<float> Sfwd(kx, ky, kz, Nx, Ny, Nz, Nc, tvec,
                             SENSEmap, fieldMap, phaseMaps);
        kspaceData = Sfwd * phantom;
    } else {
        // SENSE mode
        Gnufft<float> G(Nd, 2.0f, Nx, Ny, Nz, kx, ky, kz, ix, iy, iz);
        SENSE<float, Gnufft<float>> S(G, SENSEmap, Nd, Ni, Nc);
        kspaceData = S * phantom;
    }

    // Optionally add noise
    if (snrDB > 0.0) {
        float signalPower = (float)norm(kspaceData, 2) / std::sqrt((float)kspaceData.n_elem);
        float noiseSigma = signalPower / std::pow(10.0f, (float)snrDB / 20.0f);
        Col<std::complex<float>> noise = noiseSigma * randn<Col<std::complex<float>>>(kspaceData.n_elem);
        kspaceData += noise;
        std::cout << "  Added noise: SNR=" << snrDB << " dB, sigma=" << noiseSigma << std::endl;
    }

    // Write ISMRMRD file
    writeSyntheticISMRMRD(outputFile, Nx, Ny, Nz, Nc, nSpokesTotal, nReadout,
                          kx, ky, kz, tvec, kspaceData, SENSEmap, fieldMap, pmPtr);

    auto t1 = std::chrono::high_resolution_clock::now();
    double elapsed_ms = std::chrono::duration<double, std::milli>(t1 - t0).count();

    std::cout << "  Output: " << outputFile << std::endl;
    std::cout << "  k-space data points: " << Nd << " x " << Nc << " = " << Nd * Nc << std::endl;
    std::cout << "  Generation time: " << elapsed_ms << " ms" << std::endl;

    return 0;
}
