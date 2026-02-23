// Standalone profiling binary for large 3D pcSenseTimeSeg reconstruction.
// Run under Instruments.app or xctrace to capture os_signpost intervals.
//
// Usage:
//   xcrun xctrace record --template 'Time Profiler' \
//     --output profiles/recon3d_large.trace \
//     --launch -- ./build_native/metal_profile

#include "Core/PGIncludes.h"
#include "Operators/pcSenseTimeSeg.h"
#include "Solvers/QuadPenalty.h"
#include "Solvers/solve_pwls_pcg.hpp"

#include "Tests/Helpers/SyntheticPhantom3D.hpp"
#include "Tests/Helpers/SyntheticCoils3D.hpp"
#include "Tests/Helpers/SyntheticFieldMap3D.hpp"
#include "Tests/Helpers/SyntheticTrajectory3D.hpp"

#include <chrono>
#include <iostream>

#ifdef METAL_COMPUTE
#include "Metal/MetalVectorOps.h"
#include "Metal/MetalGridding.h"
#include "Metal/MetalNufftPipeline.h"
#endif

using namespace arma;

int main() {
    using T1 = float;
    using CxT1 = std::complex<T1>;

    const uword Nx = 120, Ny = 120, Nz = 60;
    const uword Ni = Nx * Ny * Nz;
    const uword Nc = 16;
    const uword nInterleaves = 16, nSamplesPerArm = 512;
    const uword R_xy = 2, R_z = 2;
    const uword L = 4;
    const uword niter = 5;

    std::cout << "=== ProfileRecon3D: 120x120x60 pcSenseTimeSeg ===" << std::endl;
    std::cout << "  Grid: " << Nx << "x" << Ny << "x" << Nz
              << " (" << Ni << " voxels)" << std::endl;

    // Generate trajectory
    Col<T1> kx, ky, kz;
    uword Nd_per_shot, Ns;
    stackOfSpiralsTrajectory3D<T1>(nInterleaves, nSamplesPerArm, Nx, Nz,
                                    R_xy, R_z, kx, ky, kz, Nd_per_shot, Ns);

    uword nKzAcq = Nz / R_z;
    Col<T1> tvec = stackOfSpiralsTimingVector<T1>(nSamplesPerArm, nKzAcq, Ns, (T1)0.005);

    std::cout << "  Ns=" << Ns << ", Nd_per_shot=" << Nd_per_shot
              << ", Nc=" << Nc << ", L=" << L << std::endl;

    // Generate synthetic data
    Col<CxT1> phantom = brainPhantom3D<T1>(Nx, Ny, Nz);
    Col<CxT1> SENSEmap = syntheticCoils3D<T1>(Nx, Ny, Nz, Nc);
    Col<T1> fmapHz = syntheticFieldMap3D<T1>(Nx, Ny, Nz, (T1)50.0);
    Col<T1> fmapRad = fmapHz * (T1)(2.0 * M_PI);
    Col<T1> shotPhase = randomShotPhase3D<T1>(Nx, Ny, Nz, Ns, (T1)(M_PI / 2.0), 42);

    std::cout << "  Data generated. Constructing pcSenseTimeSeg operator..." << std::endl;

    // Operator setup
    auto t0 = std::chrono::high_resolution_clock::now();
    pcSenseTimeSeg<T1> P(kx, ky, kz, Nx, Ny, Nz, Nc, tvec, L, 1,
                          SENSEmap, fmapRad, shotPhase);
    auto t1 = std::chrono::high_resolution_clock::now();
    double setupMs = std::chrono::duration<double, std::milli>(t1 - t0).count();
    std::cout << "  Setup: " << setupMs << " ms" << std::endl;

    // Forward simulation
    std::cout << "  Forward simulation..." << std::endl;
    auto t2 = std::chrono::high_resolution_clock::now();
    Col<CxT1> y = P * phantom;
    auto t3 = std::chrono::high_resolution_clock::now();
    double fwdMs = std::chrono::duration<double, std::milli>(t3 - t2).count();
    std::cout << "  Forward: " << fwdMs << " ms, data size: " << y.n_elem << std::endl;

    // Reconstruct
    Col<T1> W = ones<Col<T1>>(y.n_elem);
    QuadPenalty<T1> R(Nx, Ny, Nz, (T1)1e-4, 3);
    Col<CxT1> x0 = zeros<Col<CxT1>>(Ni);

    std::cout << "  Reconstructing (" << niter << " PCG iterations)..." << std::endl;

#ifdef METAL_COMPUTE
    metal_vecops_reset_stats();
    metal_gridding_reset_stats();
    metal_nufft_reset_stats();
#endif

    auto t4 = std::chrono::high_resolution_clock::now();
    Col<CxT1> xhat = solve_pwls_pcg<T1>(x0, P, W, y, R, niter);
    auto t5 = std::chrono::high_resolution_clock::now();
    double reconMs = std::chrono::duration<double, std::milli>(t5 - t4).count();

    T1 nrmse = norm(phantom - xhat) / norm(phantom);
    std::cout << "  Recon: " << reconMs << " ms" << std::endl;
    std::cout << "  NRMSE: " << nrmse << std::endl;
    std::cout << "  Total: " << (setupMs + fwdMs + reconMs) << " ms" << std::endl;

#ifdef METAL_COMPUTE
    std::cout << "\n=== Metal Dispatch Statistics (PCG only) ===" << std::endl;
    std::cout << "  Vector ops:  " << metal_vecops_dispatch_count()
              << " dispatches, " << metal_vecops_wait_seconds() << " s GPU wait" << std::endl;
    std::cout << "  Gridding:    " << metal_gridding_dispatch_count()
              << " dispatches, " << metal_gridding_wait_seconds() << " s GPU wait" << std::endl;
    std::cout << "  NUFFT pipe:  " << metal_nufft_dispatch_count()
              << " dispatches, " << metal_nufft_wait_seconds() << " s GPU wait" << std::endl;
    uint64_t totalDispatches = metal_vecops_dispatch_count()
                             + metal_gridding_dispatch_count()
                             + metal_nufft_dispatch_count();
    double totalWait = metal_vecops_wait_seconds()
                     + metal_gridding_wait_seconds()
                     + metal_nufft_wait_seconds();
    std::cout << "  TOTAL:       " << totalDispatches
              << " dispatches, " << totalWait << " s GPU wait" << std::endl;
    std::cout << "  Wall-clock recon: " << reconMs / 1000.0 << " s" << std::endl;
    std::cout << "  CPU active (approx): " << (reconMs / 1000.0 - totalWait) << " s" << std::endl;
    std::cout << "  Avg wait/dispatch: " << (totalWait / totalDispatches * 1e6) << " us" << std::endl;
#endif

    return 0;
}
