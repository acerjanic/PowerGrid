/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file SyntheticTrajectory3D.hpp
/// @brief Stack-of-spirals trajectory generator with SENSE acceleration,
///        timing vector, image coordinates, and per-shot phase ramps.

#ifndef POWERGRID_TESTS_SYNTHETICTRAJECTORY3D_HPP
#define POWERGRID_TESTS_SYNTHETICTRAJECTORY3D_HPP

#include <armadillo>
#include <cmath>
#include <random>

/// Generate a 3D stack-of-spirals k-space trajectory.
///
/// Each "shot" consists of one Archimedean spiral interleave replicated across
/// all acquired kz positions. With R_xy acceleration, nInterleaves/R_xy
/// interleaves are acquired. With R_z acceleration, Nz/R_z kz encodes are
/// acquired (uniformly spaced).
///
/// Data layout in output vectors: shot-major, then kz-encode, then sample.
/// pcSENSE/pcSenseTimeSeg reshape kx/ky/kz into Mat(Nd, Ns), so each shot's
/// data forms one contiguous column.
///
/// @param nInterleaves      Total interleaves for fully-sampled in-plane
/// @param nSamplesPerArm    Samples per spiral arm
/// @param Nx                In-plane matrix size (determines kmax = Nx/2)
/// @param Nz                Fully-sampled kz encodes
/// @param R_xy              In-plane acceleration factor
/// @param R_z               Through-slice acceleration factor
/// @param[out] kx           k-space x coordinates (Nd_total = Nd_per_shot * Ns)
/// @param[out] ky           k-space y coordinates
/// @param[out] kz           k-space z coordinates
/// @param[out] Nd_per_shot  Data points per shot
/// @param[out] Ns           Number of shots
template<typename T1>
void stackOfSpiralsTrajectory3D(
    arma::uword nInterleaves, arma::uword nSamplesPerArm,
    arma::uword Nx, arma::uword Nz,
    arma::uword R_xy, arma::uword R_z,
    arma::Col<T1>& kx, arma::Col<T1>& ky, arma::Col<T1>& kz,
    arma::uword& Nd_per_shot, arma::uword& Ns)
{
    Ns = nInterleaves / R_xy;
    arma::uword nKzAcq = Nz / R_z;
    Nd_per_shot = nSamplesPerArm * nKzAcq;
    arma::uword Nd_total = Nd_per_shot * Ns;

    kx.set_size(Nd_total);
    ky.set_size(Nd_total);
    kz.set_size(Nd_total);

    T1 kmax = (T1)Nx / (T1)2.0;
    T1 nTurns = (T1)Nx / ((T1)2.0 * (T1)nInterleaves);

    for (arma::uword ss = 0; ss < Ns; ss++) {
        // Rotation angle: skip by R_xy interleaves
        T1 rotAngle = (T1)2.0 * M_PI * (T1)(ss * R_xy) / (T1)nInterleaves;

        for (arma::uword kzIdx = 0; kzIdx < nKzAcq; kzIdx++) {
            // kz value: uniformly spaced with stride R_z, centered
            T1 kzVal = -(T1)Nz / (T1)2.0 + (T1)kzIdx * (T1)R_z + (T1)R_z / (T1)2.0;

            for (arma::uword s = 0; s < nSamplesPerArm; s++) {
                T1 t = (T1)s / (T1)(nSamplesPerArm - 1);
                T1 r = kmax * t;
                T1 theta = nTurns * (T1)2.0 * M_PI * t + rotAngle;

                arma::uword idx = ss * Nd_per_shot + kzIdx * nSamplesPerArm + s;
                kx(idx) = r * std::cos(theta);
                ky(idx) = r * std::sin(theta);
                kz(idx) = kzVal;
            }
        }
    }
}

/// Generate timing vector for stack-of-spirals readout.
///
/// Each shot acquires kz encodes sequentially, each taking readoutTimePerArm.
/// Timing resets to 0 for each shot (pcSENSE reshapes to Mat(Nd,Ns) columns).
///
/// @param nSamplesPerArm    Samples per spiral arm
/// @param nKzAcq            Acquired kz encodes per shot
/// @param Ns                Number of shots
/// @param readoutTimePerArm Readout duration per spiral arm (seconds)
/// @return                  Timing vector of length Nd_per_shot * Ns
template<typename T1>
arma::Col<T1> stackOfSpiralsTimingVector(
    arma::uword nSamplesPerArm, arma::uword nKzAcq,
    arma::uword Ns, T1 readoutTimePerArm)
{
    arma::uword Nd_per_shot = nSamplesPerArm * nKzAcq;
    arma::Col<T1> tvec(Nd_per_shot * Ns);

    for (arma::uword ss = 0; ss < Ns; ss++) {
        for (arma::uword kzIdx = 0; kzIdx < nKzAcq; kzIdx++) {
            T1 kzOffset = (T1)kzIdx * readoutTimePerArm;
            for (arma::uword s = 0; s < nSamplesPerArm; s++) {
                T1 t = kzOffset + (T1)s * readoutTimePerArm / (T1)(nSamplesPerArm - 1);
                arma::uword idx = ss * Nd_per_shot + kzIdx * nSamplesPerArm + s;
                tvec(idx) = t;
            }
        }
    }

    return tvec;
}

/// Generate 3D image-space coordinates matching pcSENSE Cube(Nx,Ny,Nz) layout.
///
/// @param Nx, Ny, Nz  Image dimensions
/// @param[out] ix, iy, iz  Normalized coordinates in [-0.5, 0.5)
template<typename T1>
void imageCoordinates3D(arma::uword Nx, arma::uword Ny, arma::uword Nz,
                        arma::Col<T1>& ix, arma::Col<T1>& iy, arma::Col<T1>& iz)
{
    arma::uword Ni = Nx * Ny * Nz;
    ix.set_size(Ni);
    iy.set_size(Ni);
    iz.set_size(Ni);

    // Matches pcSENSE.cpp constructor loop exactly
    for (arma::uword ii = 0; ii < Ny; ii++) {
        for (arma::uword jj = 0; jj < Nx; jj++) {
            for (arma::uword kk = 0; kk < Nz; kk++) {
                arma::uword idx = ii + jj * Nx + kk * Nx * Ny;
                ix(idx) = ((T1)jj - (T1)Nx / (T1)2.0) / (T1)Nx;
                iy(idx) = ((T1)ii - (T1)Ny / (T1)2.0) / (T1)Ny;
                iz(idx) = ((T1)kk - (T1)Nz / (T1)2.0) / (T1)Nz;
            }
        }
    }
}

/// Generate random per-shot linear phase ramps (models motion-induced phase).
///
/// Each shot gets a random linear ramp: PMap(x,y,z,ss) = gx*x + gy*y + gz*z
/// with random coefficients in [-maxPhase, maxPhase].
///
/// @param Nx, Ny, Nz  Image dimensions
/// @param Ns          Number of shots
/// @param maxPhase    Maximum phase gradient magnitude (radians)
/// @param seed        RNG seed for reproducibility
/// @return            Vectorized phase map (Ni*Ns elements, radians)
template<typename T1>
arma::Col<T1> randomShotPhase3D(
    arma::uword Nx, arma::uword Ny, arma::uword Nz,
    arma::uword Ns, T1 maxPhase, unsigned int seed = 42)
{
    arma::uword Ni = Nx * Ny * Nz;
    arma::Mat<T1> PMap(Ni, Ns, arma::fill::zeros);

    std::mt19937 rng(seed);
    std::uniform_real_distribution<T1> dist(-maxPhase, maxPhase);

    for (arma::uword ss = 0; ss < Ns; ss++) {
        T1 gx = dist(rng);
        T1 gy = dist(rng);
        T1 gz = dist(rng) * (T1)0.3;  // smaller z phase

        // Matches pcSENSE Cube loop order
        for (arma::uword ii = 0; ii < Ny; ii++) {
            T1 y = (T1)2.0 * ((T1)ii - (T1)Ny / (T1)2.0) / (T1)Ny;
            for (arma::uword jj = 0; jj < Nx; jj++) {
                T1 x = (T1)2.0 * ((T1)jj - (T1)Nx / (T1)2.0) / (T1)Nx;
                for (arma::uword kk = 0; kk < Nz; kk++) {
                    T1 z = (T1)2.0 * ((T1)kk - (T1)Nz / (T1)2.0) / (T1)Nz;
                    arma::uword idx = ii + jj * Nx + kk * Nx * Ny;
                    PMap(idx, ss) = gx * x + gy * y + gz * z;
                }
            }
        }
    }

    return arma::vectorise(PMap);
}

#endif // POWERGRID_TESTS_SYNTHETICTRAJECTORY3D_HPP
