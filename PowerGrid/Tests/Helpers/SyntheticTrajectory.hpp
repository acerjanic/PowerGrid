/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file SyntheticTrajectory.hpp
/// @brief Radial and spiral k-space trajectory generators for integration tests.

#ifndef POWERGRID_TESTS_SYNTHETICTRAJECTORY_HPP
#define POWERGRID_TESTS_SYNTHETICTRAJECTORY_HPP

#include <armadillo>
#include <cmath>

/// Generate a 2D radial (golden-angle or uniform) k-space trajectory.
///
/// Produces nSpokes spokes with nReadout samples each, spanning
/// [-Nx/2, Nx/2] in k-space. Uniform angular spacing: angle = pi * spoke / nSpokes.
///
/// @param nSpokes    Number of radial spokes
/// @param nReadout   Number of samples per spoke
/// @param Nx         Image dimension (determines kmax = Nx/2)
/// @param[out] kx    k-space x coordinates (nSpokes * nReadout elements)
/// @param[out] ky    k-space y coordinates
/// @param[out] kz    k-space z coordinates (all zeros for 2D)
template<typename T1>
void radialTrajectory2D(arma::uword nSpokes, arma::uword nReadout,
                        arma::uword Nx,
                        arma::Col<T1>& kx, arma::Col<T1>& ky, arma::Col<T1>& kz)
{
    arma::uword totalSamples = nSpokes * nReadout;
    kx.set_size(totalSamples);
    ky.set_size(totalSamples);
    kz.zeros(totalSamples);

    T1 kmax = (T1)Nx / (T1)2.0;

    for (arma::uword s = 0; s < nSpokes; s++) {
        T1 angle = M_PI * (T1)s / (T1)nSpokes;
        T1 cosA = std::cos(angle);
        T1 sinA = std::sin(angle);

        for (arma::uword r = 0; r < nReadout; r++) {
            // Linear readout from -kmax to +kmax
            T1 kr = -kmax + (T1)2.0 * kmax * (T1)r / (T1)(nReadout - 1);
            arma::uword idx = s * nReadout + r;
            kx(idx) = kr * cosA;
            ky(idx) = kr * sinA;
        }
    }
}

/// Generate a 2D Archimedean spiral k-space trajectory.
///
/// @param nInterleaves     Number of spiral interleaves
/// @param nSamplesPerArm   Samples per spiral arm
/// @param Nx               Image dimension (determines kmax = Nx/2)
/// @param[out] kx          k-space x coordinates
/// @param[out] ky          k-space y coordinates
/// @param[out] kz          k-space z coordinates (all zeros)
template<typename T1>
void spiralTrajectory2D(arma::uword nInterleaves, arma::uword nSamplesPerArm,
                        arma::uword Nx,
                        arma::Col<T1>& kx, arma::Col<T1>& ky, arma::Col<T1>& kz)
{
    arma::uword totalSamples = nInterleaves * nSamplesPerArm;
    kx.set_size(totalSamples);
    ky.set_size(totalSamples);
    kz.zeros(totalSamples);

    T1 kmax = (T1)Nx / (T1)2.0;
    // Number of turns determines how tightly the spiral winds
    T1 nTurns = (T1)Nx / (T1)(2.0 * nInterleaves);

    for (arma::uword il = 0; il < nInterleaves; il++) {
        T1 rotAngle = (T1)2.0 * M_PI * (T1)il / (T1)nInterleaves;

        for (arma::uword s = 0; s < nSamplesPerArm; s++) {
            T1 t = (T1)s / (T1)(nSamplesPerArm - 1);
            T1 r = kmax * t;
            T1 theta = nTurns * (T1)2.0 * M_PI * t + rotAngle;

            arma::uword idx = il * nSamplesPerArm + s;
            kx(idx) = r * std::cos(theta);
            ky(idx) = r * std::sin(theta);
        }
    }
}

/// Generate image-space coordinates matching PowerGrid convention.
///
/// For a 2D image of size Nx x Ny, returns ix, iy, iz vectors where:
///   ix(j,i,k) = (i - Nx/2) / Nx
///   iy(j,i,k) = (j - Ny/2) / Ny
///   iz = 0
///
/// @param Nx       Image width
/// @param Ny       Image height
/// @param[out] ix  x image coordinates (Nx*Ny elements)
/// @param[out] iy  y image coordinates
/// @param[out] iz  z image coordinates (all zeros)
template<typename T1>
void imageCoordinates2D(arma::uword Nx, arma::uword Ny,
                        arma::Col<T1>& ix, arma::Col<T1>& iy, arma::Col<T1>& iz)
{
    arma::uword Ni = Nx * Ny;
    ix.set_size(Ni);
    iy.set_size(Ni);
    iz.zeros(Ni);

    // PowerGrid convention: outer loop over y, inner loop over x,
    // but stored column-major (matches pcSENSE constructor ordering)
    for (arma::uword jj = 0; jj < Ny; jj++) {       // y
        for (arma::uword ii = 0; ii < Nx; ii++) {    // x
            arma::uword idx = jj + ii * Ny;  // column-major in Cube(Ny, Nx, 1)
            ix(idx) = ((T1)ii - (T1)Nx / (T1)2.0) / (T1)Nx;
            iy(idx) = ((T1)jj - (T1)Ny / (T1)2.0) / (T1)Ny;
        }
    }
}

/// Generate a linear timing vector for a single readout.
///
/// @param nSamples     Number of samples
/// @param readoutTime  Total readout duration in seconds (e.g. 0.01 for 10ms)
/// @return             Timing vector from 0 to readoutTime
template<typename T1>
arma::Col<T1> linearTimingVector(arma::uword nSamples, T1 readoutTime) {
    return arma::linspace<arma::Col<T1>>((T1)0, readoutTime, nSamples);
}

#endif // POWERGRID_TESTS_SYNTHETICTRAJECTORY_HPP
