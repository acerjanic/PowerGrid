/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file SyntheticCoils.hpp
/// @brief Analytical coil sensitivity map generator for integration tests.
///
/// Models Nc receive coils equally spaced around the FOV using a Biot-Savart
/// inspired 1/r sensitivity profile with smooth phase variation.

#ifndef POWERGRID_TESTS_SYNTHETICCOILS_HPP
#define POWERGRID_TESTS_SYNTHETICCOILS_HPP

#include <armadillo>
#include <cmath>
#include <complex>

/// Generate 2D coil sensitivity maps for Nc coils arranged in a ring.
///
/// Each coil has a 1/r magnitude falloff from its position and a smooth
/// linear phase gradient. The maps are normalized so that the sum-of-squares
/// across coils equals 1 at each voxel.
///
/// @param Nx          Image width
/// @param Ny          Image height
/// @param Nc          Number of coils
/// @param coilRadius  Distance of coil centers from FOV center (in FOV units,
///                    e.g. 1.5 means 1.5x the half-FOV)
/// @return            (Nx*Ny*Nc) x 1 complex column vector: Nc maps concatenated
template<typename T1>
arma::Col<std::complex<T1>> syntheticCoils2D(
    arma::uword Nx, arma::uword Ny, arma::uword Nc,
    T1 coilRadius = (T1)1.5)
{
    arma::Mat<std::complex<T1>> SMap(Nx * Ny, Nc, arma::fill::zeros);

    for (arma::uword cc = 0; cc < Nc; cc++) {
        // Coil position on a ring
        T1 angle = (T1)2.0 * M_PI * (T1)cc / (T1)Nc;
        T1 cx = coilRadius * std::cos(angle);
        T1 cy = coilRadius * std::sin(angle);

        for (arma::uword jj = 0; jj < Ny; jj++) {
            T1 y = (T1)2.0 * ((T1)jj - (T1)Ny / (T1)2.0) / (T1)Ny;
            for (arma::uword ii = 0; ii < Nx; ii++) {
                T1 x = (T1)2.0 * ((T1)ii - (T1)Nx / (T1)2.0) / (T1)Nx;

                // Distance from coil to voxel
                T1 dx = x - cx;
                T1 dy = y - cy;
                T1 dist = std::sqrt(dx * dx + dy * dy);

                // Magnitude: 1/r with saturation to avoid singularity
                T1 mag = (T1)1.0 / ((T1)1.0 + dist * dist);

                // Phase: smooth linear gradient from coil position
                T1 phase = std::atan2(dy, dx);

                SMap(ii + jj * Nx, cc) = std::complex<T1>(
                    mag * std::cos(phase),
                    mag * std::sin(phase));
            }
        }
    }

    // Normalize: sum-of-squares across coils = 1 at each voxel
    arma::Col<T1> sos(Nx * Ny, arma::fill::zeros);
    for (arma::uword cc = 0; cc < Nc; cc++) {
        sos += arma::real(SMap.col(cc) % arma::conj(SMap.col(cc)));
    }
    sos = arma::sqrt(sos);
    // Avoid division by zero
    sos.elem(arma::find(sos < (T1)1e-10)).fill((T1)1.0);

    for (arma::uword cc = 0; cc < Nc; cc++) {
        SMap.col(cc) /= arma::conv_to<arma::Col<std::complex<T1>>>::from(sos);
    }

    return arma::vectorise(SMap);
}

#endif // POWERGRID_TESTS_SYNTHETICCOILS_HPP
