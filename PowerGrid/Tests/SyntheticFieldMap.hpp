/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file SyntheticFieldMap.hpp
/// @brief Synthetic B0 field map generator for TimeSegmentation integration tests.
///
/// Generates smooth quadratic off-resonance maps typical of through-plane
/// susceptibility gradients.

#ifndef POWERGRID_TESTS_SYNTHETICFIELDMAP_HPP
#define POWERGRID_TESTS_SYNTHETICFIELDMAP_HPP

#include <armadillo>
#include <cmath>

/// Generate a smooth 2D B0 inhomogeneity field map.
///
/// Models a quadratic off-resonance pattern:
///   f(x,y) = maxHz * (x^2 + 0.5*y^2)
/// where x,y are in [-1,1]. This represents a typical through-plane
/// gradient from susceptibility effects near air-tissue interfaces.
///
/// The returned map is in Hz (not radians/second). Multiply by 2*pi
/// to convert to radians/second for TimeSegmentation.
///
/// @param Nx                 Image width
/// @param Ny                 Image height
/// @param maxOffResonanceHz  Maximum off-resonance in Hz (e.g. 100.0)
/// @return                   Vectorized Nx*Ny field map in Hz
template<typename T1>
arma::Col<T1> syntheticFieldMap2D(arma::uword Nx, arma::uword Ny,
                                  T1 maxOffResonanceHz)
{
    arma::Col<T1> fmap(Nx * Ny);

    for (arma::uword jj = 0; jj < Ny; jj++) {
        T1 y = (T1)2.0 * ((T1)jj - (T1)Ny / (T1)2.0) / (T1)Ny;
        for (arma::uword ii = 0; ii < Nx; ii++) {
            T1 x = (T1)2.0 * ((T1)ii - (T1)Nx / (T1)2.0) / (T1)Nx;
            fmap(ii + jj * Nx) = maxOffResonanceHz * (x * x + (T1)0.5 * y * y);
        }
    }

    return fmap;
}

#endif // POWERGRID_TESTS_SYNTHETICFIELDMAP_HPP
