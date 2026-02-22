/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file SyntheticFieldMap3D.hpp
/// @brief 3D B0 inhomogeneity field map generator.
///
/// Models through-slice susceptibility gradients with a quadratic+linear
/// z-dependent pattern typical of frontal sinus / ear canal regions.

#ifndef POWERGRID_TESTS_SYNTHETICFIELDMAP3D_HPP
#define POWERGRID_TESTS_SYNTHETICFIELDMAP3D_HPP

#include <armadillo>

/// Generate a 3D B0 field map with through-slice variation.
///
/// Model: f(x,y,z) = maxHz * (0.3*x^2 + 0.2*y^2 + 0.5*z^2 + 0.2*z)
/// The heavy z-weighting and linear z term create through-slice dephasing
/// that exercises TimeSegmentation correction.
///
/// @param Nx                 Image width
/// @param Ny                 Image height
/// @param Nz                 Number of slices
/// @param maxOffResonanceHz  Peak off-resonance in Hz
/// @return                   Vectorized field map in Hz (Nx*Ny*Nz elements)
template<typename T1>
arma::Col<T1> syntheticFieldMap3D(arma::uword Nx, arma::uword Ny, arma::uword Nz,
                                  T1 maxOffResonanceHz)
{
    arma::uword Ni = Nx * Ny * Nz;
    arma::Col<T1> fmap(Ni);

    // Loop order matches pcSENSE.cpp Cube(Nx,Ny,Nz) vectorization
    for (arma::uword ii = 0; ii < Ny; ii++) {
        T1 y = (T1)2.0 * ((T1)ii - (T1)Ny / (T1)2.0) / (T1)Ny;
        for (arma::uword jj = 0; jj < Nx; jj++) {
            T1 x = (T1)2.0 * ((T1)jj - (T1)Nx / (T1)2.0) / (T1)Nx;
            for (arma::uword kk = 0; kk < Nz; kk++) {
                T1 z = (T1)2.0 * ((T1)kk - (T1)Nz / (T1)2.0) / (T1)Nz;
                arma::uword idx = ii + jj * Nx + kk * Nx * Ny;

                fmap(idx) = maxOffResonanceHz *
                    ((T1)0.3 * x*x + (T1)0.2 * y*y + (T1)0.5 * z*z + (T1)0.2 * z);
            }
        }
    }

    return fmap;
}

#endif // POWERGRID_TESTS_SYNTHETICFIELDMAP3D_HPP
