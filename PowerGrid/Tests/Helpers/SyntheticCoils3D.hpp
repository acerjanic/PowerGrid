/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file SyntheticCoils3D.hpp
/// @brief 3D coil sensitivity map generator for multi-channel MRI.
///
/// Extends the 2D coil model to 3D: coils on a ring in the xy-plane with
/// 1/(1+r^2) magnitude falloff and atan2 phase variation, SoS-normalized.

#ifndef POWERGRID_TESTS_SYNTHETICCOILS3D_HPP
#define POWERGRID_TESTS_SYNTHETICCOILS3D_HPP

#include <armadillo>
#include <cmath>

/// Generate 3D coil sensitivity maps for Nc coils.
///
/// Coils are placed on a ring of radius coilRadius in the xy-plane at z=0.
/// Sensitivity magnitude: 1/(1+dist^2), phase: atan2(dy, dx).
/// Maps are SoS-normalized so that sum-of-squares = 1 at each voxel.
///
/// @param Nx         Image width
/// @param Ny         Image height
/// @param Nz         Number of slices
/// @param Nc         Number of coils
/// @param coilRadius Distance of coils from FOV center (in normalized coords)
/// @return           Vectorized (Nx*Ny*Nz*Nc) complex column (Nc maps concatenated)
template<typename T1>
arma::Col<std::complex<T1>> syntheticCoils3D(
    arma::uword Nx, arma::uword Ny, arma::uword Nz, arma::uword Nc,
    T1 coilRadius = (T1)1.5)
{
    arma::uword Ni = Nx * Ny * Nz;
    arma::Mat<std::complex<T1>> SMap(Ni, Nc, arma::fill::zeros);

    for (arma::uword cc = 0; cc < Nc; cc++) {
        T1 angle = (T1)2.0 * M_PI * (T1)cc / (T1)Nc;
        T1 cx = coilRadius * std::cos(angle);
        T1 cy = coilRadius * std::sin(angle);
        T1 cz = (T1)0.0;

        // Loop order matches pcSENSE.cpp Cube(Nx,Ny,Nz) vectorization
        for (arma::uword ii = 0; ii < Ny; ii++) {
            T1 y = (T1)2.0 * ((T1)ii - (T1)Ny / (T1)2.0) / (T1)Ny;
            for (arma::uword jj = 0; jj < Nx; jj++) {
                T1 x = (T1)2.0 * ((T1)jj - (T1)Nx / (T1)2.0) / (T1)Nx;
                for (arma::uword kk = 0; kk < Nz; kk++) {
                    T1 z = (T1)2.0 * ((T1)kk - (T1)Nz / (T1)2.0) / (T1)Nz;
                    arma::uword idx = ii + jj * Nx + kk * Nx * Ny;

                    T1 dx = x - cx;
                    T1 dy = y - cy;
                    T1 dz = z - cz;
                    T1 dist2 = dx*dx + dy*dy + dz*dz;
                    T1 mag = (T1)1.0 / ((T1)1.0 + dist2);
                    T1 phs = std::atan2(dy, dx);

                    SMap(idx, cc) = std::polar(mag, phs);
                }
            }
        }
    }

    // Normalize by sum-of-squares across coils at each voxel
    for (arma::uword i = 0; i < Ni; i++) {
        T1 sos = (T1)0.0;
        for (arma::uword cc = 0; cc < Nc; cc++) {
            T1 m = std::abs(SMap(i, cc));
            sos += m * m;
        }
        T1 norm = std::sqrt(sos);
        if (norm > (T1)1e-10) {
            for (arma::uword cc = 0; cc < Nc; cc++) {
                SMap(i, cc) /= norm;
            }
        }
    }

    return arma::vectorise(SMap);
}

#endif // POWERGRID_TESTS_SYNTHETICCOILS3D_HPP
