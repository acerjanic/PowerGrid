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
/// Two staggered rings of coils at different z-positions, mimicking a
/// cylindrical receive array. Nc/2 coils per ring with angular offset,
/// 1/(1+r^2) magnitude falloff and atan2 phase variation, SoS-normalized.

#ifndef POWERGRID_TESTS_SYNTHETICCOILS3D_HPP
#define POWERGRID_TESTS_SYNTHETICCOILS3D_HPP

#include <armadillo>
#include <cmath>

/// Generate 3D coil sensitivity maps for Nc coils.
///
/// Coils are arranged in two staggered rings on a cylinder of radius
/// coilRadius.  The lower ring (Nc/2 coils) sits at z = -ringZ and
/// the upper ring (Nc/2 coils) at z = +ringZ, offset by half an
/// angular step to provide better azimuthal coverage.
///
/// Sensitivity magnitude: 1/(1+dist^2), phase: atan2 of displacement
/// projected into each coil's local tangent plane.
///
/// Maps are SoS-normalized so that sum-of-squares = 1 at each voxel.
///
/// @param Nx         Image width
/// @param Ny         Image height
/// @param Nz         Number of slices
/// @param Nc         Number of coils (must be even)
/// @param coilRadius Distance of coils from FOV center in xy (normalized)
/// @param ringZ      z-offset of each ring from center (normalized)
/// @return           Vectorized (Nx*Ny*Nz*Nc) complex column (Nc maps concatenated)
template<typename T1>
arma::Col<std::complex<T1>> syntheticCoils3D(
    arma::uword Nx, arma::uword Ny, arma::uword Nz, arma::uword Nc,
    T1 coilRadius = (T1)1.5, T1 ringZ = (T1)0.5)
{
    arma::uword Ni = Nx * Ny * Nz;
    arma::uword coilsPerRing = Nc / 2;
    arma::Mat<std::complex<T1>> SMap(Ni, Nc, arma::fill::zeros);

    for (arma::uword cc = 0; cc < Nc; cc++) {
        // Determine which ring this coil belongs to
        arma::uword ring = cc / coilsPerRing;   // 0 = lower, 1 = upper
        arma::uword idxInRing = cc % coilsPerRing;

        // Angular position: upper ring offset by half a step
        T1 angularStep = (T1)2.0 * M_PI / (T1)coilsPerRing;
        T1 angle = angularStep * (T1)idxInRing;
        if (ring == 1) {
            angle += angularStep * (T1)0.5;  // stagger upper ring
        }

        T1 cx = coilRadius * std::cos(angle);
        T1 cy = coilRadius * std::sin(angle);
        T1 cz = (ring == 0) ? -ringZ : ringZ;

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
