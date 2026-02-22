/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file SyntheticPhantom3D.hpp
/// @brief Multi-compartment 3D brain-like phantom generator.
///
/// Generates a simplified brain phantom with skull, gray matter, white matter,
/// and lateral ventricles using nested 3D ellipsoids with overwrite semantics.

#ifndef POWERGRID_TESTS_SYNTHETICPHANTOM3D_HPP
#define POWERGRID_TESTS_SYNTHETICPHANTOM3D_HPP

#include <armadillo>
#include <cmath>

template<typename T1>
struct Ellipsoid3D {
    T1 intensity;
    T1 ax, ay, az;   // semi-axes
    T1 x0, y0, z0;   // center
};

/// Generate a 3D multi-compartment brain-like phantom.
///
/// Compartments (painted in order, later overwrites earlier):
///   1. Skull shell (0.20)
///   2. Gray matter (0.80)
///   3. White matter core (1.00)
///   4. Left lateral ventricle (0.30)
///   5. Right lateral ventricle (0.30)
///
/// @param Nx  Image width
/// @param Ny  Image height
/// @param Nz  Number of slices
/// @return    Vectorized complex column (Nx*Ny*Nz), matching pcSENSE Cube layout
template<typename T1>
arma::Col<std::complex<T1>> brainPhantom3D(arma::uword Nx, arma::uword Ny, arma::uword Nz) {
    const int nEllipsoids = 5;
    Ellipsoid3D<T1> ellipsoids[nEllipsoids] = {
        // intensity  ax    ay    az    x0     y0    z0
        {  0.20,     0.69, 0.92, 0.90, 0.0,   0.0,  0.0  },  // skull
        {  0.80,     0.62, 0.85, 0.82, 0.0,   0.0,  0.0  },  // gray matter
        {  1.00,     0.40, 0.55, 0.50, 0.0,   0.0,  0.0  },  // white matter
        {  0.30,     0.08, 0.15, 0.12, 0.15,  0.0,  0.0  },  // left ventricle
        {  0.30,     0.08, 0.15, 0.12,-0.15,  0.0,  0.0  },  // right ventricle
    };

    arma::uword Ni = Nx * Ny * Nz;
    arma::Col<std::complex<T1>> img(Ni, arma::fill::zeros);

    // Loop order matches pcSENSE.cpp Cube(Nx, Ny, Nz) layout:
    //   ii indexes rows (0..Ny-1), jj indexes cols (0..Nx-1), kk indexes slices
    //   vectorise gives linear index = ii + jj*Nx + kk*Nx*Ny
    // Note: when Nx==Ny (our use case), this is consistent with the Cube dimensions.
    for (arma::uword ii = 0; ii < Ny; ii++) {
        T1 y = (T1)2.0 * ((T1)ii - (T1)Ny / (T1)2.0) / (T1)Ny;
        for (arma::uword jj = 0; jj < Nx; jj++) {
            T1 x = (T1)2.0 * ((T1)jj - (T1)Nx / (T1)2.0) / (T1)Nx;
            for (arma::uword kk = 0; kk < Nz; kk++) {
                T1 z = (T1)2.0 * ((T1)kk - (T1)Nz / (T1)2.0) / (T1)Nz;
                arma::uword idx = ii + jj * Nx + kk * Nx * Ny;

                // Overwrite semantics: later ellipsoids replace earlier ones
                for (int e = 0; e < nEllipsoids; e++) {
                    T1 dx = (x - ellipsoids[e].x0) / ellipsoids[e].ax;
                    T1 dy = (y - ellipsoids[e].y0) / ellipsoids[e].ay;
                    T1 dz = (z - ellipsoids[e].z0) / ellipsoids[e].az;
                    if (dx*dx + dy*dy + dz*dz <= (T1)1.0) {
                        img(idx) = std::complex<T1>(ellipsoids[e].intensity, 0);
                    }
                }
            }
        }
    }

    return img;
}

#endif // POWERGRID_TESTS_SYNTHETICPHANTOM3D_HPP
