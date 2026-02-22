/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file SyntheticPhantom.hpp
/// @brief Analytical Shepp-Logan phantom generator for integration tests.
///
/// Implements the modified Shepp-Logan phantom (Toft 1996) with 10 ellipses
/// and improved contrast over the original 1974 version.

#ifndef POWERGRID_TESTS_SYNTHETICPHANTOM_HPP
#define POWERGRID_TESTS_SYNTHETICPHANTOM_HPP

#include <armadillo>
#include <cmath>

/// Shepp-Logan ellipse parameters: {intensity, a, b, x0, y0, phi (degrees)}
/// Using the modified Shepp-Logan values (Toft 1996) for better visual contrast.
template<typename T1>
struct EllipseParam {
    T1 intensity;
    T1 a;       // semi-axis x
    T1 b;       // semi-axis y
    T1 x0;      // center x
    T1 y0;      // center y
    T1 phi;     // rotation angle in degrees
};

/// Generate a 2D modified Shepp-Logan phantom image.
///
/// @param Nx  Number of pixels in x direction
/// @param Ny  Number of pixels in y direction
/// @return    Nx*Ny x 1 complex column vector (vectorized image, column-major)
template<typename T1>
arma::Col<std::complex<T1>> sheppLogan2D(arma::uword Nx, arma::uword Ny) {
    // Modified Shepp-Logan parameters (Toft 1996)
    // {intensity, a, b, x0, y0, phi_deg}
    const int nEllipses = 10;
    EllipseParam<T1> ellipses[nEllipses] = {
        { 1.0,   0.69,   0.92,   0.0,    0.0,    0.0   },  // outer skull
        {-0.8,   0.6624, 0.8740, 0.0,   -0.0184, 0.0   },  // inner skull
        {-0.2,   0.1100, 0.3100, 0.22,   0.0,   -18.0  },  // left eye
        {-0.2,   0.1600, 0.4100,-0.22,   0.0,    18.0  },  // right eye
        { 0.1,   0.2100, 0.2500, 0.0,    0.35,   0.0   },  // nose
        { 0.1,   0.0460, 0.0460, 0.0,    0.1,    0.0   },  // mouth
        { 0.1,   0.0460, 0.0460, 0.0,   -0.1,    0.0   },  // left ear
        { 0.1,   0.0460, 0.0230,-0.08,  -0.605,  0.0   },  // right ear
        { 0.1,   0.0230, 0.0230, 0.0,   -0.606,  0.0   },  // small detail 1
        { 0.1,   0.0230, 0.0460, 0.06,  -0.605,  0.0   },  // small detail 2
    };

    arma::Col<std::complex<T1>> img(Nx * Ny, arma::fill::zeros);

    for (int e = 0; e < nEllipses; e++) {
        T1 A   = ellipses[e].intensity;
        T1 a   = ellipses[e].a;
        T1 b   = ellipses[e].b;
        T1 x0  = ellipses[e].x0;
        T1 y0  = ellipses[e].y0;
        T1 phi = ellipses[e].phi * M_PI / (T1)180.0;

        T1 cosPhi = std::cos(phi);
        T1 sinPhi = std::sin(phi);

        for (arma::uword jj = 0; jj < Ny; jj++) {
            // y coordinate: maps pixel jj to [-1, 1]
            T1 y = (T1)2.0 * ((T1)jj - (T1)Ny / (T1)2.0) / (T1)Ny;
            for (arma::uword ii = 0; ii < Nx; ii++) {
                // x coordinate: maps pixel ii to [-1, 1]
                T1 x = (T1)2.0 * ((T1)ii - (T1)Nx / (T1)2.0) / (T1)Nx;

                // Rotate and translate to ellipse frame
                T1 xr = cosPhi * (x - x0) + sinPhi * (y - y0);
                T1 yr = -sinPhi * (x - x0) + cosPhi * (y - y0);

                // Check if inside ellipse
                T1 val = (xr * xr) / (a * a) + (yr * yr) / (b * b);
                if (val <= (T1)1.0) {
                    // Column-major: index = ii + jj * Nx
                    img(ii + jj * Nx) += std::complex<T1>(A, 0);
                }
            }
        }
    }

    return img;
}

#endif // POWERGRID_TESTS_SYNTHETICPHANTOM_HPP
