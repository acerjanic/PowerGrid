/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [Gfft.h]

    Synopsis    [Object that represents a uniform discrete Fourier
                    transform implemented via a fast Fourier transform (FFT).]

    Description [Forward transforms are denoted by G*data and adjoint transforms
                    are denoted by G/data. See documentation for more
                    information]

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

/// @file Gfft.h
/// @brief Uniform Cartesian FFT encoding operator.

#ifndef __PowerGrid__Gfft__h
#define __PowerGrid__Gfft__h

#ifdef OPENACC_GPU // GPU Version
#include "cufft.h"
#include "FFT/fftGPU.h"
#include "Gridding/gridding.h"
#include "openacc.h"
#elif OPENACC_MP
#include "FFT/fftGPU.h"
#include "Gridding/gridding.h"
#include "openacc.h"
#elif defined(METAL_COMPUTE) // Apple Metal path: vDSP FFT
#include "Metal/fftAccelerate.h"
#include "FFT/fftCPU.h"    // fallback for non-pow2 sizes
#include "Gridding/gridding.h"
#else // CPU Version
#include "FFT/fftCPU.h"
#include "Gridding/gridding.h"
#endif

/// @brief Uniform Cartesian FFT encoding operator.
///
/// Implements a fully sampled Cartesian Fourier encoding operator using
/// FFTW on CPU or cuFFT on GPU. This is the most efficient encoding operator
/// for Cartesian trajectories.
///
/// Use `G * x` for the forward FFT and `G / d` for the adjoint (inverse FFT).
///
/// @tparam T1  Floating-point precision type (`float` or `double`).
template <typename T1> class Gfft {
  typedef complex<T1> CxT1;

public:
  /// @brief Default constructor. Produces an uninitialized operator.
  Gfft();

  /// @brief Construct a Gfft operator for an Nx x Ny x Nz image grid.
  ///
  /// @param ix  Image x-dimension in pixels (Nx).
  /// @param iy  Image y-dimension in pixels (Ny).
  /// @param iz  Image z-dimension in pixels (Nz); use 1 for 2-D.
  Gfft(uword ix, uword iy, uword iz);

  /// @brief Image x-dimension in pixels.
  uword Nx = 0;
  /// @brief Image y-dimension in pixels.
  uword Ny = 0;
  /// @brief Image z-dimension in pixels.
  uword Nz = 0;
  #ifdef OPENACC_GPU
    /// @brief GPU stream handle (OpenACC/CUDA).
    void *stream;
    /// @brief cuFFT plan handle.
    cufftHandle *plan;
  #endif

  /// @brief Forward FFT: image -> k-space.
  ///
  /// @param d  Input image vector of length Nx*Ny*Nz.
  /// @returns  Output k-space vector of length Nx*Ny*Nz.
  Col<CxT1> operator*(const Col<CxT1> &d) const;

  /// @brief Adjoint FFT (inverse FFT): k-space -> image.
  ///
  /// @param d  Input k-space vector of length Nx*Ny*Nz.
  /// @returns  Output image vector of length Nx*Ny*Nz.
  Col<CxT1> operator/(const Col<CxT1> &d) const;
};

// Explicit Instantiations
extern template class Gfft<float>;
extern template class Gfft<double>;

#endif // __PowerGrid__Gfft__h
