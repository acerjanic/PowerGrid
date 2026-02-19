/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file fftAccelerate.h
/// @brief vDSP/Accelerate FFT interface mirroring fftCPU.h.
///
/// Provides in-place 2D and 3D complex FFT/IFFT using Apple's vDSP library
/// (part of Accelerate.framework).  Data layout is interleaved complex
/// (pairs of floats).
///
/// **Restriction:** vDSP requires power-of-2 dimensions.  If a dimension is
/// not a power of two the function delegates automatically to the FFTW-backed
/// `fftCPU` routines (which handle arbitrary sizes).  The double-precision
/// overloads always delegate to FFTW because Metal Shading Language does not
/// support `double`.

#pragma once
#ifdef METAL_COMPUTE

#include <Accelerate/Accelerate.h>
#include <armadillo>
using arma::uword;

// ---------------------------------------------------------------------------
// float overloads — vDSP with pow-of-2 check; fall back to FFTW otherwise
// ---------------------------------------------------------------------------

/// @brief In-place 2D forward FFT (float, interleaved).
/// @param d_data  Interleaved complex buffer of size 2*nx*ny floats.
/// @param nx      Size in the first (row) dimension.
/// @param ny      Size in the second (column) dimension.
void fft2dAccelerate(float* d_data, uword nx, uword ny);

/// @brief In-place 2D inverse FFT (float, interleaved).  Normalises by 1/(nx*ny).
void ifft2dAccelerate(float* d_data, uword nx, uword ny);

/// @brief In-place 3D forward FFT (float, interleaved).
void fft3dAccelerate(float* d_data, uword nx, uword ny, uword nz);

/// @brief In-place 3D inverse FFT (float, interleaved).  Normalises by 1/(nx*ny*nz).
void ifft3dAccelerate(float* d_data, uword nx, uword ny, uword nz);

// ---------------------------------------------------------------------------
// double overloads — always delegate to FFTW (vDSP double requires
// vDSP_fft2d_zipD which has limited GPU offload; also Metal uses float only)
// ---------------------------------------------------------------------------

/// @brief In-place 2D forward FFT (double, interleaved) — delegates to FFTW.
void fft2dAccelerate(double* d_data, uword nx, uword ny);

/// @brief In-place 2D inverse FFT (double, interleaved) — delegates to FFTW.
void ifft2dAccelerate(double* d_data, uword nx, uword ny);

/// @brief In-place 3D forward FFT (double, interleaved) — delegates to FFTW.
void fft3dAccelerate(double* d_data, uword nx, uword ny, uword nz);

/// @brief In-place 3D inverse FFT (double, interleaved) — delegates to FFTW.
void ifft3dAccelerate(double* d_data, uword nx, uword ny, uword nz);

#endif // METAL_COMPUTE
