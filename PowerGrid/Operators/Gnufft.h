/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [Gnufft.h]

    Synopsis    [Object that represents a non-uniform discrete Fourier
                    transform implemented via a a CPU and GPU accelerated
                    gridding implementation of the non-uniform Fast Fourier
                    Tranform (NUFFT).]

    Description [Forward transforms are denoted by G*data and adjoint transforms
                    are denoted by G/data. See documentation for more
                    information]

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

/// @file Gnufft.h
/// @brief Non-uniform Fast Fourier Transform (NUFFT) operator via Kaiser-Bessel gridding.

#ifndef PowerGrid_Gnufft_h
#define PowerGrid_Gnufft_h

#include "Core/PGIncludes.h"
#include "Gridding/gridding.h"
#include "Core/pgCol.hpp"
#include "Core/pgComplex.hpp"

using namespace arma;
using namespace std;

/// @brief Non-uniform Fast Fourier Transform operator using Kaiser-Bessel gridding.
///
/// Implements an efficient NUFFT via convolution with a Kaiser-Bessel (KB) kernel
/// onto an oversampled Cartesian grid, followed by an FFT. Supports both CPU and
/// GPU execution (via OpenACC). An LUT (look-up table) is used for fast kernel
/// evaluation.
///
/// k-space coordinates must lie within [-N/2, N/2) in each dimension.
///
/// Use `G * x` for the forward transform and `G / d` for the adjoint.
///
/// @tparam T1  Floating-point precision type (`float` or `double`).
template <typename T1>
class Gnufft {
  typedef complex<T1> CxT1;

public:
  /// @brief Default constructor. Produces an uninitialized operator.
  Gnufft();

  /// @brief Construct a Gnufft operator with the given trajectory and image dimensions.
  ///
  /// @param dataLength  Number of k-space samples (n2, output of forward transform).
  /// @param gridos      Grid oversampling factor (typically 2.0).
  /// @param nx          Image x-dimension in pixels.
  /// @param ny          Image y-dimension in pixels.
  /// @param nz          Image z-dimension in pixels (use 1 for 2-D).
  /// @param k1          k-space x-coordinates, length @p dataLength (units: cycles, range [-nx/2, nx/2)).
  /// @param k2          k-space y-coordinates, length @p dataLength.
  /// @param k3          k-space z-coordinates, length @p dataLength.
  /// @param i1          Image-space x-coordinates, length nx*ny*nz (not used in gridding, stored for reference).
  /// @param i2          Image-space y-coordinates, length nx*ny*nz.
  /// @param i3          Image-space z-coordinates, length nx*ny*nz.
  Gnufft(uword dataLength, T1 gridos, uword nx, uword ny, uword nz,
         const Col<T1> &k1, const Col<T1> &k2, const Col<T1> &k3,
         const Col<T1> &i1, const Col<T1> &i2, const Col<T1> &i3);

  /// @brief Destructor. Frees the KB kernel LUT and all internal grid buffers.
  ~Gnufft();

  /// @brief Number of image pixels (n1 = Nx * Ny * Nz).
  uword n1 = 0;
  /// @brief Number of k-space samples.
  uword n2 = 0;
  /// @brief Image x-dimension in pixels.
  uword Nx = 0;
  /// @brief Image y-dimension in pixels.
  uword Ny = 0;
  /// @brief Image z-dimension in pixels.
  uword Nz = 0;

  /// @brief k-space x-coordinates array (length n2, range [-Nx/2, Nx/2)).
  /// @brief k-space y-coordinates array (length n2, range [-Ny/2, Ny/2)).
  /// @brief k-space z-coordinates array (length n2, range [-Nz/2, Nz/2)).
  T1 *kx, *ky, *kz;
  /// @brief Image-space x-coordinates (length n1, stored for reference).
  Col<T1> ix;
  /// @brief Image-space y-coordinates (length n1, stored for reference).
  Col<T1> iy;
  /// @brief Image-space z-coordinates (length n1, stored for reference).
  Col<T1> iz;

  /// @brief Grid oversampling factor.
  T1 gridOS;
  /// @brief Kaiser-Bessel kernel look-up table.
  T1 *LUT = 0;
  /// @brief Number of entries in the KB kernel LUT.
  uword sizeLUT = 0;
  /// @brief KB kernel shape parameter (computed from gridOS and kernelWidth).
  T1 beta;
  /// @brief Kaiser-Bessel kernel half-width in grid units.
  T1 kernelWidth;
  /// @brief GPU stream handle (OpenACC/CUDA, NULL on CPU builds).
  void *stream;

  #ifdef OPENACC_GPU
    cufftHandle plan;
  #else
    void* plan;
  #endif

  #ifdef METAL_COMPUTE
    /// @brief Metal gridding context (legacy); nullptr when T1 != float or Metal unavailable.
    MetalGriddingContext* metalCtx = nullptr;
    /// @brief Full GPU NUFFT pipeline context (macOS 14+); nullptr if unavailable.
    MetalNufftPipelineContext* pipelineCtx = nullptr;
  #endif

  /// @brief Total number of elements in the (non-oversampled) image grid.
  uword imageNumElems;
  /// @brief Total number of elements in the oversampled grid.
  uword gridNumElems;

  T1 *pGridData, *pGridData_d, *pGridData_os, *pGridData_os_d;
  T1 *pSamples;
  complex<T1> *gridData, *gridData_d, *gridData_os, *gridData_os_d;
  complex<T1> *samples;

  mutable Col<CxT1> XformedData;
  mutable Col<CxT1> XformedImg;  
  mutable Col<T1> realXformedData;
  mutable Col<T1> imagXformedData;
  mutable Col<T1> realXformedImg;
  mutable Col<T1> imagXformedImg;
  /// @brief Forward NUFFT: image -> k-space via KB gridding.
  ///
  /// Applies deapodization, zero-pads to the oversampled grid, FFT-shifts,
  /// performs an FFT, and interpolates at the non-uniform k-space locations.
  ///
  /// @param d  Input image vector of length @p n1.
  /// @returns  Output k-space vector of length @p n2.
  Col<CxT1> operator*(const Col<CxT1> &d) const;

  /// @brief Adjoint NUFFT: k-space -> image via KB gridding.
  ///
  /// Spreads k-space samples onto the oversampled grid using the KB kernel,
  /// performs an inverse FFT, crops to image size, and applies deapodization.
  ///
  /// @param d  Input k-space vector of length @p n2.
  /// @returns  Output image vector of length @p n1.
  Col<CxT1> operator/(const Col<CxT1> &d) const;

  /// @brief Forward NUFFT (pgCol overload for Metal path).
  /// @param d  Input image vector of length @p n1.
  /// @returns  Output k-space vector of length @p n2.
  pgCol<pgComplex<T1>> operator*(const pgCol<pgComplex<T1>> &d) const;
  /// @brief Adjoint NUFFT (pgCol overload for Metal path).
  /// @param d  Input k-space vector of length @p n2.
  /// @returns  Output image vector of length @p n1.
  pgCol<pgComplex<T1>> operator/(const pgCol<pgComplex<T1>> &d) const;

private:
#ifdef METAL_COMPUTE
  /// Metal forward pipeline: deapodize → zero-pad → FFT → gridding.
  /// Writes result to pSamples (n2 complex elements).
  void metalForwardImpl(const T1* dataPtr) const;

  /// Metal adjoint pipeline: gridding → IFFT → crop → deapodize.
  /// Writes result to pGridData (n1 complex elements).
  void metalAdjointImpl(const T1* dataPtr) const;
#endif

public:

  /// @brief Forward spatial interpolation only (no deapodization or FFT).
  ///
  /// @param d  Input image vector of length @p n1.
  /// @returns  Spatially interpolated k-space vector of length @p n2.
  Col<CxT1> forwardSpatialInterp(const Col<CxT1> &d) const;

  /// @brief Adjoint spatial interpolation only (no deapodization or FFT).
  ///
  /// @param d  Input k-space vector of length @p n2.
  /// @returns  Output image vector of length @p n1.
  Col<CxT1> adjointSpatialInterp(const Col<CxT1> &d) const;
/*
  Col<CxT1> trimmedForwardOp(const Col<CxT1> &d,
                             const Col<CxT1> &tempInterp) const;

  // Adjoint transform operation
  Col<CxT1> trimmedAdjointOp(const Col<CxT1> &d,
                             const Col<CxT1> &tempInterp) const;
                             */
};

// Explicit Instantiation
extern template class Gnufft<float>;
extern template class Gnufft<double>;

#endif // PowerGrid_Gnufft_h
