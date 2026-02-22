/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [Gdft.h]

    Synopsis    [Object that represents a non-uniform field corrected discrete
                    Fourier tranform.]

    Description [Forward transforms are denoted by G*data and adjoint transforms
                    are denoted by G/data. See documentation for more
                    information]

    Revision    [0.2.0; Alex Cerjanic, BIOE UIUC]

    Date        [12/2/2016]

 *****************************************************************************/

/// @file GdftR2.h
/// @brief Non-uniform DFT operator with R2* field map gradient support.

#ifndef PowerGrid_GdftR2_h
#define PowerGrid_GdftR2_h

#include "Core/PGIncludes.h"
#include "FFT/ftCpuWithGrads.h"
#include "Core/pgCol.hpp"
#include "Core/pgComplex.hpp"

using namespace arma;
using namespace std;

/// @brief Non-uniform DFT operator extended with R2* (transverse relaxation) gradient maps.
///
/// Extends Gdft with R2* field map gradients (Gx, Gy, Gz), enabling first-order
/// Taylor expansion of the off-resonance phase for improved field-corrected
/// reconstruction accuracy.
///
/// Use `G * x` for the forward transform and `G / d` for the adjoint.
///
/// @tparam T1  Floating-point precision type (`float` or `double`).
template <typename T1>
class GdftR2 {
  typedef complex<T1> CxT1;

public:
  /// @brief Default constructor. Produces an uninitialized operator.
  GdftR2();

  /// @brief Construct a GdftR2 operator with trajectory, field map, and image dimensions.
  ///
  /// @param a     Number of k-space samples (n1, output size of forward transform).
  /// @param b     Number of image pixels (n2, input size of forward transform).
  /// @param k1    k-space x-coordinates, length @p a (units: cycles/FOV).
  /// @param k2    k-space y-coordinates, length @p a.
  /// @param k3    k-space z-coordinates, length @p a.
  /// @param i1    Image-space x-coordinates, length @p b (units: fraction of FOV).
  /// @param i2    Image-space y-coordinates, length @p b.
  /// @param i3    Image-space z-coordinates, length @p b.
  /// @param f1    Off-resonance field map, length @p b (units: rad/s).
  /// @param t1    Per-sample readout time vector, length @p a (units: s).
  /// @param numX  Image x-dimension (pixels).
  /// @param numY  Image y-dimension (pixels).
  /// @param numZ  Image z-dimension (pixels; use 1 for 2-D).
  GdftR2(uword a, uword b, const Col<T1> &k1, const Col<T1> &k2,
       const Col<T1> &k3, const Col<T1> &i1, const Col<T1> &i2,
       const Col<T1> &i3, const Col<T1> &f1, const Col<T1> &t1,
       const int numX, const int numY, const int numZ);

  /// @brief Number of k-space samples (output length of forward transform).
  uword n1 = 0;
  /// @brief Number of image pixels (input length of forward transform).
  uword n2 = 0;

  /// @brief k-space x-coordinates (length n1).
  Col<T1> kx;
  /// @brief k-space y-coordinates (length n1).
  Col<T1> ky;
  /// @brief k-space z-coordinates (length n1).
  Col<T1> kz;
  /// @brief Image-space x-coordinates (length n2).
  Col<T1> ix;
  /// @brief Image-space y-coordinates (length n2).
  Col<T1> iy;
  /// @brief Image-space z-coordinates (length n2).
  Col<T1> iz;
  /// @brief Off-resonance field map in rad/s (length n2).
  Col<T1> FM;
  /// @brief R2* gradient map in the x-direction (length n2).
  Col<T1> Gx;
  /// @brief R2* gradient map in the y-direction (length n2).
  Col<T1> Gy;
  /// @brief R2* gradient map in the z-direction (length n2).
  Col<T1> Gz;
  /// @brief Per-sample readout time vector in seconds (length n1).
  Col<T1> t;

  /// @brief Image x-dimension in pixels.
  /// @brief Image y-dimension in pixels.
  /// @brief Image z-dimension in pixels.
  int numX, numY, numZ;

#ifdef METAL_COMPUTE
  void* metalCtx = nullptr; // MetalDFTContext*, created for T1=float
  ~GdftR2();
#endif

  /// @brief Forward DFT with R2* gradients: image -> k-space.
  ///
  /// @param d  Input image vector of length @p n2.
  /// @returns  Output k-space vector of length @p n1.
  Col<CxT1> operator*(const Col<CxT1> &d) const;

  /// @brief Adjoint DFT with R2* gradients: k-space -> image.
  ///
  /// @param d  Input k-space vector of length @p n1.
  /// @returns  Output image vector of length @p n2.
  Col<CxT1> operator/(const Col<CxT1> &d) const;

  /// @brief Compute R2* gradient maps from the field map using finite differences.
  ///
  /// Populates @p Gx, @p Gy, and @p Gz with the spatial gradients of the field map.
  ///
  /// @param Gx  Output: x-gradient of the field map (length n2).
  /// @param Gy  Output: y-gradient of the field map (length n2).
  /// @param Gz  Output: z-gradient of the field map (length n2).
  void calcGradientMaps(Col<T1> &Gx, Col<T1> &Gy, Col<T1> &Gz);

  /// @brief Forward finite difference operator along dimension @p dim.
  ///
  /// @param d    Input vector.
  /// @param dim  Spatial dimension (0 = x, 1 = y, 2 = z).
  /// @returns    Forward difference of @p d along @p dim.
  Col<T1> Cd(const Col<T1> &d, uword dim) const;

  /// @brief Forward DFT with R2* gradients (pgCol overload for Metal path).
  /// @param d  Input image vector of length @p n2.
  /// @returns  Output k-space vector of length @p n1.
  pgCol<pgComplex<T1>> operator*(const pgCol<pgComplex<T1>> &d) const;
  /// @brief Adjoint DFT with R2* gradients (pgCol overload for Metal path).
  /// @param d  Input k-space vector of length @p n1.
  /// @returns  Output image vector of length @p n2.
  pgCol<pgComplex<T1>> operator/(const pgCol<pgComplex<T1>> &d) const;
};

extern template class GdftR2<float>;
extern template class GdftR2<double>;

#endif // PowerGrid_GdftR2_h
