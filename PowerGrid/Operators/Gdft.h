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

/// @file Gdft.h
/// @brief Non-uniform field-corrected discrete Fourier transform operator.

#ifndef PowerGrid_Gdft_h
#define PowerGrid_Gdft_h

#include "Core/PGIncludes.h"
#include "FFT/ftCpu.h"
#include "Core/pgCol.hpp"
#include "Core/pgComplex.hpp"

using namespace arma;
using namespace std;

/// @brief Non-uniform field-corrected discrete Fourier transform (DFT) operator.
///
/// Computes the MRI signal model d = G * x where x is the image and d is the
/// k-space data. Supports optional off-resonance correction via a field map FM
/// and a per-sample time vector t.
///
/// Use `G * x` for the forward transform and `G / d` for the adjoint.
///
/// @tparam T1  Floating-point precision type (`float` or `double`).
template <typename T1>
class Gdft {
  typedef complex<T1> CxT1;

public:
  /// @brief Default constructor. Produces an uninitialized operator.
  Gdft();

  /// @brief Construct a Gdft operator with the given trajectory and field map.
  ///
  /// @param a   Number of k-space samples (n1, output size of forward transform).
  /// @param b   Number of image pixels (n2, input size of forward transform).
  /// @param k1  k-space x-coordinates, length @p a (units: cycles/FOV, range [-N/2, N/2)).
  /// @param k2  k-space y-coordinates, length @p a.
  /// @param k3  k-space z-coordinates, length @p a (set to zeros for 2-D).
  /// @param i1  Image-space x-coordinates, length @p b (units: fraction of FOV, range [0, 1)).
  /// @param i2  Image-space y-coordinates, length @p b.
  /// @param i3  Image-space z-coordinates, length @p b (set to zeros for 2-D).
  /// @param f1  Off-resonance field map, length @p b (units: rad/s).
  /// @param t1  Per-sample readout time vector, length @p a (units: s).
  Gdft(uword a, uword b, const Col<T1> &k1, const Col<T1> &k2,
       const Col<T1> &k3, const Col<T1> &i1, const Col<T1> &i2,
       const Col<T1> &i3, const Col<T1> &f1, const Col<T1> &t1);

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
  /// @brief Per-sample readout time vector in seconds (length n1).
  Col<T1> t;

#ifdef METAL_COMPUTE
  void* metalCtx = nullptr; // MetalDFTContext*, created for T1=float
  ~Gdft();
#endif

  /// @brief Forward DFT: image -> k-space.
  ///
  /// Computes d[j] = sum_k x[k] * exp(-i*2*pi*(kx[j]*ix[k] + ky[j]*iy[k] + kz[j]*iz[k]) - i*FM[k]*t[j]).
  ///
  /// @param d  Input image vector of length @p n2.
  /// @returns  Output k-space vector of length @p n1.
  Col<CxT1> operator*(const Col<CxT1> &d) const;

  /// @brief Adjoint DFT: k-space -> image.
  ///
  /// Computes the conjugate transpose of the forward operator (G^H).
  ///
  /// @param d  Input k-space vector of length @p n1.
  /// @returns  Output image vector of length @p n2.
  Col<CxT1> operator/(const Col<CxT1> &d) const;

  /// @brief Forward DFT (pgCol overload for Metal path).
  /// @param d  Input image vector of length @p n2.
  /// @returns  Output k-space vector of length @p n1.
  pgCol<pgComplex<T1>> operator*(const pgCol<pgComplex<T1>> &d) const;
  /// @brief Adjoint DFT (pgCol overload for Metal path).
  /// @param d  Input k-space vector of length @p n1.
  /// @returns  Output image vector of length @p n2.
  pgCol<pgComplex<T1>> operator/(const pgCol<pgComplex<T1>> &d) const;
};

extern template class Gdft<float>;
extern template class Gdft<double>;

#endif // PowerGrid_Gdft_h
