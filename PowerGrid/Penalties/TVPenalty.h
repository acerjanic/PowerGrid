/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [TVPenalty.h]

    Synopsis    [Implementation of an approximate total variational penalty.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

/// @file TVPenalty.h
/// @brief Approximate total-variation regularization penalty.

#ifndef PowerGrid_TVPenalty_h
#define PowerGrid_TVPenalty_h

#include "Robject.h"

using namespace arma;

/// @brief Approximate total-variation (TV) regularization penalty.
///
/// Implements a smooth approximation to TV using the Fair / Huber-like
/// potential ψ(d) = |d| − δ·log(1 + |d|/δ) with smoothing parameter δ.
/// This avoids the non-differentiability of exact TV at d=0 while
/// preserving edge-preserving properties for large |d|.
///
/// @tparam T1  Floating-point precision type (`float` or `double`).
template <typename T1> class TVPenalty : public Robject<T1> {
  typedef complex<T1> CxT1;

public:
  /// @brief Default constructor.
  TVPenalty();

  // It was declared as type Mat<uword> and the 3D type was a cube. We need to
  // vectorize it before it is passed to QuadPenalty.
  // Custom Class Constructor
  /// @brief Construct an approximate TV penalty operator.
  ///
  /// @param nx            Image x-dimension.
  /// @param ny            Image y-dimension.
  /// @param nz            Image z-dimension (use 1 for 2-D).
  /// @param beta          Regularization strength β.
  /// @param delta         Smoothing parameter δ > 0; smaller values approach exact TV.
  /// @param dims2penalize Number of spatial dimensions to penalize (default 3).
  TVPenalty(uword nx, uword ny, uword nz, T1 beta, T1 delta, uword dims2penalize = 3);

  // Class Methods

  /// @brief TV weight function wpot(d) = 1 / (1 + |d|/δ).
  ///
  /// @param d  Finite-difference vector.
  /// @returns  Per-element weights, same length as @p d.
  Col<CxT1> wpot(const Col<CxT1> &d) const;

  /// @brief TV derivative dpot(d) = d / (1 + |d|/δ).
  ///
  /// @param d  Finite-difference vector.
  /// @returns  Per-element derivatives, same length as @p d.
  virtual Col<CxT1> dpot(const Col<CxT1> &d) const;

  /// @brief TV potential pot(d) = |d| − δ·log(1 + |d|/δ).
  ///
  /// @param d  Finite-difference vector.
  /// @returns  Per-element potential values, same length as @p d.
  Col<CxT1> pot(const Col<CxT1> &d) const;

private:
  /// @brief Smoothing parameter δ for the Fair potential approximation.
  T1 Delta;
};

// Explicit Instantiation
extern template class TVPenalty<double>;
extern template class TVPenalty<float>;

#endif // PowerGrid_TVPenalty_h
