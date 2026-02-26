/*
   (C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
   All rights reserved.

   See LICENSE.txt for the University of Illinois/NCSA Open Source license.

   Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
 */

/*****************************************************************************

    File Name   [QuadPenalty.h]

    Synopsis    [Quadratic Penalty using Robject interface to maintain
                    commonality with the Image Reconstruction Toolkit (IRT).]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

*****************************************************************************/

/// @file QuadPenalty.h
/// @brief Quadratic (Tikhonov) regularization penalty.

#ifndef PowerGrid_QuadPenalty_h
#define PowerGrid_QuadPenalty_h

#include "Robject.h"

using namespace arma;

/// @brief Quadratic (Tikhonov) regularization penalty.
///
/// Implements the penalty R(x) = β · ½ ‖C·x‖² where C is the finite-
/// difference operator inherited from Robject.  All three potential functions
/// correspond to the quadratic case: wpot returns all-ones, dpot returns the
/// identity, and pot returns ½|d|².
///
/// @tparam T1  Floating-point precision type (`float` or `double`).
template <typename T1> class QuadPenalty : public Robject<T1> {
typedef complex<T1> CxT1;

public:
/// @brief Default constructor.
QuadPenalty();

// It was declared as type Mat<uword> and the 3D type was a cube. We need to
// vectorize it before it is passed to QuadPenalty.
// Custom Class Constructor
/// @brief Construct a quadratic penalty operator.
///
/// @param nx            Image x-dimension.
/// @param ny            Image y-dimension.
/// @param nz            Image z-dimension (use 1 for 2-D).
/// @param beta          Regularization strength β.
/// @param dims2penalize Number of spatial dimensions to penalize (default 3).
QuadPenalty(uword nx, uword ny, uword nz, T1 beta, uword dims2penalize = 3);

// Class Methods

/// @brief Quadratic weight function: returns all-ones (constant weight).
///
/// @param d  Finite-difference vector.
/// @returns  Vector of ones, same length as @p d.
Col<CxT1> wpot(const Col<CxT1> &d) const;

/// @brief Quadratic derivative: returns @p d (identity).
///
/// @param d  Finite-difference vector.
/// @returns  @p d unchanged.
Col<CxT1> dpot(const Col<CxT1> &d) const;
/// @brief Quadratic potential: returns ½|d|².
///
/// @param d  Finite-difference vector.
/// @returns  Per-element values ½|d[k]|².
Col<CxT1> pot(const Col<CxT1> &d) const;

// pgCol overloads
pgCol<pgComplex<T1>> wpot(const pgCol<pgComplex<T1>>& d) const override;
pgCol<pgComplex<T1>> dpot(const pgCol<pgComplex<T1>>& d) const override;
pgCol<pgComplex<T1>> pot(const pgCol<pgComplex<T1>>& d) const override;
};

// Explicit Instantiation
extern template class QuadPenalty<float>;
extern template class QuadPenalty<double>;

#endif
