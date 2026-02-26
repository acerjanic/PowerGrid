/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [Robject.h]

    Synopsis    [Base class of implementing quadratic and non-quadratic
                    regularization.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

/// @file Robject.h
/// @brief Abstract base class for regularization penalty functions.

#ifndef PowerGrid_Robject_h
#define PowerGrid_Robject_h

#include "Core/PGIncludes.h"
#include "Core/pgCol.hpp"

using namespace arma;
using namespace std;

/// @brief Abstract base class for regularization penalty functions.
///
/// Provides finite-difference operators Cd (forward) and Ctd (adjoint) along
/// with high-level helpers Penalty, Gradient, and Denom. Subclasses override
/// the potential functions wpot, dpot, and pot to implement specific penalties
/// (e.g. quadratic, total variation).
///
/// @tparam T1  Floating-point precision type (`float` or `double`).
template <typename T1> class Robject {
  typedef complex<T1> CxT1;

public:
  // Robject();
  /// @brief Default constructor.
  Robject(){};
  // Class members
  /// @brief Image x-dimension (number of pixels in x).
  uword Nx;
  /// @brief Image y-dimension (number of pixels in y).
  uword Ny;
  /// @brief Image z-dimension (number of pixels in z; 1 for 2-D).
  uword Nz;
  /// @brief Voxel spacing in x (unused by default, available for subclasses).
  T1 DeltaX;
  /// @brief Voxel spacing in y.
  T1 DeltaY;
  /// @brief Voxel spacing in z.
  T1 DeltaZ;
  /// @brief Regularization strength parameter β.
  T1 Beta;
  /// @brief Number of spatial dimensions to penalize (1, 2, or 3).
  uword Dims2Penalize;

  // It was declared as type Mat<uword> and the 3D type was a cube. We need to
  // vectorize it before it is passed to QuadPenalty.
  // Custom Class Constructor
  /// @brief Construct a regularization object.
  ///
  /// @param nx            Image x-dimension.
  /// @param ny            Image y-dimension.
  /// @param nz            Image z-dimension (use 1 for 2-D).
  /// @param beta          Regularization strength β.
  /// @param dims2penalize Number of spatial dimensions to penalize (default 3).
  Robject(uword nx, uword ny, uword nz, T1 beta, uword dims2penalize = 3);

  // Class Methods - Declared virtual so they can be implemented in the base
  // classes. Also they are virtual so that if you try to call Robject, things
  // crash rather than give un results.
  /// @brief Penalty weight function w(d) = ψ'(d)/d.
  ///
  /// Default implementation returns all-ones (quadratic).
  ///
  /// @param d  Finite-difference vector.
  /// @returns  Per-element weights, same length as @p d.
  virtual Col<CxT1> wpot(const Col<CxT1> &d) const {
    return ones<Col<CxT1>>(d.n_rows);
  }

  /// @brief Penalty derivative function ψ'(d).
  ///
  /// Default implementation returns @p d (quadratic derivative).
  ///
  /// @param d  Finite-difference vector.
  /// @returns  Per-element derivatives, same length as @p d.
  virtual Col<CxT1> dpot(const Col<CxT1> &d) const { return d; }

  /// @brief Penalty potential function ψ(d).
  ///
  /// Default implementation returns ½|d|² (quadratic potential).
  ///
  /// @param d  Finite-difference vector.
  /// @returns  Per-element potential values, same length as @p d.
  virtual Col<CxT1> pot(const Col<CxT1> &d) const {
    Col<T1> temp = abs(d) % abs(d) / 2.0;
    return conv_to<Col<CxT1>>::from(temp);
  }

  /// @brief Forward finite-difference operator along dimension @p dim.
  ///
  /// @param d    Input image vector of length Nx·Ny·Nz.
  /// @param dim  Spatial dimension (0=x, 1=y, 2=z).
  /// @returns    Finite-difference output, same length as @p d.
  Col<CxT1> Cd(const Col<CxT1> &d, uword dim) const;

  /// @brief Adjoint finite-difference operator along dimension @p dim.
  ///
  /// @param d    Input vector of length Nx·Ny·Nz.
  /// @param dim  Spatial dimension (0=x, 1=y, 2=z).
  /// @returns    Adjoint output, same length as @p d.
  Col<CxT1> Ctd(const Col<CxT1> &d, uword dim) const;

  /// @brief Evaluate the total penalty value R(x) = β · Σ_dim Σ_k ψ(Cd(x)[k]).
  ///
  /// @param x  Image vector of length Nx·Ny·Nz.
  /// @returns  Scalar penalty value.
  T1 Penalty(const Col<CxT1> &x) const;

  /// @brief Compute the gradient of the penalty ∇R(x).
  ///
  /// @param x  Image vector of length Nx·Ny·Nz.
  /// @returns  Gradient vector, same length as @p x.
  Col<CxT1> Gradient(const Col<CxT1> &x) const;

  /// @brief Compute the quadratic surrogate denominator for step-size selection.
  ///
  /// @param ddir  Search direction vector, length Nx·Ny·Nz.
  /// @param x     Current image estimate, same length.
  /// @returns     Scalar denominator value.
  CxT1 Denom(const Col<CxT1> &ddir, const Col<CxT1> &x) const;

  // --- pgCol overloads (Phase 3: eliminate arma boundary in PCG solver) ---

  /// @brief Forward finite-difference operator (pgCol overload).
  pgCol<pgComplex<T1>> Cd(const pgCol<pgComplex<T1>>& d, arma::uword dim) const;

  /// @brief Adjoint finite-difference operator (pgCol overload).
  pgCol<pgComplex<T1>> Ctd(const pgCol<pgComplex<T1>>& d, arma::uword dim) const;

  /// @brief Evaluate total penalty (pgCol overload).
  T1 Penalty(const pgCol<pgComplex<T1>>& x) const;

  /// @brief Compute gradient of penalty (pgCol overload).
  pgCol<pgComplex<T1>> Gradient(const pgCol<pgComplex<T1>>& x) const;

  /// @brief Compute quadratic surrogate denominator (pgCol overload).
  pgComplex<T1> Denom(const pgCol<pgComplex<T1>>& ddir, const pgCol<pgComplex<T1>>& x) const;

  // --- pgCol virtual potential functions ---
  /// @brief Penalty weight function (pgCol overload).
  virtual pgCol<pgComplex<T1>> wpot(const pgCol<pgComplex<T1>>& d) const {
      pgCol<pgComplex<T1>> out(d.n_elem);
      out.ones();
      return out;
  }

  /// @brief Penalty derivative function (pgCol overload).
  virtual pgCol<pgComplex<T1>> dpot(const pgCol<pgComplex<T1>>& d) const {
      return d;
  }

  /// @brief Penalty potential function (pgCol overload).
  virtual pgCol<pgComplex<T1>> pot(const pgCol<pgComplex<T1>>& d) const {
      // ½|d|² for quadratic default
      pgCol<T1> mag = abs(d);
      pgCol<T1> sq = mag % mag;
      sq /= T1(2);
      return to_complex(sq);
  }
};

// Explicit Instantiation
extern template class Robject<float>;
extern template class Robject<double>;

#endif // PowerGrid_Robject_h
