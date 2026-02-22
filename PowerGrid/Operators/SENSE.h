/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [SENSE.h]

    Synopsis    [Object implementing sensitivity encoding reconstructions. ]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

/// @file SENSE.h
/// @brief Multi-coil sensitivity encoding (SENSE) operator.

#ifndef PowerGrid_SENSE_hpp
#define PowerGrid_SENSE_hpp

#include "Gdft.h"
#include "GdftR2.h"
#include "Gnufft.h"
#include "Core/PGIncludes.h"
#include "Gridding/TimeSegmentation.h"
#include "Core/pgCol.hpp"
#include "Core/pgMat.hpp"

using namespace std;
using namespace arma;

/// @brief Multi-coil sensitivity encoding (SENSE) operator.
///
/// Wraps any encoding operator @p Tobj with per-coil sensitivity maps to form
/// the full multi-coil signal model:
///
///   d_c = G * (S_c . x)   for each coil c
///
/// where S_c is the sensitivity map for coil c and . denotes element-wise
/// multiplication. All coil data is stacked into a single vector of length n1*nc.
///
/// Use `S * x` for the forward transform and `S / d` for the adjoint.
///
/// @tparam T1    Floating-point precision type (`float` or `double`).
/// @tparam Tobj  Encoding operator type (e.g., `Gdft<T1>`, `Gnufft<T1>`,
///               `TimeSegmentation<T1, Gnufft<T1>>`).
template <typename T1, typename Tobj> class SENSE {
  typedef complex<T1> CxT1;
  typedef Col<CxT1> ColCxT1;
  typedef Mat<CxT1> MatCxT1;

public:
  /// @brief Default constructor. Produces an uninitialized operator.
  SENSE();

  /// @brief Number of k-space samples per coil (output length per coil).
  uword n1 = 0;
  /// @brief Number of image pixels (input length of forward transform).
  uword n2 = 0;
  /// @brief Number of receiver coils.
  uword nc = 0;
  /// @brief Pointer to the underlying single-coil encoding operator.
  Tobj *G_obj;
  /// @brief Sensitivity maps matrix, size n2 x nc (one column per coil).
  Mat<CxT1> SMap;
  /// @brief Conjugate sensitivity maps matrix, size n2 x nc.
  Mat<CxT1> conjSMap;

  mutable Mat<CxT1> outData;
  mutable Mat<CxT1> outImg;
  mutable Mat<CxT1> coilWeightData;
  mutable Mat<CxT1> coilWeightImg;
  mutable Mat<CxT1> coilImages;

#ifdef METAL_COMPUTE
  // pgMat copies of sensitivity maps for Metal GPU dispatch (float only).
  pgMat<pgComplex<T1>> SMap_pg;
  pgMat<pgComplex<T1>> conjSMap_pg;
  mutable pgMat<pgComplex<T1>> outData_pg;
  mutable pgMat<pgComplex<T1>> outImg_pg;
#endif

  /// @brief Construct a SENSE operator from a single-coil operator and sensitivity maps.
  ///
  /// @param G        Single-coil encoding operator (e.g., Gdft or Gnufft).
  /// @param SENSEmap Sensitivity maps as a flat column vector of length n2*nc,
  ///                 ordered as [coil0_pixel0, coil0_pixel1, ..., coil1_pixel0, ...].
  /// @param a        Number of k-space samples per coil (n1).
  /// @param b        Number of image pixels (n2).
  /// @param c        Number of receiver coils (nc).
  SENSE(Tobj &G, Col<CxT1> SENSEmap, uword a, uword b, uword c);

  /// @brief Forward SENSE transform: image -> stacked multi-coil k-space.
  ///
  /// Returns a vector of length n1*nc with coil data concatenated as
  /// [coil0_data; coil1_data; ...; coil_{nc-1}_data].
  ///
  /// @param d  Input image vector of length @p n2.
  /// @returns  Output k-space vector of length n1*nc.
  Col<CxT1> operator*(const Col<CxT1> &d) const;

  /// @brief Adjoint SENSE transform: stacked multi-coil k-space -> image.
  ///
  /// Applies the adjoint of each coil operator, weights by the conjugate
  /// sensitivity map, and sums across coils.
  ///
  /// @param d  Input k-space vector of length n1*nc.
  /// @returns  Output image vector of length @p n2.
  Col<CxT1> operator/(const Col<CxT1> &d) const;

  /// @brief Forward SENSE transform (pgCol overload for Metal path).
  /// @param d  Input image vector of length @p n2.
  /// @returns  Output k-space vector of length n1*nc.
  pgCol<pgComplex<T1>> operator*(const pgCol<pgComplex<T1>> &d) const;
  /// @brief Adjoint SENSE transform (pgCol overload for Metal path).
  /// @param d  Input k-space vector of length n1*nc.
  /// @returns  Output image vector of length @p n2.
  pgCol<pgComplex<T1>> operator/(const pgCol<pgComplex<T1>> &d) const;
};

// Explicit Instantiations
extern template class SENSE<float, Gnufft<float>>;
extern template class SENSE<float, TimeSegmentation<float, Gnufft<float>>>;
extern template class SENSE<double, Gnufft<double>>;
extern template class SENSE<double, TimeSegmentation<double, Gnufft<double>>>;
extern template class SENSE<float, Gdft<float>>;
extern template class SENSE<double, Gdft<double>>;
extern template class SENSE<float, GdftR2<float>>;
extern template class SENSE<double, GdftR2<double>>;
#endif
