/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [TimeSegmentation.h]

    Synopsis    [Wrappers to the cuFFT library supporting single and double
                                precision for GPU accelerated FFTs]

    Description []

    Revision    [0.1.0; Giang-Chau Ngo, BIOE UIUC
                 0.2.0; Alex Cerjanic, BIOE UIUC]

    Date        [12/2/2016]

 *****************************************************************************/
/// @file TimeSegmentation.h
/// @brief Off-resonance correction via time-segmented interpolation.

// We are using two template types at the moment. One for the type of data to be
// processed (ie Col<cx_double>) and one for the type of G object (ie
// Gfft<Col<cx_double>>
// T1 is the data type for complex, T2 is the data type for real data

#ifndef PowerGrid_TimeSegmentation_h
#define PowerGrid_TimeSegmentation_h

#include "Operators/Gdft.h"
#include "Operators/Gnufft.h"
#include "Core/PGIncludes.h"
#include "Core/pgCol.hpp"
#include "Core/pgMat.hpp"
#include <vector>

using namespace std;

/// @brief Off-resonance correction via time-segmented interpolation.
///
/// Wraps any encoding operator @p Tobj with field-map-based off-resonance
/// correction using the time-segmentation approximation.  The field map and
/// per-sample timing vector are used to construct L time-segment sub-problems,
/// each solved with the underlying operator, then combined via interpolation
/// coefficients AA.
///
/// @tparam T1    Floating-point precision type (`float` or `double`).
/// @tparam Tobj  Underlying encoding operator type (e.g., `Gnufft<T1>`).
template <typename T1, typename Tobj> class TimeSegmentation {
  typedef complex<T1> CxT1;

public:
  /// @brief Default constructor.
  TimeSegmentation();

  /// @brief Number of k-space samples (output length of forward transform).
  uword n1;
  /// @brief Number of image pixels (input length of forward transform).
  uword n2;
  /// @brief Number of time segments L.
  int L;
  /// @brief Interpolation type: 1 = Hanning, 2 = exact min-max, 3 = approx min-max.
  uword type;
  /// @brief Number of shots; used to reduce complexity of interpolator calculation.
  uword Nshots;
  /// @brief Time-segment length tau (seconds).
  T1 tau;
  /// @brief Minimum time in the time vector (e.g., TE for spiral-out), in seconds.
  T1 T_min;
  /// @brief Pointer to the underlying single-coil encoding operator.
  Tobj *obj;
  /// @brief Off-resonance field map (rad/s), length n2.
  Col<T1> fieldMap;
  /// @brief Per-sample readout time vector (s), length n1.
  Col<T1> timeVec;
  /// @brief Interpolation coefficient matrix, size L x n1.
  Mat<CxT1> AA;
  /// @brief Cached conj(AA) for adjoint operator — avoids recomputing per call.
  Mat<CxT1> conjAA;
  CxT1 i = CxT1(0., 1.);
  /// @brief Phase modulation matrices Wo (n2 x L) for forward transform.
  Mat<CxT1> Wo;
  /// @brief Conjugate phase modulation matrices WoH (n2 x L) for adjoint transform.
  Mat<CxT1> WoH;
  Col<T1> RowOnes;
  mutable Mat<complex<T1>> outData;
  mutable Mat<complex<T1>> outImg;
  mutable Mat<complex<T1>> tempD;
  mutable Mat<complex<T1>> tempAD;

#ifdef METAL_COMPUTE
  // Page-aligned column vectors for Metal GPU dispatch (float only).
  // Stored as separate pgCols rather than pgMat so every column is page-aligned
  // (16384 bytes), enabling Metal zero-copy GPU access for all element-wise ops.
  std::vector<pgCol<pgComplex<T1>>> Wo_cols;      // L columns from Wo
  std::vector<pgCol<pgComplex<T1>>> WoH_cols;     // L columns from WoH
  std::vector<pgCol<pgComplex<T1>>> AA_cols;      // L columns from AA
  std::vector<pgCol<pgComplex<T1>>> conjAA_cols;  // L columns from conj(AA)
  mutable std::vector<pgCol<pgComplex<T1>>> tempD_cols;   // L working buffers (n2 each)
  mutable std::vector<pgCol<pgComplex<T1>>> tempAD_cols;  // L working buffers (n1 each)
#endif

  /// @brief Construct a TimeSegmentation operator.
  ///
  /// @param G          Reference to the underlying encoding operator.
  /// @param map_in     Off-resonance field map in rad/s, length n2.
  /// @param timeVec_in Per-sample readout times in seconds, length n1.
  /// @param a          Number of k-space samples (n1).
  /// @param b          Number of image pixels (n2).
  /// @param c          Number of time segments (L).
  /// @param interptype Interpolation method: 1=Hanning, 2=exact min-max, 3=approx min-max (default 1).
  /// @param shots      Number of shots for segmentation (default 1).
  TimeSegmentation(Tobj &G, Col<T1> map_in, Col<T1> timeVec_in, uword a,
                   uword b, uword c, uword interptype = 1, uword shots = 1);

  /// @brief Forward transform: image -> k-space with off-resonance correction.
  ///
  /// @param d  Input image vector of length @p n2.
  /// @returns  Output k-space vector of length @p n1.
  Col<CxT1> operator*(const Col<CxT1> &d) const;

  /// @brief Adjoint transform: k-space -> image with off-resonance correction.
  ///
  /// @param d  Input k-space vector of length @p n1.
  /// @returns  Output image vector of length @p n2.
  Col<CxT1> operator/(const Col<CxT1> &d) const;

  // pgCol overloads — avoid arma conversion overhead
  pgCol<pgComplex<T1>> operator*(const pgCol<pgComplex<T1>> &d) const;
  pgCol<pgComplex<T1>> operator/(const pgCol<pgComplex<T1>> &d) const;

  protected:

  /// @brief Compute 1-D FFT along the time dimension for interpolator construction.
  ///
  /// @param d   Input complex vector.
  /// @param KK  FFT length.
  /// @returns   FFT output vector of length @p KK.
  Col<CxT1> calcFFT1D(const Col<CxT1> &d, uword KK) const;
};

// Now we insert the explicit instantiations we need
extern template class TimeSegmentation<float, Gnufft<float>>;
extern template class TimeSegmentation<float, Gdft<float>>;
extern template class TimeSegmentation<double, Gnufft<double>>;
extern template class TimeSegmentation<double, Gdft<double>>;
#endif // PowerGrid_TimeSegmentation_h
