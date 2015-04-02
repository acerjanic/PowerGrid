// This file is distributed under a BSD 3-Clause license.
// See Ceemple_BSD_license.txt for details.

#ifndef CEEMPLEARMADILLO_H
#define CEEMPLEARMADILLO_H

#include "CeempleComplex.h"
#include "armadillo"

template <typename T>
inline arma::Mat<T> join_cols(const arma::Mat<T> &A, const double b) {
  arma::Mat<T> B;
  B = (T)b;

  return join_cols(A, B);
}

template <typename T>
inline arma::Mat<T> join_cols(const double a, const arma::Mat<T> &B) {
  arma::Mat<T> A;
  A = (T)a;

  return join_cols(A, B);
}

template <typename T>
inline arma::Mat<T> join_rows(const arma::Mat<T> &A, const double b) {
  arma::Mat<T> B;
  B = (T)b;

  return join_rows(A, B);
}

template <typename T>
inline arma::Mat<T> join_rows(const double a, const arma::Mat<T> &B) {
  arma::Mat<T> A;
  A = (T)a;

  return join_rows(A, B);
}

template <typename T> inline arma::uword size(const arma::Mat<T> &A, double d) {
  if (d == 1)
    return A.n_rows;
  else if (d == 2)
    return A.n_cols;
  else
    return 0;
}

#endif // CEEMPLEARMADILLO_H
