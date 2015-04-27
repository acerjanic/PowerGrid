// This file is distributed under a BSD 3-Clause license.
// See Ceemple_BSD_license.txt for details.

#ifndef CEEMPLECOMPLEX_H
#define CEEMPLECOMPLEX_H

#include <complex>
// complex, int
template <class T> std::complex<T> operator+(const std::complex<T> &x, int y) {
  return x + double(y);
}
template <class T> std::complex<T> operator+(int x, const std::complex<T> &y) {
  return double(x) + y;
}
template <class T> std::complex<T> operator-(const std::complex<T> &x, int y) {
  return x - double(y);
}
template <class T> std::complex<T> operator-(int x, const std::complex<T> &y) {
  return double(x) - y;
}
template <class T> std::complex<T> operator*(const std::complex<T> &x, int y) {
  return x * double(y);
}
template <class T> std::complex<T> operator*(int x, const std::complex<T> &y) {
  return double(x) * y;
}
template <class T> std::complex<T> operator/(const std::complex<T> &x, int y) {
  return x / double(y);
}
template <class T> std::complex<T> operator/(int x, const std::complex<T> &y) {
  return double(x) / y;
}

template <class T> std::complex<T> operator+(const std::complex<T> &x, float y) {
  return x + double(y);
}
template <class T> std::complex<T> operator+(float x, const std::complex<T> &y) {
  return double(x) + y;
}
template <class T> std::complex<T> operator-(const std::complex<T> &x, float y) {
  return x - double(y);
}
template <class T> std::complex<T> operator-(float x, const std::complex<T> &y) {
  return double(x) - y;
}
template <class T> std::complex<T> operator*(const std::complex<T> &x, float y) {
  return x * double(y);
}
template <class T> std::complex<T> operator*(float x, const std::complex<T> &y) {
  return double(x) * y;
}
template <class T> std::complex<T> operator/(const std::complex<T> &x, float y) {
  return x / double(y);
}
template <class T> std::complex<T> operator/(float x, const std::complex<T> &y) {
  return double(x) / y;
}

#endif // CEEMPLECOMPLEX_H
