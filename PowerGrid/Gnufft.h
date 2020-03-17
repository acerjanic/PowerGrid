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

#ifndef PowerGrid_Gnufft_h
#define PowerGrid_Gnufft_h

#include "PGIncludes.h"
#include "gridding.h"

using namespace arma;
using namespace std;

template <typename T1> // This is of type complex<double> or complex<float>, or
// any other type like float or single
class Gnufft {
  typedef complex<T1> CxT1;

public:
  // Default Class Constructor and Destructor
  Gnufft();

  //~Ggrid();
  // Class Constructor
  Gnufft(uword dataLength, T1 gridos, uword nx, uword ny, uword nz,
         const Col<T1> &k1, const Col<T1> &k2, const Col<T1> &k3,
         const Col<T1> &i1, const Col<T1> &i2, const Col<T1> &i3);

  // Class destructor to free LUT
  ~Gnufft();

  // Class variables go here. Change as necessary
  uword n1 = 0;
  uword n2 = 0;
  uword Nx = 0;
  uword Ny = 0;
  uword Nz = 0;

  T1 *kx, *ky, *kz; // kspace coordinates
  Col<T1> ix; // image space coordinates
  Col<T1> iy;
  Col<T1> iz;

  T1 gridOS;   // grid oversampling
  T1 *LUT = 0; // Lookup table for the gridding operations
  uword sizeLUT = 0;
  T1 beta; // beta factor for gridding not the same as beta in regularization!
  T1 kernelWidth; // Kaiser Bessel Kernel Support
  void *stream;
  
  #ifdef OPENACC_GPU
    cufftHandle plan;
  #else
    void* plan;
  #endif

  uword imageNumElems;
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
  // Overloaded methods for forward and adjoint transform
  // Forward transform operation using gridding
  Col<CxT1> operator*(const Col<CxT1> &d) const;
  // Adjoint transform operation
  Col<CxT1> operator/(const Col<CxT1> &d) const;

  Col<CxT1> forwardSpatialInterp(const Col<CxT1> &d) const;
  // Adjoint transform operation
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
