/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [Gdft.cpp]

    Synopsis    [Object that represents a non-uniform field corrected discrete
                    Fourier tranform.]

    Description [Forward transforms are denoted by G*data and adjoint transforms
                    are denoted by G/data. See documentation for more
                    information]

    Revision    [0.2.0; Alex Cerjanic, BIOE UIUC]

    Date        [12/2/2016]

 *****************************************************************************/
#include "GdftR2.h"

#ifdef METAL_COMPUTE
#include "Metal/MetalDFT.h"
#endif

using namespace arma;

template <typename T1>
GdftR2<T1>::GdftR2(
    uword a, uword b, const Col<T1> &k1, const Col<T1> &k2, const Col<T1> &k3,
    const Col<T1> &i1, const Col<T1> &i2, const Col<T1> &i3, const Col<T1> &f1,
    const Col<T1> &t1, const int num_x, const int num_y, const int num_z)
{

  n1 = a;
  n2 = b;
  kx = k1;
  ky = k2;
  kz = k3;
  ix = i1;
  iy = i2;
  iz = i3;
  FM = f1;
  t = t1;

  numX = num_x;
  numY = num_y;
  numZ = num_z;

  calcGradientMaps(Gx, Gy, Gz);

#ifdef METAL_COMPUTE
  if constexpr (std::is_same<T1, float>::value) {
    metalCtx = metal_dft_create_with_grads(
        kx.memptr(), ky.memptr(), kz.memptr(),
        ix.memptr(), iy.memptr(), iz.memptr(),
        FM.memptr(), t.memptr(),
        Gx.memptr(), Gy.memptr(), Gz.memptr(),
        (unsigned int)n1, (unsigned int)n2,
        (unsigned int)numX, (unsigned int)numY, (unsigned int)numZ);
  }
#endif
}

#ifdef METAL_COMPUTE
template <typename T1>
GdftR2<T1>::~GdftR2() {
  if constexpr (std::is_same<T1, float>::value) {
    if (metalCtx) {
      metal_dft_destroy(static_cast<MetalDFTContext*>(metalCtx));
      metalCtx = nullptr;
    }
  }
}
#endif

template <typename T1>
void GdftR2<T1>::calcGradientMaps(Col<T1> &gx, Col<T1> &gy, Col<T1> &gz) {

  gx = Cd(FM,0);
  gy = Cd(FM,1);
  gz = Cd(FM,2);

}

template <typename T1>
Col<T1> GdftR2<T1>::Cd(const Col<T1> &d, uword dim) const {

        Col<T1> out(numX * numY * numZ);
        out.zeros();
        uword ll, jj, kk;
        switch (dim) {
        case (uword)0:
                ll = 1;
                jj = 0;
                kk = 0;
                break;
        case (uword)1:
                ll = 0;
                jj = 1;
                kk = 0;
                break;
        case (uword)2:
                ll = 0;
                jj = 0;
                kk = 1;
                break;
        default:
                std::cout << "Warning differences along dimension greater than 3! "
                        "Undefined case!"
                     << std::endl;
        }

        //Centered differences
        uword offset = ll + jj * numY + kk * numX * numY;
        for (uword ii = offset; ii < (numY * numX * numZ - offset); ii++) {
                out(ii) = d(ii + offset) - d(ii - offset);
        }

        // 'left' edge
        for (uword ii = 0; ii < offset; ii++) {
                out(ii) = d(ii + offset) - d(ii);
        }

        // 'right' edge
        for (uword ii = (numY * numX * numZ - offset); ii < (numY * numX * numZ); ii++) {
                out(ii) = d(ii) - d(ii - offset);
        }
 
 
        return out;
}



// Overloaded methods for forward and adjoint transform
// Forward transform operation
template <typename T1>
Col<complex<T1>> GdftR2<T1>::operator*(const Col<complex<T1>> &d) const {
  RANGE()
  Col<T1> realData = real(d);
  Col<T1> imagData = imag(d);

  Col<T1> realXformedData;
  Col<T1> imagXformedData;
  realXformedData.zeros(this->n1);
  imagXformedData.zeros(this->n1);

#ifdef METAL_COMPUTE
  if constexpr (std::is_same<T1, float>::value) {
    if (metalCtx) {
      metal_dft_forward(static_cast<MetalDFTContext*>(metalCtx),
          realData.memptr(), imagData.memptr(),
          realXformedData.memptr(), imagXformedData.memptr());

      Col<complex<T1>> XformedData(this->n1);
      XformedData.set_real(realXformedData);
      XformedData.set_imag(imagXformedData);
      return XformedData.eval();
    }
  }
#endif

  ftCpuWithGrads<T1>(realXformedData.memptr(), imagXformedData.memptr(),
            realData.memptr(), imagData.memptr(),
            kx.memptr(), ky.memptr(), kz.memptr(), ix.memptr(), iy.memptr(),
            iz.memptr(), FM.memptr(), Gx.memptr(), Gy.memptr(), Gz.memptr(),
            t.memptr(), this->n1, this->n2, this->numX,
            this->numY, this->numZ);

  Col<complex<T1>> XformedData(this->n1);
  XformedData.set_real(realXformedData);
  XformedData.set_imag(imagXformedData);

  return XformedData.eval();
}

// Adjoint transform operation
template <typename T1>
Col<complex<T1>> GdftR2<T1>::operator/(const Col<complex<T1>> &d) const {
  RANGE()
  Col<T1> realData = real(d);
  Col<T1> imagData = imag(d);

  Col<T1> realXformedData;
  Col<T1> imagXformedData;
  realXformedData.zeros(this->n2);
  imagXformedData.zeros(this->n2);

#ifdef METAL_COMPUTE
  if constexpr (std::is_same<T1, float>::value) {
    if (metalCtx) {
      metal_dft_adjoint(static_cast<MetalDFTContext*>(metalCtx),
          realData.memptr(), imagData.memptr(),
          realXformedData.memptr(), imagXformedData.memptr());

      Col<complex<T1>> XformedData(this->n2);
      XformedData.set_real(realXformedData);
      XformedData.set_imag(imagXformedData);
      return XformedData.eval();
    }
  }
#endif

  iftCpuWithGrads<T1>(realXformedData.memptr(), imagXformedData.memptr(),
             realData.memptr(), imagData.memptr(),
             kx.memptr(), ky.memptr(), kz.memptr(), ix.memptr(), iy.memptr(),
             iz.memptr(), FM.memptr(), Gx.memptr(), Gy.memptr(), Gz.memptr(),
             t.memptr(), this->n1, this->n2, this->numX,
             this->numY, this->numZ);

  Col<complex<T1>> XformedData(this->n2);
  XformedData.set_real(realXformedData);
  XformedData.set_imag(imagXformedData);

  return XformedData.eval();
}
// pgCol overloads — Metal GPU for float, arma fallback for double
template <typename T1>
pgCol<pgComplex<T1>> GdftR2<T1>::
operator*(const pgCol<pgComplex<T1>>& d) const
{
#ifdef METAL_COMPUTE
    if constexpr (std::is_same<T1, float>::value) {
        if (metalCtx) {
            Col<T1> realData(n2), imagData(n2);
            const T1* src = reinterpret_cast<const T1*>(d.memptr());
            T1* rp = realData.memptr();
            T1* ip = imagData.memptr();
            for (uword j = 0; j < n2; j++) {
                rp[j] = src[2*j];
                ip[j] = src[2*j+1];
            }

            Col<T1> realOut(n1), imagOut(n1);
            metal_dft_forward(static_cast<MetalDFTContext*>(metalCtx),
                rp, ip, realOut.memptr(), imagOut.memptr());

            pgCol<pgComplex<T1>> result(n1);
            T1* dst = reinterpret_cast<T1*>(result.memptr());
            const T1* ro = realOut.memptr();
            const T1* io = imagOut.memptr();
            for (uword j = 0; j < n1; j++) {
                dst[2*j]   = ro[j];
                dst[2*j+1] = io[j];
            }
            return result;
        }
    }
#endif
    Col<complex<T1>> armaResult = this->operator*(d.getArma());
    return pgCol<pgComplex<T1>>(armaResult);
}

template <typename T1>
pgCol<pgComplex<T1>> GdftR2<T1>::
operator/(const pgCol<pgComplex<T1>>& d) const
{
#ifdef METAL_COMPUTE
    if constexpr (std::is_same<T1, float>::value) {
        if (metalCtx) {
            Col<T1> realData(n1), imagData(n1);
            const T1* src = reinterpret_cast<const T1*>(d.memptr());
            T1* rp = realData.memptr();
            T1* ip = imagData.memptr();
            for (uword j = 0; j < n1; j++) {
                rp[j] = src[2*j];
                ip[j] = src[2*j+1];
            }

            Col<T1> realOut(n2), imagOut(n2);
            metal_dft_adjoint(static_cast<MetalDFTContext*>(metalCtx),
                rp, ip, realOut.memptr(), imagOut.memptr());

            pgCol<pgComplex<T1>> result(n2);
            T1* dst = reinterpret_cast<T1*>(result.memptr());
            const T1* ro = realOut.memptr();
            const T1* io = imagOut.memptr();
            for (uword j = 0; j < n2; j++) {
                dst[2*j]   = ro[j];
                dst[2*j+1] = io[j];
            }
            return result;
        }
    }
#endif
    Col<complex<T1>> armaResult = this->operator/(d.getArma());
    return pgCol<pgComplex<T1>>(armaResult);
}

// Explicit Instantiations
template class GdftR2<float>;
template class GdftR2<double>;
