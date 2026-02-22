/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [SENSE.cpp]

    Synopsis    [Object implementing sensitivity encoding reconstructions. ]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/
#include "SENSE.h"
#include "Core/PGLog.hpp"

using namespace arma;

// We are using two template types at the moment. One for the type of data to be
// processed (ie Col<cx_double>) and one for the type of G object (ie
// Gfft<Col<cx_double>>
template <typename T1, typename Tobj>
SENSE<T1, Tobj>::SENSE(Tobj& G, Col<complex<T1>> SENSEmap, uword a, uword b,
    uword c)
{
    n1 = a;
    n2 = b;
    nc = c;
    PG_INFO("SENSE: n1={}, n2={}, nc={}", n1, n2, nc);
    G_obj = &G;
    SMap = reshape(SENSEmap, n2, nc);
    conjSMap = conj(SMap);
    outData.set_size(n1, nc);
    outImg.set_size(n2, nc);

#ifdef METAL_COMPUTE
    if constexpr (std::is_same<T1, float>::value) {
        SMap_pg = pgMat<pgComplex<T1>>(SMap);
        conjSMap_pg = pgMat<pgComplex<T1>>(conjSMap);
        outData_pg = pgMat<pgComplex<T1>>(n1, nc);
        outImg_pg = pgMat<pgComplex<T1>>(n2, nc);
    }
#endif
}

// Overloaded operators go here

// Forward transformation is *
// d is the vector of data of type T1, note it is const, so we don't modify it
// directly rather return another vector of type T1
template <typename T1, typename Tobj>
Col<complex<T1>> SENSE<T1, Tobj>::operator*(const Col<complex<T1>>& d) const
{
    RANGE("SENSE::operator*")

#ifdef METAL_COMPUTE
    if constexpr (std::is_same<T1, float>::value) {
        // Metal path: pgCol/pgMat with GPU-dispatched element-wise ops
        pgCol<pgComplex<T1>> d_pg(d);

        for (unsigned int ii = 0; ii < this->nc; ii++) {
            // Element-wise multiply: weighted = d .* SMap(:,ii)
            // SMap_pg.col(ii) returns a view — Metal cvec_mul dispatches here
            pgCol<pgComplex<T1>> weighted(d_pg);
            weighted %= SMap_pg.col(ii);

            // Forward transform (pgCol overload — no arma conversion)
            outData_pg.set_col(ii, (*this->G_obj) * weighted);
        }

        pgCol<pgComplex<T1>> outVec = vectorise(outData_pg);
        return outVec.getArma();
    }
#endif

    // Armadillo path (double, or non-Metal builds)
    outImg = this->SMap;
#pragma omp parallel for schedule(dynamic) shared(outData, d, SMap)
    for (unsigned int ii = 0; ii < this->nc; ii++) {
        outImg.unsafe_col(ii) %= d;
    }

    for (unsigned int ii = 0; ii < this->nc; ii++) {
        outData.unsafe_col(ii) = (*this->G_obj) * outImg.unsafe_col(ii);
    }

    return vectorise(outData);
}

// For the adjoint operation, we have to weight the adjoint transform of the
// coil data by the SENSE map.
template <typename T1, typename Tobj>
Col<complex<T1>> SENSE<T1, Tobj>::operator/(const Col<complex<T1>>& d) const
{
    RANGE("SENSE::operator/")

#ifdef METAL_COMPUTE
    if constexpr (std::is_same<T1, float>::value) {
        // Metal path: pgCol/pgMat with GPU-dispatched element-wise ops
        pgCol<pgComplex<T1>> d_pg(d);

        for (unsigned int ii = 0; ii < this->nc; ii++) {
            // Extract coil's k-space data slice (view, no copy)
            pgCol<pgComplex<T1>> dSlice = d_pg.subvec(ii * n1, (ii + 1) * n1 - 1);

            // Adjoint transform (pgCol overload — no arma conversion)
            outImg_pg.set_col(ii, (*this->G_obj) / dSlice);
        }

        // Weight by conjugate SENSE map: outImg(:,ii) .*= conj(SMap(:,ii))
        // Metal cvec_mul dispatches via the view's operator%=
        for (unsigned int ii = 0; ii < this->nc; ii++) {
            outImg_pg.col(ii) %= conjSMap_pg.col(ii);
        }

        // Sum across coils (row-wise sum, dim=1)
        pgCol<pgComplex<T1>> sumVec = sum(outImg_pg, 1);
        return sumVec.getArma();
    }
#endif

    // Armadillo path (double, or non-Metal builds)
    for (unsigned int ii = 0; ii < this->nc; ii++) {
        outImg.unsafe_col(ii) = (*this->G_obj) / d.subvec((ii)*n1, ((ii + 1) * n1) - 1);
    }

#pragma omp parallel for schedule(dynamic) shared(outImg, conjSMap)
    for (unsigned int ii = 0; ii < this->nc; ii++) {
        outImg.unsafe_col(ii) %= this->conjSMap.unsafe_col(ii);
    }

    return sum(outImg, 1);
}

// pgCol forward: image → k-space through coil sensitivities
template <typename T1, typename Tobj>
pgCol<pgComplex<T1>> SENSE<T1, Tobj>::
operator*(const pgCol<pgComplex<T1>>& d) const
{
#ifdef METAL_COMPUTE
    if constexpr (std::is_same<T1, float>::value) {
        for (unsigned int ii = 0; ii < this->nc; ii++) {
            pgCol<pgComplex<T1>> weighted(d);
            weighted %= SMap_pg.col(ii);
            outData_pg.set_col(ii, (*this->G_obj) * weighted);
        }
        return vectorise(outData_pg);
    }
#endif
    Col<CxT1> armaResult = this->operator*(d.getArma());
    return pgCol<pgComplex<T1>>(armaResult);
}

// pgCol adjoint: k-space → image through conjugate coil sensitivities
template <typename T1, typename Tobj>
pgCol<pgComplex<T1>> SENSE<T1, Tobj>::
operator/(const pgCol<pgComplex<T1>>& d) const
{
#ifdef METAL_COMPUTE
    if constexpr (std::is_same<T1, float>::value) {
        for (unsigned int ii = 0; ii < this->nc; ii++) {
            pgCol<pgComplex<T1>> dSlice = d.subvec(ii * n1, (ii + 1) * n1 - 1);
            outImg_pg.set_col(ii, (*this->G_obj) / dSlice);
        }

        for (unsigned int ii = 0; ii < this->nc; ii++) {
            outImg_pg.col(ii) %= conjSMap_pg.col(ii);
        }

        return sum(outImg_pg, 1);
    }
#endif
    Col<CxT1> armaResult = this->operator/(d.getArma());
    return pgCol<pgComplex<T1>>(armaResult);
}

// Explicit Instantiations
template class SENSE<float, Gnufft<float>>;
template class SENSE<float, TimeSegmentation<float, Gnufft<float>>>;
template class SENSE<double, Gnufft<double>>;
template class SENSE<double, TimeSegmentation<double, Gnufft<double>>>;
template class SENSE<float, Gdft<float>>;
template class SENSE<double, Gdft<double>>;
template class SENSE<float, GdftR2<float>>;
template class SENSE<double, GdftR2<double>>;
