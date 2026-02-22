/*
   (C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
   All rights reserved.

   See LICENSE.txt for the University of Illinois/NCSA Open Source license.

   Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
 */

/*****************************************************************************

    File Name   [pcSENSE.cpp]

    Synopsis    [Implements a phase corrected SENSE algorithm. ]

    Description []

    Revision    [0.1.0; Joseph Holtrop, BIOE UIUC]

    Date        [4/19/2016]

*****************************************************************************/
#include "pcSENSE.h"
#include "Core/PGLog.hpp"

template <typename T1> pcSENSE<T1>::~pcSENSE() {

	for (uword jj = 0; jj < Ns; jj++) {
		//delete G[jj];
		delete AObj[jj];
	}
	//delete[] G;
	delete[] AObj;

}

// Class constructor
template <typename T1>
pcSENSE<T1>::pcSENSE(Col<T1> kx, Col<T1> ky, Col<T1> kz, uword nx, uword ny,
                     uword nz, uword nc, Col<T1> t, Col<complex<T1>> SENSEmap,
                     Col<T1> FieldMap, Col<T1> ShotPhaseMap) {
        Ni = nx * ny * nz;
        Nc = nc;
        Ns = ShotPhaseMap.n_elem / Ni;
        Nd = kx.n_elem / Ns;
        PG_INFO("pcSENSE: Nd={}, Ns={}, Nc={}, Ni={}", Nd, Ns, Nc, Ni);
        SMap = reshape(SENSEmap, Ni, Nc);
        PMap = reshape(ShotPhaseMap, Ni, Ns);
        FMap = FieldMap;
        Kx = reshape(kx, Nd, Ns);
        Ky = reshape(ky, Nd, Ns);
        Kz = reshape(kz, Nd, Ns);
        Nx = nx;
        Ny = ny;
        Nz = nz;
        Tvec = reshape(t, Nd, Ns);

        Cube<T1> ix;
        ix.zeros(Nx, Ny, Nz);
        Cube<T1> iy;
        iy.zeros(Nx, Ny, Nz);
        Cube<T1> iz;
        iz.zeros(Nx, Ny, Nz);

        // generate the image space coordinates of the voxels we want to reconstruct
        // after vectorizing ix and iy the image coordinates must match the Field and
        // SENSe map image coordinates
        for (uword ii = 0; ii < Ny; ii++) { // y
                for (uword jj = 0; jj < Nx; jj++) { // x
                        for (uword kk = 0; kk < Nz; kk++) { // z
                                ix(ii, jj, kk) = ((T1)jj - (T1)Nx / 2.0) / ((T1)Nx);
                                iy(ii, jj, kk) = ((T1)ii - (T1)Ny / 2.0) / ((T1)Ny);
                                iz(ii, jj, kk) = ((T1)kk - (T1)Nz / 2.0) / ((T1)Nz);
                        }
                }
        }

        Ix = vectorise(ix);
        Iy = vectorise(iy);
        Iz = vectorise(iz);

        AObj = new Gdft<T1> *[Ns];

        // Initialize the field correction and G objects we need for this
        // reconstruction
        for (uword jj = 0; jj < Ns; jj++) {

                AObj[jj] =
                        new Gdft<T1>(Nd, Nx * Ny * Nz, Kx.col(jj), Ky.col(jj), Kz.col(jj), Ix,
                                     Iy, Iz, vectorise(FMap), vectorise(Tvec.col(jj)));
        }
        
        //Precompute some things used in the forward and adjoint operations
        expiPMap = conj(exp(-i * PMap));       // exp(+i*PMap) — forward weight
        conjExpiPMap = conj(expiPMap);          // exp(-i*PMap) — adjoint weight
        conjSMap = conj(SMap);

#ifdef METAL_COMPUTE
        if constexpr (std::is_same<T1, float>::value) {
            SMap_pg = pgMat<pgComplex<T1>>(SMap);
            conjSMap_pg = pgMat<pgComplex<T1>>(conjSMap);
            expiPMap_pg = pgMat<pgComplex<T1>>(expiPMap);
            conjExpiPMap_pg = pgMat<pgComplex<T1>>(conjExpiPMap);
        }
#endif
}

// Overloaded operators go here

// Forward transformation is *
// d is the vector of data of type T1, note it is const, so we don't modify it
// directly rather return another vector of type T1
template <typename T1>
Col<complex<T1> > pcSENSE<T1>::operator*(const Col<complex<T1> > &d) const {
        RANGE("pcSENSE::operator*")

#ifdef METAL_COMPUTE
        if constexpr (std::is_same<T1, float>::value) {
            // Metal path: pgCol/pgMat with GPU-dispatched element-wise ops
            pgCol<pgComplex<T1>> d_pg(d);
            pgMat<pgComplex<T1>> outData_pg(Nd, Ns * Nc);

            for (unsigned int ii = 0; ii < Nc; ii++) {
                for (unsigned int jj = 0; jj < Ns; jj++) {
                    // Compute weight = SMap(:,ii) .* expiPMap(:,jj)
                    pgCol<pgComplex<T1>> weight = SMap_pg.col_copy(ii);
                    weight %= expiPMap_pg.col(jj);

                    // weighted = d .* weight
                    pgCol<pgComplex<T1>> weighted(d_pg);
                    weighted %= weight;

                    // Forward transform (arma boundary)
                    Col<complex<T1>> result = (*AObj[jj]) * weighted.getArma();
                    outData_pg.set_col(jj + ii * Ns, pgCol<pgComplex<T1>>(result));
                }
            }

            pgCol<pgComplex<T1>> outVec = vectorise(outData_pg);
            return outVec.getArma();
        }
#endif

        // Armadillo path (double, or non-Metal builds)
        Mat<complex<T1> > outData = zeros<Mat<complex<T1> > >(Nd, Ns * Nc);
        for (unsigned int ii = 0; ii < Nc; ii++) {
                for (unsigned int jj = 0; jj < Ns; jj++) {
                        outData.col(jj + ii * Ns) =
                                (*AObj[jj]) * (d % (SMap.col(ii) % expiPMap.col(jj)));
                }
        }
        return vectorise(outData);
}

// For the adjoint operation, we have to weight the adjoint transform of the
// coil data by the SENSE map.
template <typename T1>
Col<complex<T1> > pcSENSE<T1>::operator/(const Col<complex<T1> > &d) const {
        RANGE("pcSENSE::operator/");

#ifdef METAL_COMPUTE
        if constexpr (std::is_same<T1, float>::value) {
            // Metal path: pgCol/pgMat with GPU-dispatched element-wise ops
            // Reshape input data into (Nd x Ns*Nc) matrix
            pgCol<pgComplex<T1>> d_pg(d);
            pgMat<pgComplex<T1>> inData_pg(d_pg, Nd, Ns * Nc);

            pgCol<pgComplex<T1>> outData_pg(Ni);
            outData_pg.zeros();

            for (unsigned int ii = 0; ii < Nc; ii++) {
                for (unsigned int jj = 0; jj < Ns; jj++) {
                    // Compute weight = conj(SMap(:,ii)) .* conj(expiPMap(:,jj))
                    pgCol<pgComplex<T1>> weight = conjSMap_pg.col_copy(ii);
                    weight %= conjExpiPMap_pg.col(jj);

                    // Adjoint transform (arma boundary)
                    pgCol<pgComplex<T1>> seg = inData_pg.col_copy(jj + ii * Ns);
                    Col<complex<T1>> adjResult = (*AObj[jj]) / seg.getArma();
                    pgCol<pgComplex<T1>> adjResult_pg(adjResult);

                    // Weight by sensitivity × phase and accumulate
                    adjResult_pg %= weight;
                    outData_pg += adjResult_pg;
                }
            }

            return outData_pg.getArma();
        }
#endif

        // Armadillo path (double, or non-Metal builds)
        Mat<complex<T1> > inData = reshape(d, Nd, Ns * Nc);
        Col<complex<T1> > outData = zeros<Col<complex<T1> > >(Ni);
        for (unsigned int ii = 0; ii < Nc; ii++) {
                for (unsigned int jj = 0; jj < Ns; jj++) {
                        outData += (conjSMap.col(ii) % conjExpiPMap.col(jj)) %
                                   ((*AObj[jj]) / inData.col(jj + ii * Ns));
                }
        }
        return vectorise(outData);
}

// pgCol forward: image → k-space through coil/phase sensitivities
template <typename T1>
pgCol<pgComplex<T1>> pcSENSE<T1>::
operator*(const pgCol<pgComplex<T1>>& d) const
{
#ifdef METAL_COMPUTE
    if constexpr (std::is_same<T1, float>::value) {
        pgMat<pgComplex<T1>> outData_pg(Nd, Ns * Nc);

        for (unsigned int ii = 0; ii < Nc; ii++) {
            for (unsigned int jj = 0; jj < Ns; jj++) {
                pgCol<pgComplex<T1>> weight = SMap_pg.col_copy(ii);
                weight %= expiPMap_pg.col(jj);

                pgCol<pgComplex<T1>> weighted(d);
                weighted %= weight;

                outData_pg.set_col(jj + ii * Ns, (*AObj[jj]) * weighted);
            }
        }

        return vectorise(outData_pg);
    }
#endif
    Col<complex<T1>> armaResult = this->operator*(d.getArma());
    return pgCol<pgComplex<T1>>(armaResult);
}

// pgCol adjoint: k-space → image through conjugate coil/phase sensitivities
template <typename T1>
pgCol<pgComplex<T1>> pcSENSE<T1>::
operator/(const pgCol<pgComplex<T1>>& d) const
{
#ifdef METAL_COMPUTE
    if constexpr (std::is_same<T1, float>::value) {
        pgMat<pgComplex<T1>> inData_pg(d, Nd, Ns * Nc);

        pgCol<pgComplex<T1>> outData_pg(Ni);
        outData_pg.zeros();

        for (unsigned int ii = 0; ii < Nc; ii++) {
            for (unsigned int jj = 0; jj < Ns; jj++) {
                pgCol<pgComplex<T1>> weight = conjSMap_pg.col_copy(ii);
                weight %= conjExpiPMap_pg.col(jj);

                pgCol<pgComplex<T1>> seg = inData_pg.col_copy(jj + ii * Ns);
                pgCol<pgComplex<T1>> adjResult = (*AObj[jj]) / seg;

                adjResult %= weight;
                outData_pg += adjResult;
            }
        }

        return outData_pg;
    }
#endif
    Col<complex<T1>> armaResult = this->operator/(d.getArma());
    return pgCol<pgComplex<T1>>(armaResult);
}

// Explicit Instantiation
template class pcSENSE<float>;
template class pcSENSE<double>;
