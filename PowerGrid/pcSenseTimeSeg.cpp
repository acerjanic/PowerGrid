/*
   (C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
   All rights reserved.

   See LICENSE.txt for the University of Illinois/NCSA Open Source license.

   Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
 */

/*****************************************************************************

    File Name   [pcSenseTimeSeg.cpp]

    Synopsis    [Implements a phase corrected SENSE algorithm. ]

    Description []

    Revision    [0.1.0; Joseph Holtrop, BIOE UIUC]

    Date        [4/19/2016]

*****************************************************************************/
#include "pcSenseTimeSeg.h"

template <typename T1>
pcSenseTimeSeg<T1>::~pcSenseTimeSeg()
{

    for (uword jj = 0; jj < Ns; jj++) {
        delete G[jj];
        delete AObj[jj];
    }
    delete[] G;
    delete[] AObj;
}

// Class constructor
template <typename T1>
pcSenseTimeSeg<T1>::pcSenseTimeSeg(Col<T1> kx, Col<T1> ky, Col<T1> kz, uword nx, uword ny,
    uword nz, uword nc, Col<T1> t, uword lSegs, uword intType, Col<complex<T1> > SENSEmap,
    Col<T1> FieldMap, Col<T1> ShotPhaseMap)
{
    Ni = nx * ny * nz;
    Nc = nc;
    Ns = ShotPhaseMap.n_elem / Ni;
    Nd = kx.n_elem / Ns;
    std::cout << "Nd = " << Nd << std::endl;
    std::cout << "Ns = " << Ns << std::endl;
    std::cout << "Nc = " << Nc << std::endl;
    std::cout << "Ni = " << Ni << std::endl;
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
    L = lSegs;
    type = intType;
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

    //AObj = new Gdft<T1> *[Ns];
    G = new Gnufft<T1>*[Ns];
    AObj = new TimeSegmentation<T1, Gnufft<T1> >*[Ns];

    // Initialize the field correction and G objects we need for this
    // reconstruction
    for (uword jj = 0; jj < Ns; jj++) {

        G[jj] = new Gnufft<T1>(Nd, 2.0, Nx, Ny, Nz, Kx.col(jj), Ky.col(jj), Kz.col(jj), Ix, Iy, Iz);
        AObj[jj] = new TimeSegmentation<T1, Gnufft<T1> >(*G[jj], vectorise(FMap),
            vectorise(Tvec.col(jj)), (uword)Nd, (uword)(Nx * Ny * Nz), (uword)L, type, 1);
    }

    // Precalculate some quantities that were generated on the fly to save time and check for errors
    shotSpecificSenseMap.zeros(Ni, Nc * Ns);
    conjShotSpecificSenseMap.zeros(Ni, Nc * Ns);
    // Precomputing some handy intermediate values
    Mat<CxT1> conjSMap = conj(SMap);
    Mat<CxT1> expiPMap = exp(-i * PMap);
    Mat<CxT1> conjExpiPMap = conj(exp(-i * PMap));

    for (unsigned int ii = 0; ii < Nc; ii++) {
        // Shot loop. Each shot has its own kspace trajectory
        for (unsigned int jj = 0; jj < Ns; jj++) {
            shotSpecificSenseMap.col(jj + ii * Ns) = SMap.col(ii) % expiPMap.col(jj);
            conjShotSpecificSenseMap.col(jj + ii * Ns) = conjSMap.col(ii) % conjExpiPMap.col(jj);
        }
    }

    //Check for NaN in the precomputed shot specific sense maps;
    if(shotSpecificSenseMap.has_nan())
        std::cout << "WARNING : shotSpecificSenseMap has NAN!! " << std::endl;

    if(conjShotSpecificSenseMap.has_nan())
        std::cout << "WARNING : conjShotSpecificSenseMap has NAN!! " << std::endl;

#ifdef METAL_COMPUTE
    if constexpr (std::is_same<T1, float>::value) {
        shotSpecificSenseMap_pg = pgMat<pgComplex<T1>>(shotSpecificSenseMap);
        conjShotSpecificSenseMap_pg = pgMat<pgComplex<T1>>(conjShotSpecificSenseMap);
    }
#endif
}

// Overloaded operators go here

// Forward transformation is *
template <typename T1>
Col<complex<T1> > pcSenseTimeSeg<T1>::operator*(const Col<complex<T1> >& d) const
{
    RANGE("pcSenseTimeSeg::operator*")
    Mat<complex<T1> > outData = zeros<Mat<complex<T1> > >(Nd, Ns * Nc);
    // Coil loop. Each coil exists for each shot, so we need to work with these.
    for (unsigned int ii = 0; ii < Nc; ii++) {

        // Shot loop. Each shot has its own kspace trajectory
        for (unsigned int jj = 0; jj < Ns; jj++) {
            outData.col(jj + ii * Ns) = (*AObj[jj]) * (d % shotSpecificSenseMap.col(jj + ii * Ns));
        }
    }
    
    if (outData.has_nan())
        std::cout << "Warning:: Output of operator* in pcSenseTimeSeg is about to return NaN" << std::endl;
    // equivalent to returning col(output) in MATLAB with IRT
    return vectorise(outData);
}

// For the adjoint operation, we have to weight the adjoint transform of the
// coil data by the SENSE map.
template <typename T1>
Col<complex<T1> > pcSenseTimeSeg<T1>::operator/(const Col<complex<T1> >& d) const
{
    RANGE("pcSenseTimeSeg::operator/");
    Mat<complex<T1> > inData = reshape(d, Nd, Ns * Nc);
    Col<complex<T1> > outData = zeros<Col<complex<T1> > >(Ni);
    // Coil Loop - for each shot we have a full set of coil data.
    for (unsigned int ii = 0; ii < Nc; ii++) {

        // Shot Loop. Each shot has it's own k-space trajectory
        for (unsigned int jj = 0; jj < Ns; jj++) {
            outData += (conjShotSpecificSenseMap.col(jj + ii * Ns)) % ((*AObj[jj]) / inData.col(jj + ii * Ns));
        }
    }

    // equivalent to returning col(output) in MATLAB with IRT
    if (outData.has_nan())
        std::cout << "Warning:: Output of operator/ in pcSenseTimeSeg is about to return NaN" << std::endl;
    return vectorise(outData);
}

// pgCol forward: image → k-space through coil/phase sensitivities + time segmentation
template <typename T1>
pgCol<pgComplex<T1>> pcSenseTimeSeg<T1>::
operator*(const pgCol<pgComplex<T1>>& d) const
{
#ifdef METAL_COMPUTE
    if constexpr (std::is_same<T1, float>::value) {
        pgMat<pgComplex<T1>> outData_pg(Nd, Ns * Nc);

        for (unsigned int ii = 0; ii < Nc; ii++) {
            for (unsigned int jj = 0; jj < Ns; jj++) {
                pgCol<pgComplex<T1>> weighted(d);
                weighted %= shotSpecificSenseMap_pg.col(jj + ii * Ns);

                outData_pg.set_col(jj + ii * Ns, (*AObj[jj]) * weighted);
            }
        }

        return vectorise(outData_pg);
    }
#endif
    Col<complex<T1>> armaResult = this->operator*(d.getArma());
    return pgCol<pgComplex<T1>>(armaResult);
}

// pgCol adjoint: k-space → image through conjugate coil/phase sensitivities + time segmentation
template <typename T1>
pgCol<pgComplex<T1>> pcSenseTimeSeg<T1>::
operator/(const pgCol<pgComplex<T1>>& d) const
{
#ifdef METAL_COMPUTE
    if constexpr (std::is_same<T1, float>::value) {
        pgMat<pgComplex<T1>> inData_pg(d, Nd, Ns * Nc);

        pgCol<pgComplex<T1>> outData_pg(Ni);
        outData_pg.zeros();

        for (unsigned int ii = 0; ii < Nc; ii++) {
            for (unsigned int jj = 0; jj < Ns; jj++) {
                pgCol<pgComplex<T1>> seg = inData_pg.col_copy(jj + ii * Ns);
                pgCol<pgComplex<T1>> adjResult = (*AObj[jj]) / seg;

                adjResult %= conjShotSpecificSenseMap_pg.col(jj + ii * Ns);
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
template class pcSenseTimeSeg<float>;
template class pcSenseTimeSeg<double>;
