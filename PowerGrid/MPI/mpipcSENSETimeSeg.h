/*
   (C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
   All rights reserved.

   See LICENSE.txt for the University of Illinois/NCSA Open Source license.

   Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
 */

/*****************************************************************************

    File Name   [mpipcSENSETimeSeg.h]

    Synopsis    [Implements a phase corrected SENSE algorithm. ]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

*****************************************************************************/

#ifndef PowerGrid_mpipcSENSETimeSeg_hpp
#define PowerGrid_mpipcSENSETimeSeg_hpp

#include "../../Support/ArmaExtensions/arma_extensions.h"
#include "../PGIncludes.h"
#include "../Gdft.h"
#include <boost/mpi.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/environment.hpp>
#include "../PowerGrid.h"

using namespace arma;

namespace bmpi = boost::mpi;
template <typename T1> class mpipcSENSETimeSeg {
typedef std::complex<T1> CxT1;

public:
mpipcSENSETimeSeg();

~mpipcSENSETimeSeg();

// Class variables go here
uword Nd = 0;   // Data size  (the size of one gdft or Gnufft object, length of
                // a single shot)
uword Ni = 0;   // Image size
uword Nc = 0;   // number of coils
uword Ns = 0;   // number of shots
Mat<CxT1> SMap;   // coil sensitivity, dimensions Image size(n1) by number of
                  // coils (nc)
Mat<T1> PMap;   // shot phase, dimensions Image size(n1) by number of shots. in
                // radians.
Col<T1> FMap;   // Fieldmap
Mat<T1> Kx;     // kspace coordinates in x direction
Mat<T1> Ky;     // kspace coordinates in y direction
Mat<T1> Kz;     // kspace coordinates in z direction
Mat<T1> Tvec;   // timing vector for a single shot (all shots assumed to have
                // same timing vector)
uword Nx;
uword Ny;
uword Nz;
Col<T1> Ix;
Col<T1> Iy;
Col<T1> Iz;
CxT1 i = CxT1(0., 1.);
uword type;   // 2 for min max time seg and 1 for Hanning
uword L;
// Gdft <T1> **AObj = NULL;
Gnufft<T1> **GObj = NULL;
TimeSegmentation <T1, Gnufft<T1>> **AObj = NULL;
// MPI Stuff
bmpi::environment *pEnv;
bmpi::communicator *pWorld;
Col<uword> shotList;
Col<uword> coilList;
std::vector<std::vector<uword> > *taskList;
// Class constructor
mpipcSENSETimeSeg(Col<T1> kx, Col<T1> ky, Col<T1> kz, uword nx, uword ny, uword nz,
           uword nc, Col<T1> t, uword L, uword type, Col<CxT1> SENSEmap, Col<T1> FieldMap,
           Col<T1> ShotPhaseMap, bmpi::environment &en,
           bmpi::communicator &wor);
// Overloaded operators go here

// Forward transformation is *
// d is the vector of data of type T1, note it is const, so we don't modify it
// directly rather return another vector of type T1
Col<CxT1> operator*(const Col<CxT1> &d) const;

// For the adjoint operation, we have to weight the adjoint transform of the
// coil data by the SENSE map.
Col<CxT1> operator/(const Col<CxT1> &d) const;
};

// Explict Instantiation
extern template class mpipcSENSETimeSeg<float>;
extern template class mpipcSENSETimeSeg<double>;

#ifdef PowerGridMPI

extern template  Col<complex<float>> reconSolve(Col<complex<float>>, mpipcSENSETimeSeg<float>&,
                               QuadPenalty<float>, Col<float>, Col<float>,
                               Col<float>, uword, uword, uword, Col<float>,
                               uword);
#endif

#endif // PowerGrid_mpipcSENSETimeSeg_hpp
