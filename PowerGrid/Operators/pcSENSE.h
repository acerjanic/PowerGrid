/*
   (C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
   All rights reserved.

   See LICENSE.txt for the University of Illinois/NCSA Open Source license.

   Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
 */

/*****************************************************************************

    File Name   [pcSENSE.h]

    Synopsis    [Implements a phase corrected SENSE algorithm. ]

    Description []

    Revision    [0.1.0; Joseph Holtrop, BIOE UIUC]

    Date        [4/19/2016]

*****************************************************************************/

/// @file pcSENSE.h
/// @brief Phase-corrected SENSE operator for multi-shot non-Cartesian acquisitions.

#ifndef PowerGrid_pcSENSE_hpp
#define PowerGrid_pcSENSE_hpp

#include "Core/PGIncludes.h"
#include "Gdft.h"
#include "Gnufft.h"
#include "Gridding/TimeSegmentation.h"

#include "Core/pgCol.hpp"
#include "Core/pgMat.hpp"

using namespace arma;
//using namespace PowerGrid;

/// @brief Phase-corrected multi-shot SENSE operator.
///
/// Extends the standard SENSE model with per-shot phase maps to correct for
/// inter-shot phase inconsistencies in multi-shot non-Cartesian acquisitions.
/// Uses a per-shot Gdft encoding operator to handle arbitrary trajectories.
///
/// @tparam T1  Floating-point precision type (`float` or `double`).
template <typename T1> class pcSENSE {
typedef std::complex<T1> CxT1;

public:
/// @brief Default constructor.
pcSENSE();

~pcSENSE();

/// @brief Number of k-space samples per shot.
uword Nd = 0;
/// @brief Number of image pixels.
uword Ni = 0;
/// @brief Number of receiver coils.
uword Nc = 0;
/// @brief Number of shots.
uword Ns = 0;
/// @brief Coil sensitivity maps, size Ni x Nc.
Mat<CxT1> SMap;
/// @brief Conjugate sensitivity maps, size Ni x Nc.
Mat<CxT1> conjSMap;
/// @brief Per-shot phase maps in radians, size Ni x Ns.
Mat<T1> PMap;
/// @brief Complex exponential of phase maps (exp(i*PMap)), size Ni x Ns.
Mat<CxT1> expiPMap;
/// @brief Conjugate complex exponential of phase maps (exp(-i*PMap)), size Ni x Ns.
Mat<CxT1> conjExpiPMap;
/// @brief Off-resonance field map (rad/s), length Ni.
Col<T1> FMap;
/// @brief k-space x-coordinates matrix, size (samples/shot) x Ns.
Mat<T1> Kx;
/// @brief k-space y-coordinates matrix, size (samples/shot) x Ns.
Mat<T1> Ky;
/// @brief k-space z-coordinates matrix, size (samples/shot) x Ns.
Mat<T1> Kz;
/// @brief Readout time vector (s), size (samples/shot) x Ns.
Mat<T1> Tvec;
/// @brief Image size in x.
uword Nx;
/// @brief Image size in y.
uword Ny;
/// @brief Image size in z.
uword Nz;
/// @brief Image-space x-coordinates, length Ni.
Col<T1> Ix;
/// @brief Image-space y-coordinates, length Ni.
Col<T1> Iy;
/// @brief Image-space z-coordinates, length Ni.
Col<T1> Iz;
CxT1 i = CxT1(0., 1.);
/// @brief Interpolation type: 1 = Hanning, 2 = exact min-max.
uword type = 1;
/// @brief Number of time segments for field-map correction.
uword L = 20;
/// @brief Array of per-shot encoding operator pointers.
Gdft<T1> **AObj = NULL;

#ifdef METAL_COMPUTE
// pgMat copies of sensitivity/phase maps for Metal GPU dispatch (float only).
pgMat<pgComplex<T1>> SMap_pg;
pgMat<pgComplex<T1>> conjSMap_pg;
pgMat<pgComplex<T1>> expiPMap_pg;
pgMat<pgComplex<T1>> conjExpiPMap_pg;
#endif
	//TimeSegmentation <T1, Gnufft<T1>> **AObj = NULL;

/// @brief Construct a phase-corrected SENSE operator.
///
/// @param kx          k-space x-coordinates (all shots concatenated), length Nd*Ns.
/// @param ky          k-space y-coordinates, length Nd*Ns.
/// @param kz          k-space z-coordinates, length Nd*Ns.
/// @param nx          Image size in x.
/// @param ny          Image size in y.
/// @param nz          Image size in z.
/// @param nc          Number of receiver coils.
/// @param t           Per-sample readout time vector (s), length Nd*Ns.
/// @param SENSEmap    Sensitivity maps as flat vector, length Ni*Nc.
/// @param FieldMap    Off-resonance field map (rad/s), length Ni.
/// @param ShotPhaseMap  Per-shot phase maps in radians, length Ni*Ns.
pcSENSE(Col<T1> kx, Col<T1> ky, Col<T1> kz, uword nx, uword ny, uword nz,
        uword nc, Col<T1> t, Col<CxT1> SENSEmap, Col<T1> FieldMap,
        Col<T1> ShotPhaseMap);

/// @brief Forward phase-corrected SENSE transform: image -> k-space.
///
/// @param d  Input image vector of length Ni.
/// @returns  Stacked k-space vector of length Nd*Ns*Nc.
Col<CxT1> operator*(const Col<CxT1> &d) const;

/// @brief Adjoint phase-corrected SENSE transform: k-space -> image.
///
/// @param d  Input k-space vector of length Nd*Ns*Nc.
/// @returns  Output image vector of length Ni.
Col<CxT1> operator/(const Col<CxT1> &d) const;

/// @brief Forward phase-corrected SENSE transform (pgCol overload for Metal path).
/// @param d  Input image vector of length Ni.
/// @returns  Stacked k-space vector of length Nd*Ns*Nc.
pgCol<pgComplex<T1>> operator*(const pgCol<pgComplex<T1>> &d) const;
/// @brief Adjoint phase-corrected SENSE transform (pgCol overload for Metal path).
/// @param d  Input k-space vector of length Nd*Ns*Nc.
/// @returns  Output image vector of length Ni.
pgCol<pgComplex<T1>> operator/(const pgCol<pgComplex<T1>> &d) const;
};

// Explicit Instantiation
extern template class pcSENSE<float>;
extern template class pcSENSE<double>;
#endif
