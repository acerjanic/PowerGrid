/*
   (C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
   All rights reserved.

   See LICENSE.txt for the University of Illinois/NCSA Open Source license.

   Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
 */

/*****************************************************************************

    File Name   [pcSenseTimeSeg.h]

    Synopsis    [Implements a phase corrected SENSE algorithm. ]

    Description []

    Revision    [0.1.0; Joseph Holtrop, BIOE UIUC]

    Date        [4/19/2016]

*****************************************************************************/

/// @file pcSenseTimeSeg.h
/// @brief Phase-corrected SENSE with time-segmented off-resonance correction.

#ifndef PowerGrid_pcSenseTimeSeg_hpp
#define PowerGrid_pcSenseTimeSeg_hpp

#include "Gdft.h"
#include "Gnufft.h"
#include "Core/PGIncludes.h"
#include "Gridding/TimeSegmentation.h"
#include "Core/pgCol.hpp"
#include "Core/pgMat.hpp"

using namespace arma;
//using namespace PowerGrid;

/// @brief Phase-corrected multi-shot SENSE with time-segmented field-map correction.
///
/// Combines per-shot phase correction (pcSENSE) with time-segmented off-resonance
/// correction (TimeSegmentation<T1, Gnufft<T1>>), providing accurate forward and
/// adjoint models for multi-shot non-Cartesian acquisitions with field inhomogeneity.
///
/// @tparam T1  Floating-point precision type (`float` or `double`).
template <typename T1>
class pcSenseTimeSeg {
    typedef std::complex<T1> CxT1;

public:
    /// @brief Default constructor.
    pcSenseTimeSeg();

    ~pcSenseTimeSeg();

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
    /// @brief Per-shot phase maps in radians, size Ni x Ns.
    Mat<T1> PMap;
    /// @brief Off-resonance field map (rad/s), length Ni.
    Col<T1> FMap;
    /// @brief k-space x-coordinates, size (samples/shot) x Ns.
    Mat<T1> Kx;
    /// @brief k-space y-coordinates, size (samples/shot) x Ns.
    Mat<T1> Ky;
    /// @brief k-space z-coordinates, size (samples/shot) x Ns.
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
    /// @brief Interpolation type: 1 = Hanning, 2 = min-max.
    uword type;
    /// @brief Number of time segments for field-map correction.
    uword L;

    /// @brief Shot-specific combined (sensitivity x phase) maps, size Ni x (Nc*Ns).
    Mat<CxT1> shotSpecificSenseMap;
    /// @brief Conjugate shot-specific maps, size Ni x (Nc*Ns).
    Mat<CxT1> conjShotSpecificSenseMap;

#ifdef METAL_COMPUTE
    pgMat<pgComplex<T1>> shotSpecificSenseMap_pg;
    pgMat<pgComplex<T1>> conjShotSpecificSenseMap_pg;
#endif

    /// @brief Array of per-shot NUFFT operator pointers.
    Gnufft<T1>** G = NULL;
    /// @brief Array of per-shot TimeSegmentation operator pointers.
    TimeSegmentation<T1, Gnufft<T1> >** AObj = NULL;

    /// @brief Construct a phase-corrected SENSE+time-segmentation operator.
    ///
    /// @param kx          k-space x-coordinates (all shots), length Nd*Ns.
    /// @param ky          k-space y-coordinates, length Nd*Ns.
    /// @param kz          k-space z-coordinates, length Nd*Ns.
    /// @param nx          Image size in x.
    /// @param ny          Image size in y.
    /// @param nz          Image size in z.
    /// @param nc          Number of receiver coils.
    /// @param t           Per-sample readout time (s), length Nd*Ns.
    /// @param L           Number of time segments.
    /// @param intType     Interpolation type: 1 = Hanning, 2 = min-max.
    /// @param SENSEmap    Sensitivity maps, flat vector length Ni*Nc.
    /// @param FieldMap    Off-resonance field map (rad/s), length Ni.
    /// @param ShotPhaseMap  Per-shot phase maps (rad), length Ni*Ns.
    pcSenseTimeSeg(Col<T1> kx, Col<T1> ky, Col<T1> kz, uword nx, uword ny, uword nz,
        uword nc, Col<T1> t, uword L, uword intType, Col<CxT1> SENSEmap, Col<T1> FieldMap,
        Col<T1> ShotPhaseMap);

    /// @brief Forward transform: image -> k-space with phase and field-map correction.
    ///
    /// @param d  Input image vector of length Ni.
    /// @returns  Stacked k-space vector of length Nd*Ns*Nc.
    Col<CxT1> operator*(const Col<CxT1>& d) const;

    /// @brief Adjoint transform: k-space -> image with phase and field-map correction.
    ///
    /// @param d  Input k-space vector of length Nd*Ns*Nc.
    /// @returns  Output image vector of length Ni.
    Col<CxT1> operator/(const Col<CxT1>& d) const;

    /// @brief Forward transform (pgCol overload for Metal path).
    /// @param d  Input image vector of length Ni.
    /// @returns  Stacked k-space vector of length Nd*Ns*Nc.
    pgCol<pgComplex<T1>> operator*(const pgCol<pgComplex<T1>> &d) const;
    /// @brief Adjoint transform (pgCol overload for Metal path).
    /// @param d  Input k-space vector of length Nd*Ns*Nc.
    /// @returns  Output image vector of length Ni.
    pgCol<pgComplex<T1>> operator/(const pgCol<pgComplex<T1>> &d) const;
};

// Explicit Instantiation
extern template class pcSenseTimeSeg<float>;
extern template class pcSenseTimeSeg<double>;
#endif
