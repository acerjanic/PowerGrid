/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [directRecon.h]

    Synopsis    [Code for calculating sum-of-squares images from fully samled 
                    non-cartesian data]

    Description []

    Revision    [0.2.0; Alex Cerjanic, BIOE UIUC]

    Date        [11/10/2019]

 *****************************************************************************/

/// @file directRecon.h
/// @brief Direct sum-of-squares reconstruction for fully sampled non-Cartesian data.

#ifndef POWERGRID_DIRECTRECON_H
#define POWERGRID_DIRECTRECON_H

#include "Core/PowerGrid.h"
#include "processIsmrmrd.hpp"
using namespace arma;

/// @brief Compute 3-D density compensation weights for non-Cartesian trajectories.
///
/// @tparam T    Floating-point precision type (`float` or `double`).
/// @param kx    k-space x-coordinates, length n1.
/// @param ky    k-space y-coordinates, length n1.
/// @param kz    k-space z-coordinates, length n1.
/// @param Ninplane  Number of in-plane k-space samples per slice.
/// @param Nz    Number of slices / z-partitions.
/// @returns     Density compensation weights vector, length n1.
template<typename T>
Col<T> calc3DDensityCompensation(Col<T> kx, Col<T> ky, Col<T> kz, uword Ninplane, uword Nz);

/// @brief Grid all coil images from an ISMRMRD dataset.
///
/// Reads raw k-space data from the ISMRMRD dataset, applies density compensation,
/// and grids each coil image using the NUFFT adjoint.
///
/// @tparam T        Floating-point precision type (`float` or `double`).
/// @param Ninplane  Number of in-plane readout samples.
/// @param Nz        Number of z-partitions / slices.
/// @param d         Pointer to the open ISMRMRD dataset.
/// @param hdr       Pointer to the parsed ISMRMRD header.
/// @param acqTrack  Pointer to acquisition tracking object.
/// @param NPhase    Phase index to reconstruct.
/// @param NEcho     Echo index to reconstruct.
/// @param NAvg      Average index to reconstruct.
/// @param NRep      Repetition index to reconstruct.
/// @returns         Stacked coil images as a flat complex column vector.
template<typename T>
Col<complex<T>> gridCoilImages(uword Ninplane, uword Nz, ISMRMRD::Dataset *d, ISMRMRD::IsmrmrdHeader *hdr,
                            acqTracking *acqTrack, uword NPhase, uword NEcho, uword NAvg, uword NRep);

/// @brief Compute a sum-of-squares combination image from individual coil images.
///
/// @tparam T         Floating-point precision type (`float` or `double`).
/// @param coilImages Stacked coil images, length Nx*Nz*NSliceMax*NCoils.
/// @param Nx         In-plane image size (pixels per slice).
/// @param Nz         Number of z-partitions.
/// @param NSliceMax  Maximum number of slices.
/// @param NCoils     Number of receiver coils.
/// @returns          Sum-of-squares combined image, length Nx*Nz*NSliceMax.
template<typename T>
Col<T> calcSumOfSquaresImage(Col<complex<T>> coilImages, uword Nx, uword Nz, uword NSliceMax, uword NCoils);

extern template Col<float> calc3DDensityCompensation<float>(Col<float>, Col<float>, Col<float>, uword, uword);
extern template Col<double> calc3DDensityCompensation<double>(Col<double>, Col<double>, Col<double>, uword, uword);

extern template Col<complex<float>> gridCoilImages<float>( uword, uword, ISMRMRD::Dataset*, ISMRMRD::IsmrmrdHeader*, acqTracking*,  uword,  uword,  uword,  uword);
extern template Col<complex<double>> gridCoilImages<double>( uword, uword, ISMRMRD::Dataset*, ISMRMRD::IsmrmrdHeader*, acqTracking*,  uword,  uword,  uword,  uword);

extern template Col<float> calcSumOfSquaresImage<float>(Col<complex<float>> coilImages, uword Nx, uword Nz, uword NSliceMax, uword NCoils);
extern template Col<double> calcSumOfSquaresImage<double>(Col<complex<double>> coilImages, uword Nx, uword Nz, uword NSliceMax, uword NCoils);
#endif // DIRECTRECON

