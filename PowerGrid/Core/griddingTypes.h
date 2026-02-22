/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [griddingTypes.hpp]

    Synopsis    [Support times.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

/// @file griddingTypes.h
/// @brief POD structs used by the gridding kernels.

#ifndef PowerGrid_griddingTypes_h
#define PowerGrid_griddingTypes_h

/// @brief Gridding configuration parameters.
///
/// @tparam T1  Floating-point precision type.
template<typename T1>
struct parameters{
    /// @brief Total number of k-space samples.
    int numSamples;
    /// @brief Image dimensions [x, y, z].
    int imageSize[3];
    /// @brief Oversampled grid dimensions [x, y, z].
    int gridSize[3];
    /// @brief Grid oversampling factor (gridSize = gridOS x imageSize).
    T1 gridOS;
    /// @brief Kaiser-Bessel kernel width in grid cells.
    T1 kernelWidth;
    /// @brief Bin size for sorted-gridding acceleration.
    int binsize;
    /// @brief Non-zero to use precomputed LUT; zero for direct kernel evaluation.
    int useLUT;
    /// @brief Non-zero to synchronise after each OpenACC kernel.
    int sync;
};

/// @brief Single k-space sample with trajectory coordinates and density compensation.
///
/// @tparam T1  Floating-point precision type.
template <typename T1>
struct ReconstructionSample{
    /// @brief Real part of the k-space sample value.
    T1 real;
    /// @brief Imaginary part of the k-space sample value.
    T1 imag;
    /// @brief k-space x-coordinate.
    T1 kX;
    /// @brief k-space y-coordinate.
    T1 kY;
    /// @brief k-space z-coordinate.
    T1 kZ;
    /// @brief Sample density compensation weight.
    T1 sdc;
    /// @brief Per-sample readout time (s).
    T1 t;
    /// @brief Padding for alignment.
    T1 dummy;
};

#endif
