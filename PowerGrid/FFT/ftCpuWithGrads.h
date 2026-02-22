/*
(C) Copyright 2010-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [ftCpuWithGrads.h]

    Synopsis    [The CPU and OpenACC annotated version of the discrete Fourier
                transform and inverse discrete Fourier transform.]

    Description []

    Revision    [0.1; Initial build; Yue Zhuo, BIOE UIUC]
    Revision    [0.1.1; Add OpenMP, Code cleaning; Xiao-Long Wu, ECE UIUC]
    Revision    [1.0a; Further optimization, Code cleaning, Adding more
 comments;
                 Xiao-Long Wu, ECE UIUC, Jiading Gai, Beckman Institute]
    Revision    [1.1; Remove OpenMP and add OpenACC annotations for GPU
                  acceleration]
    Date        [4/19/2016]

 *****************************************************************************/

/// @file ftCpuWithGrads.h
/// @brief Non-uniform DFT kernels extended with R2* field-map gradient support.

#ifndef FTWithGrads_CPU_H
#define FTWithGrads_CPU_H

/*---------------------------------------------------------------------------*/
/*  Included library headers                                                 */
/*---------------------------------------------------------------------------*/
#include "Core/PGIncludes.h"
/*---------------------------------------------------------------------------*/
/*  Namespace declared - begin                                               */
/*---------------------------------------------------------------------------*/

// namespace uiuc_mri {

/*---------------------------------------------------------------------------*/
/*  Function prototypes                                                      */
/*---------------------------------------------------------------------------*/
/*
    void
ftCpu(T1 *kdata_r, T1 *kdata_i,
      const T1 *idata_r, const T1 *idata_i,
      const DataTraj *ktraj, const DataTraj *itraj,
      const T1 *fm, const T1 *t,
      const int num_k, const int num_i
      );

    void
iftCpu(T1 *idata_r, T1 *idata_i,
       const T1 *kdata_r, const T1 *kdata_i,
       const DataTraj *ktraj, const DataTraj *itraj,
       const T1 *fm, const T1 *t,
       const int num_k, const int num_i
       );
*/
/*---------------------------------------------------------------------------*/
/*  Included library headers                                                 */
/*---------------------------------------------------------------------------*/

#include <math.h>
#include <stdio.h>
#include <string.h>
#ifdef _OPENACC
    #include "openacc.h"
#endif

//#include <tools.h>
//#include <structures.h>

/*---------------------------------------------------------------------------*/
/*  Namespace declared - begin                                               */
/*---------------------------------------------------------------------------*/

// namespace uiuc_mri {

/*---------------------------------------------------------------------------*/
/*  Function definitions                                                     */
/*---------------------------------------------------------------------------*/

/*===========================================================================*/
/*                                                                           */
/*  Synopsis    [CPU version of the sin function.]                           */
/*                                                                           */
/*  Description [This function is used to avoid additional computations when */
/*      the data values are too small.]                                      */
/*                                                                           */
/*===========================================================================*/

/// @brief Normalized sinc helper function used internally.
///
/// @tparam T1  Floating-point type.
/// @param x    Input value.
/// @returns    sinc(x) = sin(x)/x, with sinc(0) = 1.
template <typename T1> T1 sincPG(T1 x);

/// @brief Non-uniform forward DFT kernel with R2* field-map gradients (CPU / OpenACC).
///
/// Extends ftCpu by incorporating spatially-varying R2* decay via gradient
/// field maps Gx, Gy, Gz.  Computes:
///   kdata[j] = Sigma_k image[k] * exp(-i*2*pi*k*x - i*FM[k]*t[j] - R2*(k)*t[j])
/// where R2*(k) is interpolated from the gradient maps.
///
/// @tparam T1      Floating-point precision type.
/// @param kdata_r  Output k-space real part array, length num_k.
/// @param kdata_i  Output k-space imaginary part array, length num_k.
/// @param idata_r  Input image real part array, length num_i.
/// @param idata_i  Input image imaginary part array, length num_i.
/// @param kx       k-space x-coordinates, length num_k.
/// @param ky       k-space y-coordinates, length num_k.
/// @param kz       k-space z-coordinates, length num_k.
/// @param ix       Image-space x-coordinates, length num_i.
/// @param iy       Image-space y-coordinates, length num_i.
/// @param iz       Image-space z-coordinates, length num_i.
/// @param FM       Off-resonance field map (rad/s), length num_i.
/// @param Gx       R2* gradient in x, length num_x.
/// @param Gy       R2* gradient in y, length num_y.
/// @param Gz       R2* gradient in z, length num_z.
/// @param t        Per-sample readout time (s), length num_k.
/// @param num_k    Number of k-space samples.
/// @param num_i    Number of image pixels.
/// @param num_x    Number of R2* gradient samples in x.
/// @param num_y    Number of R2* gradient samples in y.
/// @param num_z    Number of R2* gradient samples in z.
template <typename T1>
void ftCpuWithGrads(T1 *kdata_r, T1 *kdata_i, const T1 *idata_r, const T1 *idata_i,
           const T1 *kx, const T1 *ky, const T1 *kz, const T1 *ix, const T1 *iy,
           const T1 *iz, const T1 *FM, const T1 *Gx, const T1 *Gy, const T1 *Gz,
           const T1 *t, const int num_k,
           const int num_i, const int num_x, const int num_y, const int num_z);

// Explicit Instantiations
extern template void ftCpuWithGrads<float>(float *, float *, const float *,
                                  const float *, const float *, const float *,
                                  const float *, const float *, const float *,
                                  const float *, const float *, const float *,
                                  const float *, const float *, const float *,
                                  const int, const int, const int, const int,
                                  const int);
extern template void ftCpuWithGrads<double>(double *, double *, const double *,
                                   const double *, const double *,
                                   const double *, const double *,
                                   const double *, const double *,
                                   const double *, const double *,
                                   const double *, 
                                   const double *, const double *,
                                   const double *, const int, const int,
                                   const int, const int, const int);
/*===========================================================================*/
/*                                                                           */
/*  Synopsis    [CPU kernel of the Inverse Fourier Transformation (IFT).] */
/*                                                                           */
/*  Description [] */
/*                                                                           */
/*===========================================================================*/
/// @brief Non-uniform adjoint DFT kernel with R2* field-map gradients (CPU / OpenACC).
///
/// Adjoint of ftCpuWithGrads; applies conjugate phase and R2* weighting.
///
/// @tparam T1      Floating-point precision type.
/// @param idata_r  Output image real part array, length num_i.
/// @param idata_i  Output image imaginary part array, length num_i.
/// @param kdata_r  Input k-space real part array, length num_k.
/// @param kdata_i  Input k-space imaginary part array, length num_k.
/// @param kx       k-space x-coordinates, length num_k.
/// @param ky       k-space y-coordinates, length num_k.
/// @param kz       k-space z-coordinates, length num_k.
/// @param ix       Image-space x-coordinates, length num_i.
/// @param iy       Image-space y-coordinates, length num_i.
/// @param iz       Image-space z-coordinates, length num_i.
/// @param FM       Off-resonance field map (rad/s), length num_i.
/// @param Gx       R2* gradient in x, length num_x.
/// @param Gy       R2* gradient in y, length num_y.
/// @param Gz       R2* gradient in z, length num_z.
/// @param t        Per-sample readout time (s), length num_k.
/// @param num_k    Number of k-space samples.
/// @param num_i    Number of image pixels.
/// @param num_x    Number of R2* gradient samples in x.
/// @param num_y    Number of R2* gradient samples in y.
/// @param num_z    Number of R2* gradient samples in z.
template <typename T1>
void iftCpuWithGrads(T1 *idata_r, T1 *idata_i, const T1 *kdata_r, const T1 *kdata_i,
            const T1 *kx, const T1 *ky, const T1 *kz, const T1 *ix,
            const T1 *iy, const T1 *iz, const T1 *FM, const T1 *Gx,
            const T1 *Gy, const T1 *Gz, const T1 *t,
            const int num_k, const int num_i, const int num_x, const int num_y,
            const int num_z);

// Explicit Instantiations
extern template void iftCpuWithGrads<float>(float *, float *, const float *,
                                   const float *, const float *, const float *,
                                   const float *, const float *, const float *,
                                   const float *, const float *, const float *,
                                   const float *, const float *, const float *,
                                   const int, const int, const int, const int,
                                   const int);
extern template void iftCpuWithGrads<double>(double *, double *, const double *,
                                    const double *, const double *,
                                    const double *, const double *,
                                    const double *, const double *,
                                    const double *, const double *,
                                    const double *,
                                    const double *, const double *,
                                    const double *, const int, const int,
                                    const int, const int, const int);
/*---------------------------------------------------------------------------*/
/*  Namespace declared - end                                                 */
/*---------------------------------------------------------------------------*/

//}
//}

#endif // FTWithGrads_CPU_H
