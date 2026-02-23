/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [gridding.h]

    Synopsis    [Implementation of the forward and adjoint non-uniform Fast
                                Fourier Transform (NUFFT) on CPU and GPU via
 OpenACC.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

/// @file gridding.h
/// @brief Non-uniform FFT gridding kernels (forward and adjoint) for CPU and GPU.

#ifndef PowerGrid_gridding_hpp
#define PowerGrid_gridding_hpp

#include <cstdlib>

#ifdef OPENACC_GPU // GPU Version
    #include "cufft.h"
    #include "FFT/fftGPU.h"
    #include "FFT/fftCPU.h"
    #include "griddingSupport.h"
    #include "Core/griddingTypes.h"
    #include "openacc.h"
    #define CFTHandle cufftHandle
#elif _OPENACC
    #include "FFT/fftCPU.h"
    #include "griddingSupport.h"
    #include "Core/griddingTypes.h"
    #define CFTHandle void
    #include "openacc.h"
#elif defined(METAL_COMPUTE) // Apple Metal path
    #include "Metal/MetalGridding.h"
    #include "Metal/MetalNufftPipeline.h"
    #include "FFT/fftAccelerate.h"
    #include "FFT/fftCPU.h"        // fallback for non-pow2 and double
    #include "griddingSupport.h"
    #include "Core/griddingTypes.h"
    #define CFTHandle MetalGriddingContext*
#else // CPU version
    #include "FFT/fftCPU.h"
    #include "griddingSupport.h"
    #include "Core/griddingTypes.h"
    #define CFTHandle void
#endif

using namespace arma;

/// @brief 2-D Kaiser-Bessel adjoint gridding (k-space -> oversampled grid).
///
/// @tparam T1       Floating-point precision type.
/// @param n         Number of k-space samples.
/// @param params    Gridding parameters (image/grid size, oversampling, kernel).
/// @param beta      Kaiser-Bessel kernel shape parameter beta.
/// @param sample    Array of ReconstructionSample structs (k-space value + coordinates).
/// @param LUT       Precomputed Kaiser-Bessel lookup table.
/// @param sizeLUT   Number of entries in the LUT.
/// @param gridData  Output oversampled grid (interleaved real/imag), length 2*gridX*gridY.
/// @returns         0 on success.
template <typename T1>
int gridding_adjoint_2D(unsigned int n, parameters<T1> params, T1 beta,
                        ReconstructionSample<T1> *__restrict sample,
                        const T1 *LUT, const uword sizeLUT,
                        T1 *__restrict gridData);

/// @brief 3-D Kaiser-Bessel adjoint gridding (k-space -> oversampled grid).
///
/// @tparam T1       Floating-point precision type.
/// @param n         Number of k-space samples.
/// @param params    Gridding parameters (image/grid size, oversampling, kernel).
/// @param beta      Kaiser-Bessel kernel shape parameter beta.
/// @param sample    Array of ReconstructionSample structs.
/// @param LUT       Precomputed Kaiser-Bessel lookup table.
/// @param sizeLUT   Number of entries in the LUT.
/// @param gridData  Output oversampled grid (interleaved real/imag), length 2*gridX*gridY*gridZ.
/// @returns         0 on success.
template <typename T1>
int gridding_adjoint_3D(unsigned int n, parameters<T1> params, T1 beta,
                        ReconstructionSample<T1> *__restrict sample,
                        const T1 *LUT, const uword sizeLUT,
                        T1 *gridData);

/// @brief 2-D Kaiser-Bessel forward gridding (oversampled grid -> k-space samples).
///
/// @tparam T1        Floating-point precision type.
/// @param n          Number of k-space samples.
/// @param params     Gridding parameters.
/// @param kx         k-space x-coordinates, length n.
/// @param ky         k-space y-coordinates, length n.
/// @param beta       Kaiser-Bessel kernel shape parameter beta.
/// @param pSamples   Output k-space samples array, length 2*n (interleaved).
/// @param LUT        Precomputed Kaiser-Bessel lookup table.
/// @param sizeLUT    Number of entries in the LUT.
/// @param pGridData  Input oversampled grid (interleaved real/imag).
/// @returns          0 on success.
template <typename T1>
int gridding_forward_2D(unsigned int n, parameters<T1> params, const T1 *kx,
                        const T1 *ky, T1 beta, T1 *__restrict pSamples,
                        const T1 *LUT, const uword sizeLUT,
                        T1 *__restrict pGridData);

/// @brief 3-D Kaiser-Bessel forward gridding (oversampled grid -> k-space samples).
///
/// @tparam T1        Floating-point precision type.
/// @param n          Number of k-space samples.
/// @param params     Gridding parameters.
/// @param kx         k-space x-coordinates, length n.
/// @param ky         k-space y-coordinates, length n.
/// @param kz         k-space z-coordinates, length n.
/// @param beta       Kaiser-Bessel kernel shape parameter beta.
/// @param pSamples   Output k-space samples array, length 2*n (interleaved).
/// @param LUT        Precomputed Kaiser-Bessel lookup table.
/// @param sizeLUT    Number of entries in the LUT.
/// @param pGridData  Input oversampled grid (interleaved real/imag).
/// @returns          0 on success.
template <typename T1>
int gridding_forward_3D(unsigned int n, parameters<T1> params, const T1 *kx,
                        const T1 *ky, const T1 *kz, T1 beta,
                        T1 *__restrict pSamples, const T1 *LUT,
                        const uword sizeLUT, T1 *__restrict pGridData);

/// @brief Full adjoint NUFFT pipeline: k-space -> image (grid, IFFT, crop, deapodize).
///
/// @tparam T1                Floating-point precision type.
/// @param numK_per_coil      Number of k-space samples per coil.
/// @param kx                 k-space x-coordinates.
/// @param ky                 k-space y-coordinates.
/// @param kz                 k-space z-coordinates.
/// @param dIn                Input k-space data (interleaved real/imag).
/// @param Nx                 Image size in x.
/// @param Ny                 Image size in y.
/// @param Nz                 Image size in z.
/// @param gridOS             Grid oversampling factor.
/// @param kernelWidth        Kaiser-Bessel kernel width.
/// @param beta               Kaiser-Bessel kernel shape parameter beta.
/// @param LUT                Precomputed Kaiser-Bessel lookup table.
/// @param sizeLUT            Number of LUT entries.
/// @param stream             OpenACC/CUDA stream (NULL for CPU).
/// @param plan               cuFFT plan pointer (NULL for CPU).
/// @param pGridData_crop_deAp  Scratch: deapodized cropped image.
/// @param pGridData_crop_d     Scratch: cropped grid data (GPU).
/// @param pGridData            Scratch: full oversampled grid.
/// @param pGridData_d          Scratch: GPU-side oversampled grid.
template <typename T1>
void computeFH_CPU_Grid(int numK_per_coil, const T1 *__restrict kx,
                        const T1 *__restrict ky, const T1 *__restrict kz,
                        const T1 *__restrict dIn,
                        int Nx, int Ny, int Nz, T1 gridOS,
                        const T1 kernelWidth, const T1 beta, const T1 *LUT,
                        const uword sizeLUT, void *stream, CFTHandle *plan,
                        T1 *pGridData_crop_deAp, T1 *pGridData_crop_d,
                        T1 *pGridData, T1 *pGridData_d);

/// @brief Full forward NUFFT pipeline: image -> k-space (deapodize, zero-pad, FFT, grid).
///
/// @tparam T1                Floating-point precision type.
/// @param numK_per_coil      Number of k-space samples per coil.
/// @param kx                 k-space x-coordinates.
/// @param ky                 k-space y-coordinates.
/// @param kz                 k-space z-coordinates.
/// @param dIn                Input image data (interleaved real/imag).
/// @param Nx                 Image size in x.
/// @param Ny                 Image size in y.
/// @param Nz                 Image size in z.
/// @param gridOS             Grid oversampling factor.
/// @param kernelWidth        Kaiser-Bessel kernel width.
/// @param beta               Kaiser-Bessel kernel shape parameter beta.
/// @param LUT                Precomputed Kaiser-Bessel lookup table.
/// @param sizeLUT            Number of LUT entries.
/// @param stream             OpenACC/CUDA stream (NULL for CPU).
/// @param plan               cuFFT plan pointer (NULL for CPU).
/// @param pGridData          Scratch: full oversampled grid (CPU).
/// @param pGridData_d        Scratch: GPU-side oversampled grid.
/// @param pGridData_os       Scratch: zero-padded oversampled grid (CPU).
/// @param pGridData_os_d     Scratch: GPU-side zero-padded grid.
/// @param pSamples           Output k-space samples array (interleaved real/imag).
template <typename T1>
void computeFd_CPU_Grid(int numK_per_coil, const T1 *__restrict kx,
                        const T1 *__restrict ky, const T1 *__restrict kz,
                        const T1 *__restrict dIn,
                        int Nx, int Ny, int Nz, T1 gridOS,
                        const T1 kernelWidth, const T1 beta, const T1 *LUT,
                        const uword sizeLUT, void *stream, CFTHandle *plan,
                        T1 *pGridData, T1 *pGridData_d, T1 *pGridData_os,
                        T1 *pGridData_os_d, T1 *pSamples);

// Explicit Instantiations
extern template int gridding_adjoint_2D<float>(unsigned int, parameters<float>,
                                               float,
                                               ReconstructionSample<float> *,
                                               const float *, const uword,
                                               float *);
extern template int gridding_adjoint_2D<double>(unsigned int,
                                                parameters<double>, double,
                                                ReconstructionSample<double> *,
                                                const double *, const uword,
                                                double *);
extern template int gridding_adjoint_3D<float>(unsigned int, parameters<float>,
                                               float,
                                               ReconstructionSample<float> *,
                                               const float *, const uword,
                                               float *);
extern template int gridding_adjoint_3D<double>(unsigned int,
                                                parameters<double>, double,
                                                ReconstructionSample<double> *,
                                                const double *, const uword,
                                                double *);
extern template int gridding_forward_2D<float>(unsigned int, parameters<float>,
                                               const float *, const float *,
                                               float beta, float *,
                                               const float *, const uword,
                                               float *);
extern template int
gridding_forward_2D<double>(unsigned int, parameters<double>, const double *,
                            const double *, double beta, double *,
                            const double *, const uword, double *);
extern template int gridding_forward_3D<float>(unsigned int, parameters<float>,
                                               const float *, const float *,
                                               const float *, float beta,
                                               float *, const float *,
                                               const uword, float *);
extern template int
gridding_forward_3D<double>(unsigned int, parameters<double>, const double *,
                            const double *, const double *, double beta,
                            double *, const double *, const uword,
                            double *);
extern template void
computeFH_CPU_Grid<float>(int, const float *, const float *, const float *,
                          const float *, int, int, int,
                          float gridOS, const float,
                          const float, const float *, const uword, void *,
                          CFTHandle *, float *, float *, float *, float *);
extern template void
computeFH_CPU_Grid<double>(int, const double *, const double *, const double *,
                           const double *, int, int, int,
                           double gridOS, const double,
                           const double, const double *, const uword, void *,
                           CFTHandle *, double *, double *, double *, double *);
extern template void 
computeFd_CPU_Grid<float>(int, const float *,
                                        const float *, const float *, const float *,
                                        int, int, int, float, const float, const float,
                                        const float *, const uword, void *, CFTHandle *,
                                        float *, float *, float *, float *, float *);
extern template void 
computeFd_CPU_Grid<double>(int, const double *,
                                        const double *, const double *, const double *,
                                        int, int, int, double, const double, const double,
                                        const double *, const uword, void *, CFTHandle *,
                                        double *, double *, double *, double *, double *);

#endif
