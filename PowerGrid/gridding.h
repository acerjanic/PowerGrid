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

#ifndef PowerGrid_gridding_hpp
#define PowerGrid_gridding_hpp

#include <cstdlib>

#ifdef OPENACC_GPU // GPU Version
    #include "cufft.h"
    #include "fftGPU.h"
    #include "fftCPU.h"
    #include "griddingSupport.h"
    #include "griddingTypes.h"
    #include "openacc.h"
    #define CFTHandle cufftHandle
#elif _OPENACC
    #include "fftCPU.h"
    #include "griddingSupport.h"
    #include "griddingTypes.h"
    #define CFTHandle void
    #include "openacc.h"
#else // CPU version
    #include "fftCPU.h"
    #include "griddingSupport.h"
    #include "griddingTypes.h"
    #define CFTHandle void
#endif

using namespace arma;

// 2D adjoint gridding on CPU
template <typename T1>
int gridding_adjoint_2D(unsigned int n, parameters<T1> params, T1 beta,
                        ReconstructionSample<T1> *__restrict sample,
                        const T1 *LUT, const uword sizeLUT,
                        T1 *__restrict gridData);
// 3D adjoint gridding on CPU
template <typename T1>
int gridding_adjoint_3D(unsigned int n, parameters<T1> params, T1 beta,
                        ReconstructionSample<T1> *__restrict sample,
                        const T1 *LUT, const uword sizeLUT,
                        T1 *gridData);

// 2D forward gridding on CPU
template <typename T1>
int gridding_forward_2D(unsigned int n, parameters<T1> params, const T1 *kx,
                        const T1 *ky, T1 beta, T1 *__restrict pSamples,
                        const T1 *LUT, const uword sizeLUT,
                        T1 *__restrict pGridData);

// 3D forward gridding on CPU
template <typename T1>
int gridding_forward_3D(unsigned int n, parameters<T1> params, const T1 *kx,
                        const T1 *ky, const T1 *kz, T1 beta,
                        T1 *__restrict pSamples, const T1 *LUT,
                        const uword sizeLUT, T1 *__restrict pGridData);

// Calculates the gridded adjoint transform
template <typename T1>
void computeFH_CPU_Grid(int numK_per_coil, const T1 *__restrict kx,
                        const T1 *__restrict ky, const T1 *__restrict kz,
                        const T1 *__restrict dIn,
                        int Nx, int Ny, int Nz, T1 gridOS,
                        const T1 kernelWidth, const T1 beta, const T1 *LUT,
                        const uword sizeLUT, void *stream, CFTHandle *plan,
                        T1 *pGridData_crop_deAp, T1 *pGridData_crop_d,
                        T1 *pGridData, T1 *pGridData_d);

// Calculates the gridded forward fourier transform
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
