/*
(C) Copyright 2010 The Board of Trustees of the University of Illinois.
All rights reserved.

Developed by:

                     IMPACT & MRFIL Research Groups
                University of Illinois, Urbana Champaign

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to
deal with the Software without restriction, including without limitation the
rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
sell copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimers.

Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimers in the documentation
and/or other materials provided with the distribution.

Neither the names of the IMPACT Research Group, MRFIL Research Group, the
University of Illinois, nor the names of its contributors may be used to
endorse or promote products derived from this Software without specific
prior written permission.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS WITH
THE SOFTWARE.
*/

/*****************************************************************************

    File Name   [ftCpu.h]

    Synopsis    [The CPU version of Fourier transform and inverse Fourier
        transform.]

    Description []

    Revision    [0.1; Initial build; Yue Zhuo, BIOE UIUC]
    Revision    [0.1.1; Add OpenMP, Code cleaning; Xiao-Long Wu, ECE UIUC]
    Revision    [1.0a; Further optimization, Code cleaning, Adding more comments;
                 Xiao-Long Wu, ECE UIUC, Jiading Gai, Beckman Institute]
    Date        [10/27/2010]

 *****************************************************************************/

#ifndef FT_CPU_H
#define FT_CPU_H

/*---------------------------------------------------------------------------*/
/*  Included library headers                                                 */
/*---------------------------------------------------------------------------*/
// Numeric constants according to the precision type.
#ifdef ENABLE_DOUBLE_PRECISION
#define MRI_PI               3.1415926535897932384626433832795029
#define MRI_NN               64
#define MRI_DELTAZ           0.003
#define MRI_ZERO             0.0
#define MRI_ONE              1.0
#define MRI_NEG_ONE         -1.0
#define MRI_POINT_FIVE       0.5
#define MRI_SMOOTH_FACTOR    0.0000001
#else
#define MRI_PI               3.1415926535897932384626433832795029f
#define MRI_NN               64
#define MRI_DELTAZ           0.003f
#define MRI_ZERO             0.0f
#define MRI_ONE              1.0f
#define MRI_NEG_ONE         -1.0f
#define MRI_POINT_FIVE       0.5f
#define MRI_SMOOTH_FACTOR    0.000001f
#endif
/*---------------------------------------------------------------------------*/
/*  Namespace declared - begin                                               */
/*---------------------------------------------------------------------------*/

//namespace uiuc_mri {

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

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <accelmath.h>

//#include <tools.h>
//#include <structures.h>

/*---------------------------------------------------------------------------*/
/*  Namespace declared - begin                                               */
/*---------------------------------------------------------------------------*/

//namespace uiuc_mri {

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

template <typename T1>
T1 sinc_cpu(T1 x)

{
    if (fabs(x) < 0.0001) {
        return 1.0;
    }
    
    return sin(MRI_PI * x) / (MRI_PI * x);
}

/*===========================================================================*/
/*                                                                           */
/*  Synopsis    [CPU kernel of the Fourier Transformation (FT).]             */
/*                                                                           */
/*  Description []                                                           */
/*                                                                           */
/*===========================================================================*/
template <typename T1>
void
ftCpu(T1 *kdata_r, T1 *kdata_i,
      const T1 *idata_r, const T1 *idata_i, const T1 *kx,
      const T1 *ky, const T1 *kz,
      const T1 *ix, const T1 *iy, const T1 *iz,
      const T1 *FM, const T1 *t,
      const int num_k, const int num_i
      )
{
    //startMriTimer(getMriTimer()->timer_ftCpu);
    
    T1 sumr = 0, sumi = 0, tpi = 0, kzdeltaz = 0, kziztpi = 0,
    expr = 0, cosexpr = 0, sinexpr = 0, t_tpi = 0,
    kx_N = 0, ky_N = 0, kxtpi = 0, kytpi = 0;
    int i = 0, j = 0;
    
    //--------------------------------------------------------------------
    //                         Initialization
    //--------------------------------------------------------------------
    tpi = 2 * MRI_PI;

    //kzdeltaz = kz[0] * MRI_DELTAZ;
    //kziztpi  = kz[0] * iz[0] * tpi;


    
    //--------------------------------------------------------------------
    //                    Fourier Transform:       Gx=G * x
    //--------------------------------------------------------------------
    // NON-conjugate transpose of G
#if 0 //USE_OPENMP // FIXME: We can choose either this or the inner loop.
#pragma omp parallel for
#endif
#pragma acc parallel loop
    for (i = 0; i < num_k; i++) { // i is the time point in k-space
        sumr = 0.0;
        sumi = 0.0;
        
        //t_tpi = t[i] / tpi;
        //kx_N = kx[i] / MRI_NN; // MRI_NN for normalized kx
        //ky_N = ky[i] / MRI_NN;
        kxtpi = kx[i] * tpi;
        kytpi = ky[i] * tpi;
        kztpi = kz[i]  * tpi;
        
#if USE_OPENMP
#pragma omp parallel for default(none) reduction(+:sumr, sumi) \
private(expr, cosexpr, sinexpr) \
shared(i, kx_N, ky_N, kzdeltaz, t_tpi, fm, kxtpi, kytpi, \
kziztpi, kx_i, ky_i, kz_i, t, idata_r, idata_i)
#endif
        T1 myti = t[i];
        for (j = 0; j < num_i; j++) { // j is the pixel point in image-space
            expr = (kxtpi * ix[j] + kytpi * iy[j] + kztpi * iz[j] +
                    (FM[j] * myti));

            cosexpr = cos(expr); sinexpr = sin(expr);

            //cosexpr = cosf(expr); sinexpr = sinf(expr);

            sumr += (cosexpr * idata_r[j]) + (sinexpr * idata_i[j]);
            sumi += (-sinexpr * idata_r[j]) + (cosexpr * idata_i[j]);
            
        }
        kdata_r[i] = sumr; // Real part
        kdata_i[i] = sumi; // Imaginary part
        //cout << "kdata[" << i << "] = " << sumr << " + 1j*" << sumi << endl ;
    }
    
    //stopMriTimer(getMriTimer()->timer_ftCpu);
}

/*===========================================================================*/
/*                                                                           */
/*  Synopsis    [CPU kernel of the Inverse Fourier Transformation (IFT).]    */
/*                                                                           */
/*  Description []                                                           */
/*                                                                           */
/*===========================================================================*/
template <typename T1>
void
iftCpu(T1 *idata_r, T1 *idata_i,
       const T1 *kdata_r, const T1 *kdata_i,
       const T1 *kx, const T1 *ky, const T1 *kz,
       const T1 *ix, const T1 *iy, const T1 *iz,
       const T1 *FM, const T1 *t,
       const int num_k, const int num_i
       )
{
    
    T1 sumr = 0, sumi = 0, tpi = 0, kzdeltaz = 0, kziztpi = 0,
    expr = 0, cosexpr = 0, sinexpr = 0,
    itraj_x_tpi = 0, itraj_y_tpi = 0;
    int i = 0, j = 0;
    
    //--------------------------------------------------------------------
    //                         Initialization
    //--------------------------------------------------------------------
    tpi = MRI_PI * 2.0;

    //kzdeltaz = kz[0] * MRI_DELTAZ;
    //kziztpi = kz[0] * iz[0] * tpi;

    
    //--------------------------------------------------------------------
    //               Inverse Fourier Transform:     x=(G^H) * Gx
    //--------------------------------------------------------------------
    // the conjugate transpose of G
#if 0 //USE_OPENMP // FIXME: We can choose either this or the inner loop.
#pragma omp parallel for
#endif
#pragma acc parallel loop
    for (j = 0; j < num_i; j++) { // j is the pixel points in image-space
        sumr = 0.0;
        sumi = 0.0;
        
        itraj_x_tpi = ix[j] * tpi;
        itraj_y_tpi = iy[j] * tpi;
        itraj_z_tpi = iz[j] * tpi;
        
        
#if USE_OPENMP
#pragma omp parallel for default(none) reduction(+:sumr, sumi) \
private(expr, cosexpr, sinexpr) \
shared(j, kx, ky, kz, t, kzdeltaz, itraj_x_tpi, itraj_y_tpi, kziztpi, \
fm, kdata_r, kdata_i)
#endif
	T1 myfmj = FM[j];
        for (i = 0; i < num_k; i++) { // i is the time points in k-space
            expr = (kx[i] * itraj_x_tpi +
                    ky[i] * itraj_y_tpi + kz[i] * itraj_z_tpi +
                    (myfmj * t[i]));

            cosexpr = cos(expr); sinexpr = sin(expr);

            //cosexpr = cosf(expr); sinexpr = sinf(expr);

            sumr += (cosexpr * kdata_r[i]) - (sinexpr * kdata_i[i]);
            sumi += (sinexpr * kdata_r[i]) + (cosexpr * kdata_i[i]);
            
        }
        
        idata_r[j] = sumr;  // Real part
        idata_i[j] = sumi;  // Imaginary part
    }
    
    //stopMriTimer(getMriTimer()->timer_iftCpu);
}

/*---------------------------------------------------------------------------*/
/*  Namespace declared - end                                                 */
/*---------------------------------------------------------------------------*/

//}
//}

#endif // FT_CPU_H

