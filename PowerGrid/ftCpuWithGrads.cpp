/*
(C) Copyright 2010-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [ftCpuWithGrads.cpp]

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
    Date        [06/15/2018]

 *****************************************************************************/

/*---------------------------------------------------------------------------*/
/*  Included library headers                                                 */
/*---------------------------------------------------------------------------*/
#include "ftCpuWithGrads.h"

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
inline T1 sincPG(T1 x)

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
void ftCpuWithGrads(T1 *kdata_r, T1 *kdata_i, const T1 *idata_r, const T1 *idata_i,
           const T1 *kx, const T1 *ky, const T1 *kz, const T1 *ix, const T1 *iy,
           const T1 *iz, const T1 *FM, const T1 *Gx, const T1 *Gy, const T1 *Gz, const T1 *t, const int num_k,
           const int num_i, const int num_x, const int num_y, const int num_z) {

  T1 sumr = 0, sumi = 0, expr = 0,
      kxtpi = 0, bfunc = 0,
     kytpi = 0, kztpi = 0;

     float cosexpr = 0, sinexpr = 0;
  //T1   kzdeltaz = 0, kziztpi = 0, kx_N = 0, ky_N = 0, t_tpi = 0;
  int i = 0, j = 0;
  //tpi = 2 * MRI_PI;

// NON-conjugate transpose of G
#pragma acc kernels copyin(kx[0:num_k], ky[0:num_k], kz[0:num_k],ix[0:num_i], \
  iy[0:num_i], iz[0:num_i], FM[0:num_i], t[0:num_k], idata_r[0:num_i], idata_i[0:num_i]) \
  copyout(kdata_r[0:num_k], kdata_i[0:num_k], Gx[0:num_i], Gy[0:num_i], Gz[0:num_i])
  {

#pragma acc loop independent gang
    for (i = 0; i < num_k; i++) { // i is the time point in k-space
      sumr = 0.0;
      sumi = 0.0;

      kxtpi = kx[i] * 2 * MRI_PI;
      kytpi = ky[i] * 2 * MRI_PI;
      kztpi = kz[i] * 2 * MRI_PI;

      T1 myti = t[i];

#pragma acc loop vector(128)
      for (j = 0; j < num_i; j++) { // j is the pixel point in image-space
        expr = (kxtpi * ix[j] + kytpi * iy[j] + kztpi * iz[j] + (FM[j] * myti));
        
        sinexpr = sinf(expr);
        cosexpr = cosf(expr);
        //my_sincPGosf(expr, &sinexpr, &cosexpr);
        bfunc = sincPG(kx[i]/num_x + Gx[j] * myti) * sincPG(ky[i]/num_y + Gy[j] * myti) * sincPG(kz[i]/num_z + Gz[j] * myti);

        sumr += bfunc * ((cosexpr * idata_r[j]) + (sinexpr * idata_i[j]));
        sumi += bfunc * ((-sinexpr * idata_r[j]) + (cosexpr * idata_i[j]));
      }
      kdata_r[i] = sumr; // Real part
      kdata_i[i] = sumi; // Imaginary part
    }
  }
}

/*===========================================================================*/
/*                                                                           */
/*  Synopsis    [CPU kernel of the Inverse Fourier Transformation (IFT).]    */
/*                                                                           */
/*  Description []                                                           */
/*                                                                           */
/*===========================================================================*/
template <typename T1>
void iftCpuWithGrads(T1 *idata_r, T1 *idata_i, const T1 *kdata_r, const T1 *kdata_i,
            const T1 *kx, const T1 *ky, const T1 *kz, const T1 *ix,
            const T1 *iy, const T1 *iz, const T1 *FM, const T1 *Gx, const T1 *Gy, const T1 *Gz,const T1 *t,
            const int num_k, const int num_i, const int num_x, const int num_y, const int num_z) {

  T1 sumr = 0, sumi = 0, expr = 0, tpi = 0, 
     cosexpr = 0, sinexpr = 0, itraj_x_tpi = 0, itraj_y_tpi = 0,
     itraj_z_tpi = 0, bfunc = 0;
  //T1 kzdeltaz = 0, kziztpi = 0;
  int i = 0, j = 0;

  //--------------------------------------------------------------------
  //                         Initialization
  //--------------------------------------------------------------------
  tpi = MRI_PI * 2.0;

// kzdeltaz = kz[0] * MRI_DELTAZ;
// kziztpi = kz[0] * iz[0] * tpi;

//--------------------------------------------------------------------
//               Inverse Fourier Transform:     x=(G^H) * Gx
//--------------------------------------------------------------------
// the conjugate transpose of G
#if 0 // USE_OPENMP // FIXME: We can choose either this or the inner loop.
#pragma omp parallel for
#endif
#pragma acc kernels copyin(kx[0:num_k], ky[0:num_k], kz[0:num_k],ix[0:num_i], \
  iy[0:num_i], iz[0:num_i], FM[0:num_i], t[0:num_k], kdata_r[0:num_k], kdata_i[0:num_k]) \
  copyout( idata_r[0:num_i], idata_i[0:num_i], Gx[0:num_i], Gy[0:num_i], Gz[0:num_i])
  {
#pragma acc loop independent gang
    for (j = 0; j < num_i; j++) { // j is the pixel points in image-space
      sumr = 0.0;
      sumi = 0.0;

      itraj_x_tpi = ix[j] * tpi;
      itraj_y_tpi = iy[j] * tpi;
      itraj_z_tpi = iz[j] * tpi;
/*
#if USE_OPENMP
#pragma omp parallel for default(none) reduction(+ : sumr, sumi) private(      \
    expr, cosexpr, sinexpr) shared(j, kx, ky, kz, t, kzdeltaz, itraj_x_tpi,    \
                                   itraj_y_tpi, kziztpi, fm, kdata_r, kdata_i)
#endif
*/
      T1 myfmj = FM[j];


#pragma acc loop vector(128)
      for (i = 0; i < num_k; i++) { // i is the time points in k-space
        expr = (kx[i] * itraj_x_tpi + ky[i] * itraj_y_tpi +
                kz[i] * itraj_z_tpi + (myfmj * t[i]));

        // cosexpr = COS(expr); sinexpr = SIN(expr);

        sinexpr = sinf(expr);
        cosexpr = cosf(expr);
        bfunc = sincPG(kx[i]/num_x + Gx[j] * t[i]) * sincPG(ky[i]/num_y + Gy[j] * t[i]) * sincPG(kz[i]/num_z + Gz[j] * t[i]);
        sumr += bfunc * ((cosexpr * kdata_r[i]) - (sinexpr * kdata_i[i]));
        sumi += bfunc * ((sinexpr * kdata_r[i]) + (cosexpr * kdata_i[i]));
      }

      idata_r[j] = sumr; // Real part
      idata_i[j] = sumi; // Imaginary part
    }
  }
  // stopMriTimer(getMriTimer()->timer_iftCpu);
}

// Explicit Instantiations
template void ftCpuWithGrads<float>(float *, float *, const float *, const float *,
                           const float *, const float *, const float *,
                           const float *, const float *, const float *,
                           const float *, const float *, const float *,
                           const float *, const float *, const int, const int,
                           const int, const int, const int);
template void ftCpuWithGrads<double>(double *, double *, const double *, const double *,
                            const double *, const double *, const double *,
                            const double *, const double *, const double *,
                            const double *, const double *, const double *,
                            const double *, const double *, const int,
                            const int, const int, const int, const int);
template void iftCpuWithGrads<float>(float *, float *, const float *, const float *,
                            const float *, const float *, const float *,
                            const float *, const float *, const float *,
                            const float *, const float *, const float *,
                            const float *, const float *, const int, const int,
                            const int, const int, const int);
template void iftCpuWithGrads<double>(double *, double *, const double *, const double *,
                             const double *, const double *, const double *,
                             const double *, const double *, const double *,
                             const double *, const double *, const double *,
                             const double *, const double *, const int,
                             const int, const int, const int, const int);
//}
//}
