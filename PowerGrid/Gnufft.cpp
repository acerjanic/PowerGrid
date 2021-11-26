/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [Gnufft.cpp]

    Synopsis    [Object that represents a non-uniform discrete Fourier
                    transform implemented via a a CPU and GPU accelerated
                    gridding implementation of the non-uniform Fast Fourier
                    Tranform (NUFFT).]

    Description [Forward transforms are denoted by G*data and adjoint transforms
                    are denoted by G/data. See documentation for more
                    information]

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

#include "Gnufft.h"
#include <chrono>
template <typename T1>
Gnufft<T1>::Gnufft(
    uword dataLength, T1 gridos, uword nx, uword ny, uword nz,
    const Col<T1> &k1, const Col<T1> &k2, const Col<T1> &k3, const Col<T1> &i1,
    const Col<T1> &i2,
    const Col<T1> &i3) // Change these arguments as you need to setup the object
{
  RANGE(__PRETTY_FUNCTION__)
  // cout << "Entering the class constructor for Ggrid" << endl;
  n1 = nx * ny * nz;
  n2 = dataLength;
  Nx = nx;
  Ny = ny;
  Nz = nz;
  ix = i1;
  iy = i2;
  iz = i3;
  //kx = k1;
  //ky = k2;
  //kz = k3;
  gridOS = gridos;

  imageNumElems = Nx * Ny;
  gridNumElems = gridOS * Nx * gridOS * Ny;

  if (Nz > 1) {
    imageNumElems = imageNumElems * Nz;
    gridNumElems = gridNumElems * gridOS * Nz;
  }

  kx = new T1[n2];
  ky = new T1[n2];
  kz = new T1[n2];

  for (int i = 0; i < n2; i++) {
    if (abs(k1(i)) > ((Nx / (T1)2.0) * 1.05) || abs(k2(i)) > ((Ny / (T1)2.0) * 1.05) ||
        abs(k3(i)) > ((Nz / (T1)2.0) * 1.05)) {
      
      cout << "Warning :k-space trajectory out of range [-N/2,N/2]:\n      "
            << "gridding requires that k-space should be contained within the "
            << "window -N/2 to N/2" << endl;
      cout << "kx = " << k1(i) << " ky = " << k2(i)<< " kz = " << k3(i)
           << " i = " << i << endl;
      }
      
      kx[i] = k1(i);
      ky[i] = k2(i);
      kz[i] = k3(i);
  }

  // Set Beta
  kernelWidth = 4.0;
  beta = MRI_PI *
         std::sqrt((gridOS - 0.5) * (gridOS - 0.5) *
                       (kernelWidth * kernelWidth * 4.0) / (gridOS * gridOS) -
                   0.8);
  #ifdef OPENACC_GPU
  stream = acc_get_cuda_stream(acc_async_sync);
  cufftCreate(&plan);
  if(Nz ==1) {
    if (cufftPlan2d(&plan, gridOS*Nx, gridOS*Ny, CUFFT_C2C) != CUFFT_SUCCESS) {
            cout <<  "CUFFT error: Plan creation failed" << endl;
    }
  } else {
    if (cufftPlan3d(&plan, gridOS*Nz, gridOS*Ny, gridOS*Nx, CUFFT_C2C) != CUFFT_SUCCESS) {
            cout << "CUFFT error: Plan creation failed" << endl;
          }
  }
  cufftSetStream(plan, (cudaStream_t)stream);
  #endif

  gridData = new complex<T1>[imageNumElems];
	gridData_d = new complex<T1>[imageNumElems];
  gridData_os = new complex<T1>[gridNumElems];
  gridData_os_d = new complex<T1>[gridNumElems];
  samples = new complex<T1>[n2];
  pGridData = reinterpret_cast<T1 *>(gridData);
  pGridData_d = reinterpret_cast<T1 *>(gridData_d);
  pGridData_os = reinterpret_cast<T1 *>(gridData_os);
  pGridData_os_d = reinterpret_cast<T1 *>(gridData_os_d);
  pSamples = reinterpret_cast<T1 *>(samples);
  // Deal with the LUT
  // Generating Look-Up Table
  // cout << "Calculating look up table" << endl;
  calculateLUT(beta, kernelWidth, LUT, sizeLUT);

#pragma acc enter data copyin(LUT[0 : sizeLUT], kx[0:n2], ky[0:n2], \
  kz[0:n2]) create(pGridData[0:2*imageNumElems], pGridData_d[0:2*imageNumElems], \
  pGridData_os[0:2*gridNumElems], pGridData_os_d[0:2*gridNumElems], pSamples[0:2*n2])
}

// Class destructor to free LUT
template <typename T1> Gnufft<T1>::~Gnufft() {
  RANGE(__PRETTY_FUNCTION__)
  #ifdef OPENACC_GPU
    cufftDestroy(plan);
  #endif
  #pragma acc exit data delete(pGridData[0:2*imageNumElems], \
   pGridData_d[0:2*imageNumElems], pGridData_os[0:2*gridNumElems], \
   pGridData_os_d[0:2*gridNumElems], kx[0:n2], ky[0:n2], kz[0:n2], \
   pSamples[0:2*n2])

  delete[] samples;
  delete[] gridData;
  delete[] gridData_d;
  delete[] gridData_os;
  delete[] gridData_os_d;
  delete[] kx;
  delete[] ky;
  delete[] kz;

  if (LUT) {
#pragma acc exit data delete (LUT)
    free(LUT);
  }
}

// Overloaded methods for forward and adjoint transform
// Forward transform operation using gridding
template <typename T1>
inline Col<complex<T1>> Gnufft<T1>::
operator*(const Col<complex<T1>> &d) const // Don't change these arguments
{
  RANGE(__PRETTY_FUNCTION__)

  const T1 *dataPtr = reinterpret_cast<const T1 *>(d.memptr());

  // T2 gridOS = 2.0;
  // cout << "About to call the forward gridding routine." << endl;
  #ifdef OPENACC_GPU
    cufftHandle *nPlan = const_cast<cufftHandle *>(&plan);
  #else
    void* nPlan = NULL;
  #endif
  computeFd_CPU_Grid(n2, kx, ky, kz, dataPtr,
                         Nx, Ny, Nz, gridOS,
                         kernelWidth, beta, LUT, sizeLUT,
                         stream, nPlan, pGridData, pGridData_d, pGridData_os,
                         pGridData_os_d, pSamples);

  Col<CxT1> temp(reinterpret_cast<CxT1 *>(pSamples), n2, false, true);
  return temp; // Return a vector of type T1
}

// Adjoint transform operation
template <typename T1>
inline Col<complex<T1>> Gnufft<T1>::operator/(const Col<complex<T1>> &d) const {
  // uword dataLength = n2;
  // Let's trim the operations to avoid data overhead and transfers
  // Basically if we know that the data points are zero, they have no impact
  // on the transform
  RANGE(__PRETTY_FUNCTION__)
  uword dataLength = this->n2;

  //Col<T1> realData = real(d).eval();
  //Col<T1> imagData = imag(d).eval();

  const T1 *dataPtr = reinterpret_cast<const T1 *>(d.memptr());

  // Process data here, like calling a brute force transform, dft...
  // I assume you create the pointers to the arrays where the transformed data
  // will be stored
  // realXformedDataPtr and imagXformedDataPtr and they are of type float*

  // T2 gridOS = 2.0;
  #ifdef OPENACC_GPU
  cufftHandle *nPlan = const_cast<cufftHandle *>(&plan);
  #else
  void *nPlan = NULL;
  #endif
  computeFH_CPU_Grid(dataLength, kx, ky, kz,
                         dataPtr, Nx, Ny, Nz, gridOS,
                         kernelWidth,
                         beta, LUT, sizeLUT, stream, nPlan, pGridData,
                         pGridData_d, pGridData_os, pGridData_os_d);
  
  Col<CxT1> temp(reinterpret_cast<CxT1 *>(pGridData), n1, false, true);
  return temp; // Return a vector of type T1
}

template <typename T1>
inline Col<complex<T1>> Gnufft<T1>::forwardSpatialInterp(const Col<complex<T1>> &d) const {
  uword dataLength = this->n2;

  const T1 *dataPtr = reinterpret_cast<const T1 *>(d.memptr());

  parameters<T1> params;
  params.sync = 0;
  params.binsize = 128;
  params.useLUT = 1;
  params.kernelWidth = kernelWidth;
  params.gridOS = gridOS;
  params.imageSize[0] = Nx; // gridSize is gridOS times larger than imageSize.
  params.imageSize[1] = Ny;
  params.imageSize[2] = Nz;
  params.gridSize[0] = std::ceil(gridOS * (T1)Nx);
  params.gridSize[1] = std::ceil(gridOS * (T1)Ny);
  if (params.gridSize[0] % 2) // 3D case, gridOS is adjusted on the z dimension:
    params.gridSize[0] += 1; // That why we need to make sure here that the xy
  if (params.gridSize[1] % 2) // dimensions have even sizes.
    params.gridSize[1] += 1;
  params.gridSize[2] = (Nz == 1) ? Nz : (std::ceil(gridOS * (T1)Nz)); // 2D or 3D
  params.numSamples = dataLength;

  unsigned int n = params.numSamples;
  
  memcpy(pGridData_os, dataPtr, sizeof(T1) * 2 * gridNumElems);
  #pragma acc update device(pGridData_os[0:2 * gridNumElems])

  //#pragma acc parallel loop present(pSamples [0:2 * n])
    for (int ii = 0; ii < 2 * n; ii++) {
        pSamples[ii] = (T1)0.0;
    }
    #pragma acc update device(pSamples [0:2 * n])

    if (Nz == 1) {
      gridding_forward_2D(n, params, kx,
                         ky, beta, pSamples, LUT,
                          sizeLUT, pGridData_os);
    } else {
      gridding_forward_3D(n, params, kx,
                      ky, kz, beta, pSamples, LUT,
                      sizeLUT, pGridData_os);
    }

  #pragma acc update host(pSamples [0:2 * n])

  Col<CxT1> temp(reinterpret_cast<CxT1 *>(pSamples), n, false, true);
  return temp; // Return a vector of type T1

}

template <typename T1>
inline Col<complex<T1>> Gnufft<T1>::adjointSpatialInterp(const Col<complex<T1>> &d) const {

  uword dataLength = this->n2;

  const T1 *dataPtr = reinterpret_cast<const T1 *>(d.memptr());

  parameters<T1> params;
  params.sync = 0;
  params.binsize = 128;
  params.useLUT = 1;
  params.kernelWidth = kernelWidth;
  params.gridOS = gridOS;
  params.imageSize[0] = Nx; // gridSize is gridOS times larger than imageSize.
  params.imageSize[1] = Ny;
  params.imageSize[2] = Nz;
  params.gridSize[0] = std::ceil(gridOS * (T1)Nx);
  params.gridSize[1] = std::ceil(gridOS * (T1)Ny);
  if (params.gridSize[0] % 2) // 3D case, gridOS is adjusted on the z dimension:
    params.gridSize[0] += 1; // That why we need to make sure here that the xy
  if (params.gridSize[1] % 2) // dimensions have even sizes.
    params.gridSize[1] += 1;
  params.gridSize[2] = (Nz == 1) ? Nz : (std::ceil(gridOS * (T1)Nz)); // 2D or 3D
  params.numSamples = dataLength;

  unsigned int n = params.numSamples;

  ReconstructionSample<T1>* samples; // Input Data
  // allocate samples
  samples = (ReconstructionSample<T1>*)malloc(
  params.numSamples * sizeof(ReconstructionSample<T1>));

  if (samples == NULL) {
    printf("ERROR: Unable to allocate memory for input data\n");
    exit(1);
  }

  //
  for (int i = 0; i < params.numSamples; i++) {

    samples[i].kX = kx[i];
    samples[i].kY = ky[i];
    samples[i].kZ = kz[i];

    samples[i].real = dataPtr[2 * i];
    samples[i].imag = dataPtr[2 * i + 1];

    //samples[i].sdc = (T1)1.0;
    // samples[i].t = t[i];
  }

  #pragma acc enter data copyin(samples [0:n])

  //#pragma acc parallel loop 
    for (int i = 0; i < gridNumElems; i++) {
        pGridData_os[2 * i] = (T1)0.0;
        pGridData_os[2 * i + 1] = (T1)0.0;
    }
  
  #pragma acc update device(pGridData_os[0:2*gridNumElems])

  if (Nz == 1) {
    gridding_adjoint_2D(n, params, beta, samples,
                        LUT, sizeLUT, pGridData_os);
  } else {
    gridding_adjoint_3D(n, params, beta, samples,
                        LUT, sizeLUT, pGridData_os);
  }

  #pragma acc update host(pGridData_os[0:2*gridNumElems])
  Col<CxT1> temp(reinterpret_cast<CxT1 *>(pGridData_os), gridNumElems, false, true);
  return temp; // Return a vector of type T1
}

// 2D adjoint gridding on CPU
template <typename T1>
int Gnufft<T1>::gridding_adjoint_2D(unsigned int n, parameters<T1> params, T1 beta,
    ReconstructionSample<T1>* __restrict sample,
    const T1* LUT, const uword sizeLUT,
    T1* __restrict pGData) const
{
    RANGE(__PRETTY_FUNCTION__)
    unsigned int NxL, NxH;
    unsigned int NyL, NyH;

    unsigned int nx;
    unsigned int ny;

    int idx;

    T1 w;

    T1 shiftedKx, shiftedKy;
    T1 distX, kbX, distY, kbY;


    T1 kernelWidth = params.kernelWidth;
    T1 gridOS = params.gridOS;

    unsigned int Nx = params.imageSize[0];
    unsigned int Ny = params.imageSize[1];
    int gridNumElems = params.gridSize[0] * params.gridSize[1];

    T1 kernelScale = ((T1)4.0 / (kernelWidth * kernelWidth) * (T1)sizeLUT);

#pragma acc parallel loop independent gang vector present(LUT [0:sizeLUT], \
    pGData [0:gridNumElems * 2], sample [0:n]) copyin(params, params.gridSize[0:3])
    for (int i = 0; i < n; i++) {
        ReconstructionSample<T1> pt = sample[i];

        shiftedKx = (gridOS) * (pt.kX + ((T1)Nx) / (T1)2.0);
        shiftedKy = (gridOS) * (pt.kY + ((T1)Ny) / (T1)2.0);

        NxL = (int)(std::max((T1)0.0, std::ceil(shiftedKx - kernelWidth * (gridOS) / (T1)2.0)));
        NxH = (int)(std::min((gridOS * (T1)Nx - (T1)1.0),
            std::floor(shiftedKx + kernelWidth * (gridOS) / (T1)2.0)));

        NyL = (int)(std::max((T1)0.0, std::ceil(shiftedKy - kernelWidth * (gridOS) / (T1)2.0)));
        NyH = (int)(std::min((gridOS * (T1)Ny - (T1)1.0),
            std::floor(shiftedKy + kernelWidth * (gridOS) / (T1)2.0)));

#pragma acc loop seq
        for (nx = NxL; nx <= NxH; ++nx) {
            int k0;
            distX = (shiftedKx - ((T1)nx)) / (gridOS);

            k0 = (int)(distX * distX * kernelScale);
            if (k0 >= sizeLUT)
                kbX = (T1)0.0;
            else
                kbX = LUT[k0];

#pragma acc loop seq
            for (ny = NyL; ny <= NyH; ++ny) {
                distY = (shiftedKy - ((T1)ny)) / (gridOS);

                k0 = (int)(distY * distY * kernelScale);

                if (k0 >= sizeLUT)
                    kbY = (T1)0.0;
                else
                    kbY = LUT[k0];

                w = kbX * kbY;

                /* grid data */
                idx = 2 * (ny + (nx)*params.gridSize[1]);

#pragma acc atomic update
                pGData[idx] += w * pt.real;

#pragma acc atomic update
                pGData[idx + 1] += w * pt.imag;
            }
        }
    }

    // PowerGrid uses: y->x->z because we are column major same as MATLAB...

    return 1;
}

// 3D adjoint gridding on CPU
template <typename T1>
int Gnufft<T1>::gridding_adjoint_3D(unsigned int n, parameters<T1> params, T1 beta,
    ReconstructionSample<T1>* __restrict sample,
    const T1* LUT, const uword sizeLUT,
    T1* pGData) const
{
    RANGE(__PRETTY_FUNCTION__)
    int NxL, NxH;
    int NyL, NyH;
    int NzL, NzH;

    int nx;
    int ny;
    int nz;

    int idx;

    T1 w;

    T1 shiftedKx, shiftedKy, shiftedKz;
    T1 distX, kbX, distY, kbY, distZ, kbZ;

    T1 kernelWidth = params.kernelWidth;
    T1 gridOS = params.gridOS;

    int Nx = params.imageSize[0];
    int Ny = params.imageSize[1];
    int Nz = params.imageSize[2];
    int gridNumElems = params.gridSize[0] * params.gridSize[1] * params.gridSize[2];


#pragma acc parallel loop gang vector pcopy(LUT [0:sizeLUT], pGData [0:gridNumElems * 2]) \
    pcopyin(params, params.gridSize [0:3], sample [0:n])
    for (int i = 0; i < n; i++) {
        ReconstructionSample<T1> pt = sample[i];

        shiftedKx = (gridOS) * (pt.kX + ((T1)Nx) / (T1)2.0);
        shiftedKy = (gridOS) * (pt.kY + ((T1)Ny) / (T1)2.0);
        shiftedKz = (gridOS) * (pt.kZ + ((T1)Nz) / (T1)2.0);

        NxL = (int)(std::max((T1)0.0, std::ceil(shiftedKx - kernelWidth * (gridOS) / (T1)2.0)));
        NxH = (int)(std::min((gridOS * (T1)Nx - (T1)1.0),
            std::floor(shiftedKx + kernelWidth * (gridOS) / (T1)2.0)));

        NyL = (int)(std::max((T1)0.0, std::ceil(shiftedKy - kernelWidth * (gridOS) / (T1)2.0)));
        NyH = (int)(std::min((gridOS * (T1)Ny - (T1)1.0),
            std::floor(shiftedKy + kernelWidth * (gridOS) / (T1)2.0)));

        NzL = (int)(std::max((T1)0.0, std::ceil(shiftedKz - kernelWidth * (gridOS) / (T1)2.0)));
        NzH = (int)(std::min((gridOS * (T1)Nz - (T1)1.0),
            std::floor(shiftedKz + kernelWidth * (gridOS) / (T1)2.0)));
#pragma acc loop seq
        for (nz = NzL; nz <= NzH; ++nz) {
            int k0;
            distZ = std::abs(shiftedKz - ((T1)nz)) / (gridOS);

            k0 = (int)((distZ * distZ * (T1)4.0 / (kernelWidth * kernelWidth)) * (T1)sizeLUT);
            if (k0 >= sizeLUT)
                kbZ = (T1)0.0;
            else
                kbZ = LUT[k0];

#pragma acc loop seq
            for (nx = NxL; nx <= NxH; ++nx) {
                distX = std::abs(shiftedKx - ((T1)nx)) / (gridOS);

                k0 = (int)((distX * distX * (T1)4.0 / (kernelWidth * kernelWidth)) * (T1)sizeLUT);
                if (k0 >= sizeLUT)
                    kbX = (T1)0.0;
                else
                    kbX = LUT[k0];

#pragma acc loop seq
                for (ny = NyL; ny <= NyH; ++ny) {
                    distY = std::abs(shiftedKy - ((T1)ny)) / (gridOS);

                    k0 = (int)((distY * distY * (T1)4.0 / (kernelWidth * kernelWidth)) * (T1)sizeLUT);
                    if (k0 >= sizeLUT)
                        kbY = (T1)0.0;
                    else
                        kbY = LUT[k0];

                    w = kbX * kbY * kbZ;

                    /* grid data */
                    idx = ny + (nx)*params.gridSize[1] + (nz)*params.gridSize[0] * params.gridSize[1];
#pragma acc atomic update
                    pGData[2 * idx] += w * pt.real;

#pragma acc atomic update
                    pGData[2 * idx + 1] += w * pt.imag;
                }
            }
        }
    }

    return 1;
}

// 2D forward gridding on CPU
template <typename T1>
int Gnufft<T1>::gridding_forward_2D(unsigned int n, parameters<T1> params, const T1* kx,
    const T1* ky, T1 beta, T1* __restrict pSamples,
    const T1* LUT, const uword sizeLUT,
    T1* __restrict pGridData) const
{
    RANGE(__PRETTY_FUNCTION__)
    int NxL, NxH;
    int NyL, NyH;

    int nx;
    int ny;

    int idx;
    T1 w;
    T1 shiftedKx, shiftedKy;
    T1 distX, kbX, distY, kbY;

    T1 kernelWidth = params.kernelWidth;
    T1 gridOS = params.gridOS;

    int Nx = params.imageSize[0];
    int Ny = params.imageSize[1];
    int gridNumElems = params.gridSize[0] * params.gridSize[1];
    
    T1 kernelScale = ((T1)4.0 / (kernelWidth * kernelWidth) * (T1)sizeLUT);

#pragma acc parallel loop gang vector present(kx [0:n], ky [0:n], pSamples [0:n * 2], \
    LUT [0:sizeLUT], pGridData [0:gridNumElems * 2]) copyin(params, params.gridSize[0:3])
    for (int i = 0; i < n; i++) {

        shiftedKx = (gridOS) * (kx[i] + ((T1)Nx) / (T1)2.0);
        shiftedKy = (gridOS) * (ky[i] + ((T1)Ny) / (T1)2.0);

        NxL = (int)(std::max((T1)0.0, std::ceil(shiftedKx - kernelWidth * (gridOS) / (T1)2.0)));
        NxH = (int)(std::min((gridOS * (T1)Nx - (T1)1.0),
            std::floor(shiftedKx + kernelWidth * (gridOS) / (T1)2.0)));

        NyL = (int)(std::max((T1)0.0, std::ceil(shiftedKy - kernelWidth * (gridOS) / (T1)2.0)));
        NyH = (int)(std::min((gridOS * (T1)Ny - (T1)1.0),
            std::floor(shiftedKy + kernelWidth * (gridOS) / (T1)2.0)));

#pragma acc loop seq
        for (nx = NxL; nx <= NxH; ++nx) {
            int k0;
            distX = (shiftedKx - ((T1)nx)) / (gridOS);

            k0 = (int)(distX * distX * kernelScale);
            if (k0 >= sizeLUT)
                kbX = (T1)0.0;
            else
                kbX = LUT[k0];

#pragma acc loop seq
            for (ny = NyL; ny <= NyH; ++ny) {
                distY = (shiftedKy - ((T1)ny)) / (gridOS);

                k0 = (int)(distY * distY * kernelScale);
                if (k0 >= sizeLUT)
                    kbY = (T1)0.0;
                else
                    kbY = LUT[k0];

                w = kbX * kbY;

                /* grid data */
                idx = ny + (nx)*params.gridSize[1];

                //#pragma acc atomic update
                pSamples[2 * i] += w * pGridData[2 * idx];

                //#pragma acc atomic update
                pSamples[2 * i + 1] += w * pGridData[2 * idx + 1];
            }
        }
    }

    // re-arrange dimensions and output
    // Nady uses: x->y->z
    // IMPATIENT uses: z->x->y
    // PowerGrid uses: x->y->z because we are column major same as Nady...
    // Nope! XX So we need to convert from (x->y->z)-order to (z->x->y)-order
    // int gridNumElems = params.gridSize[0] * params.gridSize[1];

    return 1;
}

// 3D forward gridding on CPU
template <typename T1>
int Gnufft<T1>::gridding_forward_3D(unsigned int n, parameters<T1> params, const T1* kx,
    const T1* ky, const T1* kz, T1 beta,
    T1* __restrict pSamples, const T1* LUT,
    const uword sizeLUT, T1* __restrict pGridData) const
{
    RANGE(__PRETTY_FUNCTION__)
    int NxL, NxH;
    int NyL, NyH;
    int NzL, NzH;

    int nx;
    int ny;
    int nz;

    int idx;
    //T1 *pSamples;
    //T1 *pGridData;
    T1 w;

    T1 shiftedKx, shiftedKy, shiftedKz;
    T1 distX, kbX, distY, kbY, distZ, kbZ;

    T1 kernelWidth = params.kernelWidth;
    // T1 beta = 18.5547;
    T1 gridOS = params.gridOS;

    int Nx = params.imageSize[0];
    int Ny = params.imageSize[1];
    int Nz = params.imageSize[2];
    //int imageNumElems =
    //    params.imageSize[0] * params.imageSize[1] * params.imageSize[2];
    int gridNumElems = params.gridSize[0] * params.gridSize[1] * params.gridSize[2];
    //pGridData = reinterpret_cast<T1 *>(gridData);
    // Jiading GAI
    // float t0 = t[0];

#pragma acc parallel loop gang vector present(LUT [0:sizeLUT],    \
    pGridData [0:gridNumElems * 2], kx [0:n], ky [0:n], kz [0:n], \
    pSamples [0:n * 2]) copyin(params, params.gridSize[0:3])
    for (int i = 0; i < n; i++) {
        // complex<T1> pt = sample[i];

        // Jiading GAI
        // float atm = hanning_d(t[i], tau, l, t0);//a_l(t_m)

        shiftedKx = (gridOS) * (kx[i] + ((T1)Nx) / (T1)2.0);
        shiftedKy = (gridOS) * (ky[i] + ((T1)Ny) / (T1)2.0);
        shiftedKz = (gridOS) * (kz[i] + ((T1)Nz) / (T1)2.0);

        NxL = (int)(std::max((T1)0.0, std::ceil(shiftedKx - kernelWidth * (gridOS) / (T1)2.0)));
        NxH = (int)(std::min((gridOS * (T1)Nx - (T1)1.0),
            std::floor(shiftedKx + kernelWidth * ((T1)gridOS) / (T1)2.0)));

        NyL = (int)(std::max((T1)0.0, std::ceil(shiftedKy - kernelWidth * (gridOS) / (T1)2.0)));
        NyH = (int)(std::min((gridOS * (T1)Ny - (T1)1.0),
            std::floor(shiftedKy + kernelWidth * ((T1)gridOS) / (T1)2.0)));

        NzL = (int)(std::max((T1)0.0, std::ceil(shiftedKz - kernelWidth * (gridOS) / 2.0f)));
        NzH = (int)(std::min((gridOS * (T1)Nz - (T1)1.0),
            std::floor(shiftedKz + kernelWidth * ((T1)gridOS) / (T1)2.0)));

#pragma acc loop seq
        for (nz = NzL; nz <= NzH; ++nz) {
            int k0;
            distZ = std::abs(shiftedKz - ((T1)nz)) / (gridOS);

            //if (params.useLUT) {
            // kbZ = kernel_value_LUT(distZ, LUT, sizeLUT, kernelWidth);
            k0 = (int)((distZ * distZ * (T1)4.0 / (kernelWidth * kernelWidth)) * (T1)sizeLUT);
            if (k0 >= sizeLUT)
                kbZ = (T1)0.0;
            else
                kbZ = LUT[k0];

#pragma acc loop seq
            for (nx = NxL; nx <= NxH; ++nx) {
                distX = std::abs(shiftedKx - ((T1)nx)) / (gridOS);

                k0 = (int)((distX * distX * (T1)4.0 / (kernelWidth * kernelWidth)) * (T1)sizeLUT);
                if (k0 >= sizeLUT)
                    kbX = (T1)0.0;
                else
                    kbX = LUT[k0];

#pragma acc loop seq
                for (ny = NyL; ny <= NyH; ++ny) {
                    distY = std::abs(shiftedKy - ((T1)ny)) / (gridOS);

                    k0 = (int)((distY * distY * (T1)4.0 / (kernelWidth * kernelWidth)) * (T1)sizeLUT);
                    if (k0 >= sizeLUT)
                        kbY = (T1)0.0;
                    else
                        kbY = LUT[k0];

                    w = kbX * kbY * kbZ;

                    /* grid data */
                    idx = ny + (nx)*params.gridSize[1] + (nz)*params.gridSize[0] * params.gridSize[1];

//#pragma acc atomic update
                    pSamples[2 * i] += w * pGridData[2 * idx];

//#pragma acc atomic update
                    pSamples[2 * i + 1] += w * pGridData[2 * idx + 1];
                }
            }
        }
    }

    return 1;
}

// Calculates the gridded adjoint transform
template <typename T1>
void Gnufft<T1>::computeFH_CPU_Grid(int numK_per_coil, const T1* __restrict kx,
    const T1* __restrict ky, const T1* __restrict kz,
    const T1* __restrict dIn,
    int Nx, int Ny, int Nz, T1 gridOS,
    const T1 kernelWidth, const T1 beta, const T1* LUT,
    const uword sizeLUT, void* stream, CFTHandle* plan,
    T1* pGridData_crop_deAp, T1* pGridData_crop_d,
    T1* pGridData, T1* pGridData_d) const
{

    RANGE(__PRETTY_FUNCTION__)
    parameters<T1> params;
    params.sync = 0;
    params.binsize = 128;

    params.useLUT = 1;
    params.kernelWidth = kernelWidth;
    params.gridOS = gridOS;
    params.imageSize[0] = Nx; // gridSize is gridOS times larger than imageSize.
    params.imageSize[1] = Ny;
    params.imageSize[2] = Nz;
    params.gridSize[0] = std::ceil(gridOS * (T1)Nx);
    params.gridSize[1] = std::ceil(gridOS * (T1)Ny);
    if (params.gridSize[0] % 2) // 3D case, gridOS is adjusted on the z dimension:
        params.gridSize[0] += 1; // That why we need to make sure here that the xy
    if (params.gridSize[1] % 2) // dimensions have even sizes.
        params.gridSize[1] += 1;
    params.gridSize[2] = (Nz == 1) ? Nz : (std::ceil(gridOS * (T1)Nz)); // 2D or 3D
    params.numSamples = numK_per_coil;

    ReconstructionSample<T1>* samples; // Input Data
    // allocate samples
    samples = (ReconstructionSample<T1>*)malloc(
        params.numSamples * sizeof(ReconstructionSample<T1>));

    if (samples == NULL) {
        printf("ERROR: Unable to allocate memory for input data\n");
        exit(1);
    }
    unsigned int n = params.numSamples;
    //
    for (int i = 0; i < params.numSamples; i++) {

        samples[i].kX = kx[i];
        samples[i].kY = ky[i];
        samples[i].kZ = kz[i];

        samples[i].real = dIn[2 * i];
        samples[i].imag = dIn[2 * i + 1];

        //samples[i].sdc = (T1)1.0;
        // samples[i].t = t[i];
    }
    // grid_size in xy-axis has to be divisible-by-two:
    //       (required by the cropImageRegion)
    // grid_size in z-axis has to be divisible-by-four:
    //       (required by the function gridding_GPU_3D(.))
    if (1 == Nz) {
        // round grid size (xy-axis) to the next divisible-by-two.
        gridOS = (T1)2.0 * std::ceil((gridOS * (T1)Nx) / (T1)2.0) / (T1)Nx;
    } else {
        // round grid size (z-axis) to the next divisible-by-four.
        gridOS = (T1)4.0 * std::ceil((gridOS * (T1)Nz) / (T1)4.0) / (T1)Nz;
    }

    int gridNumElems = params.gridSize[0] * params.gridSize[1] * params.gridSize[2];

    int imageNumElems = params.imageSize[0] * params.imageSize[1] * params.imageSize[2];

    // Have to set 'gridData' to zero.
    // Because they will be involved in accumulative operations
    // inside gridding functions.
#pragma acc enter data copyin(samples [0:n])

#pragma acc parallel loop
    for (int i = 0; i < gridNumElems; i++) {
        pGridData[2 * i] = (T1)0.0;
        pGridData[2 * i + 1] = (T1)0.0;
    }
    // Gridding with CPU - adjoint
    if (Nz == 1) {
        gridding_adjoint_2D(n, params, beta, samples, LUT, sizeLUT, pGridData);
    } else {
        gridding_adjoint_3D(n, params, beta, samples, LUT, sizeLUT, pGridData);
    }

    if (Nz == 1) {
        ifftshift2<T1>(pGridData_d, pGridData, params.gridSize[0],
            params.gridSize[1]);
    } else {
        ifftshift3<T1>(pGridData_d, pGridData, params.gridSize[0],
            params.gridSize[1], params.gridSize[2]);
    }

//#ifdef _OPENACC // We're on GPU
#ifdef OPENACC_GPU // We're on GPU
    // Inside this region the device data pointer will be used for cuFFT

#pragma acc host_data use_device(pGridData_d)
    {
        acc_set_cuda_stream(1, stream);

        if (Nz == 1) {
            ifft2dGPU(pGridData_d, params.gridSize[0], params.gridSize[1], stream, plan);
        } else {
            ifft3dGPU(pGridData_d, params.gridSize[0], params.gridSize[1],
                params.gridSize[2], stream, plan);
        }
    }

#else // We're on CPU so we'll use FFTW

    // Launch FFT on the CPU
    if (Nz == 1) {
        ifft2dCPU(pGridData_d, params.gridSize[0], params.gridSize[1]);
    } else {
        ifft3dCPU(pGridData_d, params.gridSize[0], params.gridSize[1],
            params.gridSize[2]);
    }

#endif
    /* 
    //#pragma acc update device(pGridData_d[0:2*gridNumElems])
    if (Nz == 1) {
        normalize_fft2d<T1>(pGridData, pGridData_d, params.gridSize[0],
            params.gridSize[1]);
    } else {
        normalize_fft3d<T1>(pGridData, pGridData_d, params.gridSize[0],
            params.gridSize[1], params.gridSize[2]);
    }
    */
    if (Nz == 1) {
        fftshift2<T1>(pGridData, pGridData_d, params.gridSize[0],
            params.gridSize[1]);
    } else {
        fftshift3<T1>(pGridData, pGridData_d, params.gridSize[0],
            params.gridSize[1], params.gridSize[2]);
    }

    if (Nz == 1) {
        crop_center_region2d<T1>(pGridData_crop_d, pGridData, params.imageSize[0],
            params.imageSize[1], params.gridSize[0],
            params.gridSize[1]);
    } else {
        crop_center_region3d<T1>(pGridData_crop_d, pGridData, params.imageSize[0],
            params.imageSize[1], params.imageSize[2],
            params.gridSize[0], params.gridSize[1],
            params.gridSize[2]);
    }

    // deapodization
    if (Nz == 1) {
        deapodization2d<T1>(pGridData_crop_deAp, pGridData_crop_d, Nx, Ny,
            kernelWidth, beta, params.gridOS);
    } else {
        deapodization3d<T1>(pGridData_crop_deAp, pGridData_crop_d, Nx, Ny, Nz,
            kernelWidth, beta, params.gridOS);
    }

#pragma acc exit data delete (samples [0:n])
#pragma acc update self(pGridData_crop_deAp [0:2 * imageNumElems])
    free(samples);
}

// Calculates the gridded forward fourier transform
template <typename T1>
void Gnufft<T1>::computeFd_CPU_Grid(int numK_per_coil, const T1* __restrict kx,
    const T1* __restrict ky, const T1* __restrict kz,
    const T1* __restrict dIn,
    int Nx, int Ny, int Nz, T1 gridOS,
    const T1 kernelWidth, const T1 beta, const T1* LUT,
    const uword sizeLUT, void* stream, CFTHandle* plan,
    T1* pGridData, T1* pGridData_d, T1* pGridData_os,
    T1* pGridData_os_d, T1* pSamples) const
{

    RANGE(__PRETTY_FUNCTION__)
    parameters<T1> params;
    params.sync = 0;
    params.binsize = 128;

    params.useLUT = 1;
    params.kernelWidth = kernelWidth;
    params.gridOS = gridOS;
    params.imageSize[0] = Nx; // gridSize is gridOS times larger than imageSize.
    params.imageSize[1] = Ny;
    params.imageSize[2] = Nz;
    params.gridSize[0] = std::ceil(gridOS * (T1)Nx);
    params.gridSize[1] = std::ceil(gridOS * (T1)Ny);
    if (params.gridSize[0] % 2) // 3D case, gridOS is adjusted on the z dimension:
        params.gridSize[0] += 1; // That why we need to make sure here that the xy
    if (params.gridSize[1] % 2) // dimensions have even sizes.
        params.gridSize[1] += 1;
    params.gridSize[2] = (Nz == 1) ? Nz : (std::ceil(gridOS * (T1)Nz)); // 2D or 3D
    params.numSamples = numK_per_coil;

    unsigned int n = params.numSamples;

    ///*
    // grid_size in xy-axis has to be divisible-by-two:
    //       (required by the cropImageRegion)
    // grid_size in z-axis has to be divisible-by-four:
    //       (required by the function gridding_GPU_3D(.))

    if (1 == Nz) {
        // round grid size (xy-axis) to the next divisible-by-two.
        gridOS = (T1)2.0 * std::ceil((gridOS * (T1)Nx) / (T1)2.0) / (T1)Nx;
    } else {
        // round grid size (z-axis) to the next divisible-by-four.
        gridOS = (T1)4.0 * std::ceil((gridOS * (T1)Nz) / (T1)4.0) / (T1)Nz;
    }
    //

    int imageNumElems = params.imageSize[0] * params.imageSize[1] * params.imageSize[2];

    memcpy(pGridData, dIn, sizeof(T1) * 2 * imageNumElems);
#pragma acc update device(pGridData [0:2 * imageNumElems])

#pragma acc parallel loop present(pSamples [0:2 * n])
    for (int ii = 0; ii < 2 * n; ii++) {
        pSamples[ii] = (T1)0.0;
    }

    // deapodization
    if (Nz == 1) {
        deapodization2d<T1>(pGridData_d, pGridData, Nx, Ny, kernelWidth, beta,
            params.gridOS);
    } else {
        deapodization3d<T1>(pGridData_d, pGridData, Nx, Ny, Nz, kernelWidth, beta,
            params.gridOS);
    }

    // zero pad
    if (Nz == 1) {
        zero_pad2d<T1>(pGridData_os, pGridData_d, Nx, Ny, params.gridOS);
    } else {

        zero_pad3d<T1>(pGridData_os, pGridData_d, Nx, Ny, Nz, params.gridOS);
    }

    if (Nz == 1) {
        fftshift2<T1>(pGridData_os_d, pGridData_os, params.gridSize[0],
            params.gridSize[1]);
    } else {
        fftshift3<T1>(pGridData_os_d, pGridData_os, params.gridSize[0],
            params.gridSize[1], params.gridSize[2]);
    }

// ifftn(gridData)
//#ifdef _OPENACC // We're on GPU 
#ifdef OPENACC_GPU
    // Inside this region the device data pointer will be used
    // cout << "about to reach openacc region in forward transform" << endl;

#pragma acc host_data use_device(pGridData_os_d)
    {

        // Query OpenACC for CUDA stream
        //void *stream = acc_get_cuda_stream(acc_async_sync);
        acc_set_cuda_stream(1, stream);
        // Launch FFT on the GPU
        if (Nz == 1) {
            fft2dGPU(pGridData_os_d, params.gridSize[0], params.gridSize[1], stream, plan);
            //fft2dCPU(pGridData_os_d, params.gridSize[0], params.gridSize[1]);

        } else {
            fft3dGPU(pGridData_os_d, params.gridSize[0], params.gridSize[1],
                params.gridSize[2], stream, plan);
            //fft3dCPU(pGridData_os_d, params.gridSize[0], params.gridSize[1],
            //params.gridSize[2]);
        }
    }
#else // We're on CPU
    if (Nz == 1) {
        fft2dCPU(pGridData_os_d, params.gridSize[0], params.gridSize[1]);
    } else {
        fft3dCPU(pGridData_os_d, params.gridSize[0], params.gridSize[1],
            params.gridSize[2]);
    }
#endif
    // ifftshift(gridData):
    //#pragma acc update device(pGridData_os_d[0:2*gridNumElems])
    if (Nz == 1) {
        ifftshift2<T1>(pGridData_os, pGridData_os_d, params.gridSize[0],
            params.gridSize[1]);
    } else {
        ifftshift3<T1>(pGridData_os, pGridData_os_d, params.gridSize[0],
            params.gridSize[1], params.gridSize[2]);
    }

    // Gridding with CPU - forward
    if (Nz == 1) {
        gridding_forward_2D(n, params, kx, ky, beta, pSamples, LUT, sizeLUT,
            pGridData_os);
    } else {
        gridding_forward_3D(n, params, kx, ky, kz, beta, pSamples, LUT, sizeLUT,
            pGridData_os);
    }

    // deallocate samples

#pragma acc update host(pSamples [0:2 * n])
    /*
	for (int ii = 0; ii < n; ii++) {
		outR_d[ii] = pSamples[2*ii];
		outI_d[ii] = pSamples[2*ii+1];
	}
*/
}



// Explicit Instantiation
template class Gnufft<float>;
template class Gnufft<double>;
