//
//  gridding.hpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 4/8/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef PowerGrid_gridding_hpp
#define PowerGrid_gridding_hpp

#include <cstdlib>

#ifdef _OPENACC //GPU Version
#include "openacc.h"
#include "accelmath.h"
#include "fftGPU.hpp"
#include "cufft.h"

#define MIN(a,b) std::min(a,b)
#define MAX(a,b) std::max(a,b)
#define SQRT(a) std::sqrt(a)
#define CEIL(a) std::ceil(a)
#define FLOOR(a) std::floor(a)
#define ABS(a) std::abs(a)
#else //CPU version

#include "fftCPU.hpp"

#define MIN(a, b) std::min(a,b)
#define MAX(a, b) std::max(a,b)
#define SQRT(a) std::sqrt(a)
#define CEIL(a) std::ceil(a)
#define FLOOR(a) std::floor(a)
#define ABS(a) std::abs(a)
#endif

using namespace arma;

// 2D adjoint gridding on CPU
template<typename T1>
int
gridding_Gold_2D(unsigned int n, parameters <T1> params, T1 beta, ReconstructionSample <T1>* __restrict sample,
		const T1* LUT, const uword sizeLUT,
		complex <T1>* __restrict gridData)
{

	unsigned int NxL, NxH;
	unsigned int NyL, NyH;
	//unsigned int NzL, NzH;

	unsigned int nx;
	unsigned int ny;
	//unsigned int nz;

	int idx;

	T1 w;

	T1 shiftedKx, shiftedKy/*, shiftedKz*/;
	T1 distX, kbX, distY, kbY/*, distZ, kbZ*/;
	T1* __restrict pGData;

	T1 kernelWidth = params.kernelWidth;
	//T1 beta = 18.5547;
	T1 gridOS = params.gridOS;

	unsigned int Nx = params.imageSize[0];
	unsigned int Ny = params.imageSize[1];
	//unsigned int Nz = params.imageSize[2];
	int gridNumElems = params.gridSize[0]*params.gridSize[1];

	pGData = reinterpret_cast<T1*>(gridData);
	//Jiading GAI
	//float t0 = t[0];

#pragma acc parallel loop gang vector pcopyin(LUT[0:sizeLUT]) pcopy(pGData[0:gridNumElems*2])
	for (int i = 0; i<n; i++) {
		ReconstructionSample <T1> pt = sample[i];

		//Jiading GAI
		//float atm = hanning_d(t[i], tau, l, t0);//a_l(t_m)

		shiftedKx = (gridOS)*(pt.kX+((T1) Nx)/(T1) 2.0);
		shiftedKy = (gridOS)*(pt.kY+((T1) Ny)/(T1) 2.0);
		//shiftedKz = ((float)gridOS)*(pt.kZ+((float)Nz)/2);

		//if(shiftedKx < 0.0f)
		//   shiftedKx = 0.0f;
		//if(shiftedKx > ((float)gridOS)*((float)Nx));
		//   shiftedKx = ((float)gridOS)*((float)Nx);
		//if(shiftedKy < 0.0f)
		//   shiftedKy = 0.0f;
		//if(shiftedKy > ((float)gridOS)*((float)Ny));
		//   shiftedKy = ((float)gridOS)*((float)Ny);


		NxL = (int) (MAX((T1) 0.0, CEIL(shiftedKx-kernelWidth*(gridOS)/(T1) 2.0)));
		NxH = (int) (MIN((gridOS*(T1) Nx-(T1) 1.0), FLOOR(shiftedKx+kernelWidth*(gridOS)/(T1) 2.0)));

		NyL = (int) (MAX((T1) 0.0, CEIL(shiftedKy-kernelWidth*(gridOS)/(T1) 2.0)));
		NyH = (int) (MIN((gridOS*(T1) Ny-(T1) 1.0), FLOOR(shiftedKy+kernelWidth*(gridOS)/(T1) 2.0)));

		//NzL = (int)(fmax(0.0f,ceil(shiftedKz - kernelWidth*((float)gridOS)/2)));
		//NzH = (int)(fmin((float)(gridOS*Nz-1),floor(shiftedKz + kernelWidth*((float)gridOS)/2)));

#pragma acc loop independent seq
		for (nx = NxL; nx<=NxH; ++nx) {
			int k0;
			distX = ABS(shiftedKx-((T1) nx))/(gridOS);
			if (params.useLUT) {
				k0 = (int) ((distX*distX*(T1) 4.0/(kernelWidth*kernelWidth))*(T1) sizeLUT);
				if (k0>=sizeLUT)
					kbY = (T1) 0.0;
				else
					kbY = LUT[k0];
				kbX = kernel_value_LUT(distX, LUT, sizeLUT, kernelWidth);
				//cout << "KBX = " << kbX << endl;

			}
			else {
				kbX = bessi0(beta*SQRT((T1) 1.0-((T1) 2.0*distX/kernelWidth)*((T1) 2.0*distX/kernelWidth)))/
						kernelWidth;

			}

			if (kbX!=kbX) {//if kbX = NaN
				kbX = 0;
			}
#pragma acc loop seq

			for (ny = NyL; ny<=NyH; ++ny) {
				distY = ABS(shiftedKy-((T1) ny))/(gridOS);
				if (params.useLUT) {

					k0 = (int) ((distY*distY*(T1) 4.0/(kernelWidth*kernelWidth))*(T1) sizeLUT);

					if (k0>=sizeLUT)
						kbY = (T1) 0.0;
					else
						kbY = LUT[k0];

				}
				else {
					kbY = bessi0(beta*SQRT((T1) 1.0-((T1) 2.0*distY/kernelWidth)*((T1) 2.0*distY/kernelWidth)))/
							kernelWidth;
				}

				if (kbY!=kbY) {//if kbY = NaN
					kbY = (T1) 0.0;
				}

				/* kernel weighting value */
				//if (params.useLUT){
				//    w = kbX * kbY;
				//} else

				w = kbX*kbY;

				/* grid data */
				idx = ny+(nx)*params.gridSize[1]/* + (nz)*gridOS*Nx*gridOS*Ny*/;

#pragma acc atomic update
				pGData[2*idx] += w*pt.real;
				// atomicAdd(pGData+2*idx, w*pt.real);

#pragma acc atomic update
				pGData[2*idx+1] += w*pt.imag;
				// atomicAdd(pGData+2*idx+1, w*pt.imag);

				//gridData[idx].y += (w*pt.imag*atm);
				//gridData[idx].x += (w*pt.real*atm);
				//gridData[idx].y += (w*pt.imag*atm);
				//gridData[idx].real(gridData[idx].real()+w*pt.real);
				//gridData[idx].imag(gridData[idx].imag()+w*pt.imag);
				/* estimate sample density */
				//#pragma acc atomic update
				//sampleDensity[idx] += w;
				//atomicAdd(sampleDensity+idx, w);
			}
		}
	}

	// re-arrange dimensions and output
	// Nady uses: x->y->z
	// IMPATIENT uses: z->x->y
	// PowerGrid uses: y->x->z because we are column major same as MATLAB...
	// Nope! XX So we need to convert from (x->y->z)-order to (z->x->y)-order
	//int gridNumElems = params.gridSize[0] * params.gridSize[1];

	//complex<T1> *gridData_reorder = (complex<T1>*) malloc(gridNumElems, sizeof(typename complex<T1>));
	//
	//for(int x=0;x<params.gridSize[0];x++)
	//    for(int y=0;y<params.gridSize[1];y++)
	//    {
	//        int lindex_nady      = x + y*params.gridSize[0];
	//        int lindex_impatient = y + x*params.gridSize[0];
	//
	//        gridData_reorder[lindex_impatient] = gridData[lindex_nady];
	//    }
	//memcpy((void*)gridData,(void*)gridData_reorder,gridNumElems*sizeof(typename complex<T1>));
	//
	//free(gridData_reorder);

	return 1;
}

// 3D adjoint gridding on CPU
template<typename T1>
int
gridding_Gold_3D(unsigned int n, parameters <T1> params, T1 beta, ReconstructionSample <T1>* __restrict sample,
		const T1* LUT, const uword sizeLUT,
		complex <T1>* gridData)
{
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
	T1* pGData;

	T1 kernelWidth = params.kernelWidth;
	//T1 beta = 18.5547;
	T1 gridOS = params.gridOS;

	int Nx = params.imageSize[0];
	int Ny = params.imageSize[1];
	int Nz = params.imageSize[2];
	int gridNumElems = params.gridSize[0]*params.gridSize[1]*params.gridSize[2];
	//Jiading GAI
	//float t0 = t[0];
	pGData = reinterpret_cast<T1*>(gridData);

#pragma acc parallel loop gang vector pcopyin(LUT[0:sizeLUT]) pcopy(pGData[0:gridNumElems*2])
	for (int i = 0; i<n; i++) {
		ReconstructionSample <T1> pt = sample[i];

		//Jiading GAI
		//float atm = hanning_d(t[i], tau, l, t0);//a_l(t_m)

		shiftedKx = (gridOS)*(pt.kX+((T1) Nx)/(T1) 2.0);
		shiftedKy = (gridOS)*(pt.kY+((T1) Ny)/(T1) 2.0);
		shiftedKz = (gridOS)*(pt.kZ+((T1) Nz)/(T1) 2.0);

		//	if(shiftedKx < 0.0f)
		//	   shiftedKx = 0.0f;
		//	if(shiftedKx > ((float)gridOS)*((float)Nx));
		//	   shiftedKx = ((float)gridOS)*((float)Nx);
		//	if(shiftedKy < 0.0f)
		//	   shiftedKy = 0.0f;
		//	if(shiftedKy > ((float)gridOS)*((float)Ny));
		//	   shiftedKy = ((float)gridOS)*((float)Ny);
		//	if(shiftedKz < 0.0f)
		//	   shiftedKz = 0.0f;
		//	if(shiftedKz > ((float)gridOS)*((float)Nz));
		//	   shiftedKz = ((float)gridOS)*((float)Nz);

		NxL = (int) (MAX((T1) 0.0, CEIL(shiftedKx-kernelWidth*(gridOS)/(T1) 2.0)));
		NxH = (int) (MIN((gridOS*(T1) Nx-(T1) 1.0), FLOOR(shiftedKx+kernelWidth*(gridOS)/(T1) 2.0)));

		NyL = (int) (MAX((T1) 0.0, CEIL(shiftedKy-kernelWidth*(gridOS)/(T1) 2.0)));
		NyH = (int) (MIN((gridOS*(T1) Ny-(T1) 1.0), FLOOR(shiftedKy+kernelWidth*(gridOS)/(T1) 2.0)));

		NzL = (int) (MAX((T1) 0.0, CEIL(shiftedKz-kernelWidth*(gridOS)/(T1) 2.0)));
		NzH = (int) (MIN((gridOS*(T1) Nz-(T1) 1.0), FLOOR(shiftedKz+kernelWidth*(gridOS)/(T1) 2.0)));
#pragma acc loop independent seq
		for (nz = NzL; nz<=NzH; ++nz) {
			int k0;
			distZ = ABS(shiftedKz-((T1) nz))/(gridOS);
			if (params.useLUT) {
				// kbZ = kernel_value_LUT(distZ, LUT, sizeLUT, kernelWidth);
				k0 = (int) ((distZ*distZ*(T1) 4.0/(kernelWidth*kernelWidth))*(T1) sizeLUT);
				if (k0>=sizeLUT)
					kbZ = (T1) 0.0;
				else
					kbZ = LUT[k0];
			}
			else {
				kbZ = bessi0(beta*SQRT((T1) 1.0-((T1) 2.0*distZ/kernelWidth)*
						((T1) 2.0*distZ/kernelWidth)))/kernelWidth;
			}

#pragma acc loop seq
			for (nx = NxL; nx<=NxH; ++nx) {
				distX = ABS(shiftedKx-((T1) nx))/(gridOS);
				if (params.useLUT) {
//                    kbX = kernel_value_LUT(distX, LUT, sizeLUT, kernelWidth);
					k0 = (int) ((distX*distX*(T1) 4.0/(kernelWidth*kernelWidth))*(T1) sizeLUT);
					if (k0>=sizeLUT)
						kbX = (T1) 0.0;
					else
						kbX = LUT[k0];
				}
				else {
					kbX = bessi0(beta*SQRT((T1) 1.0-((T1) 2.0*distX/kernelWidth)*
							((T1) 2.0*distX/kernelWidth)))/kernelWidth;
				}

#pragma acc loop seq
				for (ny = NyL; ny<=NyH; ++ny) {
					distY = ABS(shiftedKy-((T1) ny))/(gridOS);
					if (params.useLUT) {
//                      kbY = kernel_value_LUT(distY, LUT, sizeLUT, kernelWidth);
						k0 = (int) ((distY*distY*(T1) 4.0/(kernelWidth*kernelWidth))*(T1) sizeLUT);
						if (k0>=sizeLUT)
							kbY = (T1) 0.0;
						else
							kbY = LUT[k0];
					}
					else {
						kbY = bessi0(beta*SQRT((T1) 1.0-((T1) 2.0*distY/kernelWidth)*
								((T1) 2.0*distY/kernelWidth)))/kernelWidth;
					}

					w = kbX*kbY*kbZ;

					/* grid data */
					idx = ny+(nx)*params.gridSize[1]+(nz)*params.gridSize[0]*params.gridSize[1];
#pragma acc atomic update
					pGData[2*idx] += w*pt.real;

#pragma acc atomic update
					pGData[2*idx+1] += w*pt.imag;

					//gridData[idx].x += (w*pt.real*atm);
					//gridData[idx].y += (w*pt.imag*atm);
					//gridData[idx].real(gridData[idx].real()+w*pt.real);
					//gridData[idx].imag(gridData[idx].imag()+w*pt.imag);

					/* estimate sample density */
					//#pragma acc atomic update
					//sampleDensity[idx] += w;
				}
			}
		}
	}

	// re-arrange dimensions and output
	// Nady uses: x->y->z
	// IMPATIENT uses: z->x->y
	// So we need to convert from (x->y->z)-order to (z->x->y)-order
	/*
	int gridNumElems = params.gridSize[0] * params.gridSize[1] * params.gridSize[2];
	complex<T1> *gridData_reorder = (complex<T1>*) malloc(gridNumElems, sizeof(typename complex<T1>));

	for(int x=0;x<params.gridSize[0];x++)
		for(int y=0;y<params.gridSize[1];y++)
			for(int z=0;z<params.gridSize[2];z++)
			{
				int lindex_nady = x + y*params.gridSize[0] + z*params.gridSize[0]*params.gridSize[1];
				int lindex_impatient = z + x*params.gridSize[2] + y*params.gridSize[0]*params.gridSize[2];

				gridData_reorder[lindex_impatient] = gridData[lindex_nady];
			}
	memcpy((void*)gridData,(void*)gridData_reorder,gridNumElems*sizeof(typename complex<T1>));

	free(gridData_reorder);
	*/
	return 1;
}

// 2D forward gridding on CPU
template<typename T1>
int
gridding_Silver_2D(unsigned int n, parameters <T1> params, const T1* kx, const T1* ky, T1 beta,
		complex <T1>* __restrict sample,
		const T1* LUT, const uword sizeLUT,
		complex <T1>* __restrict gridData)
{

	int NxL, NxH;
	int NyL, NyH;
	//unsigned int NzL, NzH;

	int nx;
	int ny;
	//unsigned int nz;

	int idx;
	T1* pSamples;
	T1* pGridData;
	T1 w;
	T1 sampleReal;
	T1 sampleImag;
	T1 shiftedKx, shiftedKy/*, shiftedKz*/;
	T1 distX, kbX, distY, kbY/*, distZ, kbZ*/;

	T1 kernelWidth = params.kernelWidth;
	//T1 beta = 18.5547;
	T1 gridOS = params.gridOS;

	int Nx = params.imageSize[0];
	int Ny = params.imageSize[1];
	int imageNumElems = params.imageSize[0]*params.imageSize[1];
	int gridNumElems = params.gridSize[0]*params.gridSize[1];

	//unsigned int Nz = params.imageSize[2];
	pSamples = reinterpret_cast<T1*>(sample);
	pGridData = reinterpret_cast<T1*>(gridData);
	//Jiading GAI
	//float t0 = t[0];

#pragma acc parallel loop gang vector pcopyin(LUT[0:sizeLUT])  pcopy(pSamples[0:n*2]) pcopy(pGridData[0:gridNumElems*2])
	for (int i = 0; i<n; i++) {
		//complex<T1> pt = sample[i];

		//Jiading GAI
		//float atm = hanning_d(t[i], tau, l, t0);//a_l(t_m)

		shiftedKx = (gridOS)*(kx[i]+((T1) Nx)/(T1) 2.0);
		shiftedKy = (gridOS)*(ky[i]+((T1) Ny)/(T1) 2.0);

		//shiftedKz = ((float)gridOS)*(pt.kZ+((float)Nz)/2);

		//if(shiftedKx < 0.0f)
		//   shiftedKx = 0.0f;
		//if(shiftedKx > ((float)gridOS)*((float)Nx));
		//   shiftedKx = ((float)gridOS)*((float)Nx);
		//if(shiftedKy < 0.0f)
		//   shiftedKy = 0.0f;
		//if(shiftedKy > ((float)gridOS)*((float)Ny));
		//   shiftedKy = ((float)gridOS)*((float)Ny);


		NxL = (int) (MAX((T1) 0.0, CEIL(shiftedKx-kernelWidth*(gridOS)/(T1) 2.0)));
		NxH = (int) (MIN((gridOS*(T1) Nx-(T1) 1.0), FLOOR(shiftedKx+kernelWidth*(gridOS)/(T1) 2.0)));

		NyL = (int) (MAX((T1) 0.0, CEIL(shiftedKy-kernelWidth*(gridOS)/(T1) 2.0)));
		NyH = (int) (MIN((gridOS*(T1) Ny-(T1) 1.0), FLOOR(shiftedKy+kernelWidth*(gridOS)/(T1) 2.0)));

		//NzL = (int)(fmax(0.0f,ceil(shiftedKz - kernelWidth*((float)gridOS)/2)));
		//NzH = (int)(fmin((float)(gridOS*Nz-1),floor(shiftedKz + kernelWidth*((float)gridOS)/2)));

#pragma acc loop independent seq
		for (nx = NxL; nx<=NxH; ++nx) {
			int k0;
			distX = ABS(shiftedKx-((T1) nx))/(gridOS);
			if (params.useLUT) {
				//kbX = kernel_value_LUT(distX, LUT, sizeLUT, kernelWidth);

				k0 = (int) ((distX*distX*(T1) 4.0/(kernelWidth*kernelWidth))*(T1) sizeLUT);
				if (k0>=sizeLUT)
					kbX = (T1) 0.0;
				else
					kbX = LUT[k0];

			}
			else {
				kbX = bessi0(beta*SQRT((T1) 1.0-((T1) 2.0*distX/kernelWidth)*((T1) 2.0*distX/kernelWidth)))/kernelWidth;

			}

			if (kbX!=kbX) {//if kbX = NaN
				kbX = 0;
			}

#pragma acc loop seq
			for (ny = NyL; ny<=NyH; ++ny) {
				distY = ABS(shiftedKy-((T1) ny))/(gridOS);
				if (params.useLUT) {
					//kbY = kernel_value_LUT(distY, LUT, sizeLUT, kernelWidth);

					k0 = (int) ((distY*distY*(T1) 4.0/(kernelWidth*kernelWidth))*(T1) sizeLUT);
					if (k0>=sizeLUT)
						kbY = (T1) 0.0;
					else
						kbY = LUT[k0];

				}
				else {
					kbY = bessi0(beta*SQRT((T1) 1.0-((T1) 2.0*distY/kernelWidth)*((T1) 2.0*distY/kernelWidth)))
							/kernelWidth;
				}

				if (kbY!=kbY) {//if kbY = NaN
					kbY = 0;
				}

				/* kernel weighting value */
				//if (params.useLUT){
				//    w = kbX * kbY;
				//} else {
				w = kbX*kbY;
				//}
				/* grid data */
				idx = ny+(nx)*params.gridSize[1]/* + (nz)*gridOS*Nx*gridOS*Ny*/;
				//gridData[idx].x += (w*pt.real*atm);
				//gridData[idx].y += (w*pt.imag*atm);

#pragma acc atomic update
				pSamples[2*i] += w*pGridData[2*idx];

#pragma acc atomic update
				pSamples[2*i+1] += w*pGridData[2*idx+1];

				//sample[i].real(sample[i].real()+w*gridData[idx].real());
				//sample[i].imag(sample[i].imag()+w*gridData[idx].imag());
				/* estimate sample density */

				//#pragma acc atomic update
				//sampleDensity[i] += w;
			}
		}
	}

	// re-arrange dimensions and output
	// Nady uses: x->y->z
	// IMPATIENT uses: z->x->y
	// PowerGrid uses: x->y->z because we are column major same as Nady...
	// Nope! XX So we need to convert from (x->y->z)-order to (z->x->y)-order
	//int gridNumElems = params.gridSize[0] * params.gridSize[1];

	//complex<T1> *gridData_reorder = (complex<T1>*) malloc(gridNumElems, sizeof(typename complex<T1>));
	//
	//for(int x=0;x<params.gridSize[0];x++)
	//    for(int y=0;y<params.gridSize[1];y++)
	//    {
	//        int lindex_nady      = x + y*params.gridSize[0];
	//        int lindex_impatient = y + x*params.gridSize[0];
	//
	//        gridData_reorder[lindex_impatient] = gridData[lindex_nady];
	//    }
	//memcpy((void*)gridData,(void*)gridData_reorder,gridNumElems*sizeof(typename complex<T1>));
	//
	//free(gridData_reorder);

	return 1;
}

// 3D forward gridding on CPU
template<typename T1>
int
gridding_Silver_3D(unsigned int n, parameters <T1> params, const T1* kx, const T1* ky, const T1* kz, T1 beta,
		complex <T1>* __restrict sample,
		const T1* LUT, const uword sizeLUT,
		complex <T1>* __restrict gridData)
{
	int NxL, NxH;
	int NyL, NyH;
	int NzL, NzH;

	int nx;
	int ny;
	int nz;

	int idx;
	T1* pSamples;
	T1* pGridData;
	T1 w;

	T1 shiftedKx, shiftedKy, shiftedKz;
	T1 distX, kbX, distY, kbY, distZ, kbZ;

	T1 kernelWidth = params.kernelWidth;
	//T1 beta = 18.5547;
	T1 gridOS = params.gridOS;

	int Nx = params.imageSize[0];
	int Ny = params.imageSize[1];
	int Nz = params.imageSize[2];
	int imageNumElems = params.imageSize[0]*params.imageSize[1]*params.imageSize[2];
	int gridNumElems = params.gridSize[0]*params.gridSize[1]*params.gridSize[2];
	pSamples = reinterpret_cast<T1*>(sample);
	pGridData = reinterpret_cast<T1*>(gridData);
	//Jiading GAI
	//float t0 = t[0];

#pragma acc parallel loop gang vector pcopyin(LUT[0:sizeLUT])  pcopy(pSamples[0:n*2]) pcopy(pGridData[0:gridNumElems*2])
	for (int i = 0; i<n; i++) {
		//complex<T1> pt = sample[i];

		//Jiading GAI
		//float atm = hanning_d(t[i], tau, l, t0);//a_l(t_m)

		shiftedKx = (gridOS)*(kx[i]+((T1) Nx)/(T1) 2.0);
		shiftedKy = (gridOS)*(ky[i]+((T1) Ny)/(T1) 2.0);
		shiftedKz = (gridOS)*(kz[i]+((T1) Nz)/(T1) 2.0);

		//	if(shiftedKx < 0.0f)
		//	   shiftedKx = 0.0f;
		//	if(shiftedKx > ((float)gridOS)*((float)Nx));
		//	   shiftedKx = ((float)gridOS)*((float)Nx);
		//	if(shiftedKy < 0.0f)
		//	   shiftedKy = 0.0f;
		//	if(shiftedKy > ((float)gridOS)*((float)Ny));
		//	   shiftedKy = ((float)gridOS)*((float)Ny);
		//	if(shiftedKz < 0.0f)
		//	   shiftedKz = 0.0f;
		//	if(shiftedKz > ((float)gridOS)*((float)Nz));
		//	   shiftedKz = ((float)gridOS)*((float)Nz);


		NxL = (int) (MAX((T1) 0.0, CEIL(shiftedKx-kernelWidth*(gridOS)/(T1) 2.0)));
		NxH = (int) (MIN((gridOS*(T1) Nx-(T1) 1.0), FLOOR(shiftedKx+kernelWidth*((T1) gridOS)/(T1) 2.0)));

		NyL = (int) (MAX((T1) 0.0, CEIL(shiftedKy-kernelWidth*(gridOS)/(T1) 2.0)));
		NyH = (int) (MIN((gridOS*(T1) Ny-(T1) 1.0), FLOOR(shiftedKy+kernelWidth*((T1) gridOS)/(T1) 2.0)));

		NzL = (int) (MAX((T1) 0.0, CEIL(shiftedKz-kernelWidth*(gridOS)/2.0f)));
		NzH = (int) (MIN((gridOS*(T1) Nz-(T1) 1.0), FLOOR(shiftedKz+kernelWidth*((T1) gridOS)/(T1) 2.0)));

#pragma acc loop independent seq
		for (nz = NzL; nz<=NzH; ++nz) {
			int k0;
			distZ = ABS(shiftedKz-((T1) nz))/(gridOS);
			if (params.useLUT) {
				//kbZ = kernel_value_LUT(distZ, LUT, sizeLUT, kernelWidth);
				k0 = (int) ((distZ*distZ*(T1) 4.0/(kernelWidth*kernelWidth))*(T1) sizeLUT);
				if (k0>=sizeLUT)
					kbZ = (T1) 0.0;
				else
					kbZ = LUT[k0];
				//cout << "KBX = " << kbX << endl;

			}
			else {
				kbZ = bessi0(beta*SQRT((T1) 1.0-((T1) 2.0*distZ/kernelWidth)*((T1) 2.0*distZ/kernelWidth)))/
						kernelWidth;

			}

			if (kbZ!=kbZ) {//if kbY = NaN
				kbZ = 0;
			}

#pragma acc loop seq
			for (nx = NxL; nx<=NxH; ++nx) {
				distX = ABS(shiftedKx-((T1) nx))/(gridOS);
				if (params.useLUT) {
					//kbX = kernel_value_LUT(distX, LUT, sizeLUT, kernelWidth);
					k0 = (int) ((distX*distX*(T1) 4.0/(kernelWidth*kernelWidth))*(T1) sizeLUT);
					if (k0>=sizeLUT)
						kbX = (T1) 0.0;
					else
						kbX = LUT[k0];
					//cout << "KBX = " << kbX << endl;

				}
				else {
					kbX = bessi0(beta*SQRT((T1) 1.0-((T1) 2.0*distX/kernelWidth)*((T1) 2.0*distX/kernelWidth)))/
							kernelWidth;

				}

				if (kbX!=kbX) {//if kbX = NaN
					kbX = 0;
				}

#pragma acc loop seq
				for (ny = NyL; ny<=NyH; ++ny) {
					distY = ABS(shiftedKy-((T1) ny))/(gridOS);
					if (params.useLUT) {
						//kbY = kernel_value_LUT(distY, LUT, sizeLUT, kernelWidth);

						k0 = (int) ((distY*distY*(T1) 4.0/(kernelWidth*kernelWidth))*(T1) sizeLUT);
						if (k0>=sizeLUT)
							kbY = (T1) 0.0;
						else
							kbY = LUT[k0];

					}
					else {
						kbY = bessi0(
								beta*SQRT(1.0-(2.0*distY/kernelWidth)*(2.0*distY/kernelWidth)))/
								kernelWidth;
					}

					if (kbY!=kbY) {//if kbY = NaN
						kbY = 0;
					}
					/* kernel weighting value */
					//if (params.useLUT){
					//    w = kbX * kbY;
					//} else {

					w = kbX*kbY*kbZ;

					/* grid data */
					idx = ny+(nx)*params.gridSize[1]+(nz)*params.gridSize[0]*params.gridSize[1];
					//gridData[idx].x += (w*pt.real*atm);
					//gridData[idx].y += (w*pt.imag*atm);
#pragma acc atomic update
					pSamples[2*i] += w*pGridData[2*idx];

#pragma acc atomic update
					pSamples[2*i+1] += w*pGridData[2*idx+1];
					//sample[i].real(sample[i].real()+w*gridData[idx].real());
					//sample[i].imag(sample[i].imag()+w*gridData[idx].imag());


					/* estimate sample density */

					//#pragma acc atomic update
					//sampleDensity[i] += w;
				}
			}
		}
	}

	// re-arrange dimensions and output
	// Nady uses: x->y->z
	// IMPATIENT uses: z->x->y
	// So we need to convert from (x->y->z)-order to (z->x->y)-order
	/*
	int gridNumElems = params.gridSize[0] * params.gridSize[1] * params.gridSize[2];
	complex<T1> *gridData_reorder = (complex<T1>*) malloc(gridNumElems, sizeof(typename complex<T1>));

	for(int x=0;x<params.gridSize[0];x++)
		for(int y=0;y<params.gridSize[1];y++)
			for(int z=0;z<params.gridSize[2];z++)
			{
				int lindex_nady = x + y*params.gridSize[0] + z*params.gridSize[0]*params.gridSize[1];
				int lindex_impatient = z + x*params.gridSize[2] + y*params.gridSize[0]*params.gridSize[2];

				gridData_reorder[lindex_impatient] = gridData[lindex_nady];
			}
	memcpy((void*)gridData,(void*)gridData_reorder,gridNumElems*sizeof(typename complex<T1>));

	free(gridData_reorder);
	*/
	return 1;
}

// Calculates the gridded adjoint transform
template<typename T1>
void
computeFH_CPU_Grid(
		int numK_per_coil, const T1* __restrict kx, const T1* __restrict ky, const T1* __restrict kz,
		const T1* __restrict dR, const T1* __restrict dI, int Nx, int Ny, int Nz,
		T1 gridOS, T1* __restrict outR_d, T1* __restrict outI_d, const T1 kernelWidth, const T1 beta,
		const T1* LUT, const uword sizeLUT)
{

	/*
	 *  Based on Eqn. (5) of Beatty's gridding paper:
	 *  "Rapid Gridding Reconstruction With a Minimal Oversampling Ratio"
	 *
	 *  Note that Beatty use their kernel width to be equal twice the window
	 *  width parameter used by Jackson et al.
	 */
	/*
	T1 kernelWidth = 4.0;
	T1 beta = MRI_PI * std::sqrt( (gridOS - 0.5) * (gridOS - 0.5) *
						   (kernelWidth * kernelWidth*4.0) /
						   (gridOS * gridOS) - 0.8
						   );
	*/

	parameters <T1> params;
	params.sync = 0;
	params.binsize = 128;

	params.useLUT = 1;
	params.kernelWidth = kernelWidth;
	params.gridOS = gridOS;
	params.imageSize[0] = Nx;//gridSize is gridOS times larger than imageSize.
	params.imageSize[1] = Ny;
	params.imageSize[2] = Nz;
	params.gridSize[0] = CEIL(gridOS*(T1) Nx);
	params.gridSize[1] = CEIL(gridOS*(T1) Ny);
	if (params.gridSize[0]%2)//3D case, gridOS is adjusted on the z dimension:
		params.gridSize[0] += 1;//That why we need to make sure here that the xy
	if (params.gridSize[1]%2)//dimensions have even sizes.
		params.gridSize[1] += 1;
	params.gridSize[2] = (Nz==1) ? Nz : (CEIL(gridOS*(T1) Nz));// 2D or 3D
	params.numSamples = numK_per_coil;

	//T1 *sampleDensity;
	complex <T1>* gridData_d;

	ReconstructionSample <T1>* samples; //Input Data
	//allocate samples
	//cout << "Allocating samples" << endl;
	samples = (ReconstructionSample <T1>*) malloc(params.numSamples*sizeof(ReconstructionSample<T1>));

	if (samples==NULL) {
		printf("ERROR: Unable to allocate memory for input data\n");
		exit(1);
	}
	//cout << "Finished samples" << endl;
	unsigned int n = params.numSamples;
	//
	for (int i = 0; i<params.numSamples; i++) {
		if (ABS(kx[i])>(Nx/(T1) 2.0) ||
				ABS(ky[i])>(Ny/(T1) 2.0) ||
				ABS(kz[i])>(Nz/(T1) 2.0)
				) {

			printf("\nError:k-space trajectory out of range [-N/2,N/2]:\n      gridding requires that k-space should be contained within the window -N/2 to N/2.\n");
			cout << "kx = " << kx[i] << " ky = " << ky[i] << " kz = " << kz[i] << " i = " << i << endl;
			exit(1);
		}
		else {

			samples[i].kX = kx[i];
			samples[i].kY = ky[i];
			samples[i].kZ = kz[i];

			samples[i].real = dR[i];
			samples[i].imag = dI[i];

			samples[i].sdc = (T1) 1.0;
			//samples[i].t = t[i];
		}
	}    ///*
	// grid_size in xy-axis has to be divisible-by-two:
	//       (required by the cropImageRegion)
	// grid_size in z-axis has to be divisible-by-four:
	//       (required by the function gridding_GPU_3D(.))
	if (1==Nz) {
		//round grid size (xy-axis) to the next divisible-by-two.
		gridOS = (T1) 2.0*CEIL((gridOS*(T1) Nx)/(T1) 2.0)/(T1) Nx;
	}
	else {
		//round grid size (z-axis) to the next divisible-by-four.
		gridOS = (T1) 4.0*CEIL((gridOS*(T1) Nz)/(T1) 4.0)/(T1) Nz;
	}
	// */
	int gridNumElems = params.gridSize[0]*
			params.gridSize[1]*
			params.gridSize[2];

	int imageNumElems = params.imageSize[0]*
			params.imageSize[1]*
			params.imageSize[2];

	complex <T1>* gridData = new complex<T1>[gridNumElems];

	// Have to set 'gridData' and 'sampleDensity' to zero.
	// Because they will be involved in accumulative operations
	// inside gridding functions.
	//#pragma acc parallel loop
	for (int i = 0; i<gridNumElems; i++) {
		gridData[i].real((T1) 0.0);
		gridData[i].imag((T1) 0.0);
		//sampleDensity[i] = (T1)0.0;
	}
	gridData_d = new complex<T1>[gridNumElems];
	complex <T1>* gridData_crop_d = new complex<T1>[imageNumElems];
	complex <T1>* gridData_crop_deAp = new complex<T1>[imageNumElems];
	T1* pGridData_crop_d = reinterpret_cast<T1*>(gridData_crop_d);
	T1* pGridData_crop_deAp = reinterpret_cast<T1*>(gridData_crop_deAp);
	T1* pGridData_d = reinterpret_cast<T1*>(gridData_d);
	T1* pGridData = reinterpret_cast<T1*>(gridData);
#pragma acc enter data copyin(pGridData[0:2*gridNumElems]) create(pGridData_d[0:2*gridNumElems],pGridData_crop_d[0:2*imageNumElems],pGridData_crop_deAp[0:2*imageNumElems], outR_d[0:imageNumElems], outI_d[0:imageNumElems])
	// Gridding with CPU - gold
	if (Nz==1) {
		gridding_Gold_2D<T1>(n, params, beta, samples, LUT, sizeLUT,
				gridData);
	}
	else {
		gridding_Gold_3D<T1>(n, params, beta, samples, LUT, sizeLUT,
				gridData);
	}

	if (Nz==1) {
		ifftshift2<T1>(pGridData_d, pGridData, params.gridSize[0], params.gridSize[1]);
	}
	else {
		ifftshift3<T1>(pGridData_d, pGridData, params.gridSize[0], params.gridSize[1], params.gridSize[2]);
	}


#ifdef _OPENACC // We're on GPU
	// Inside this region the device data pointer will be used for cuFFT

//#pragma acc data copy(pGridData_d[0:2*gridNumElems])
//	{
#pragma acc host_data use_device(pGridData_d)
	{
	   // Query OpenACC for CUDA stream
	   void *stream = acc_get_cuda_stream(acc_async_sync);

	   // Launch FFT on the GPU
	   if (Nz == 1) {
		   ifft2dGPU(pGridData_d, params.gridSize[0], params.gridSize[1], stream);
	   } else {
		   ifft3dGPU(pGridData_d, params.gridSize[0], params.gridSize[1], params.gridSize[2], stream);
	   }
	}


#else // We're on CPU so we'll use FFTW

	// Launch FFT on the GPU
	if (Nz==1) {
		ifft2dCPU(pGridData_d, params.gridSize[0], params.gridSize[1]);
	}
	else {
		ifft3dCPU(pGridData_d, params.gridSize[0], params.gridSize[1], params.gridSize[2]);
	}

#endif
	//cout << "Got through the update device directive" << endl;
	if (Nz==1) {
		fftshift2<T1>(pGridData, pGridData_d, params.gridSize[0],
				params.gridSize[1]);
	}
	else {
		fftshift3<T1>(pGridData, pGridData_d, params.gridSize[0],
				params.gridSize[1], params.gridSize[2]);
	}


	if (Nz==1) {
		crop_center_region2d<T1>(pGridData_crop_d, pGridData,
				params.imageSize[0], params.imageSize[1],
				params.gridSize[0], params.gridSize[1]);
	}
	else {
		crop_center_region3d<T1>(pGridData_crop_d, pGridData,
				params.imageSize[0], params.imageSize[1], params.imageSize[2],
				params.gridSize[0], params.gridSize[1], params.gridSize[2]);
	}
	// deapodization
	if (Nz==1) {
		deapodization2d<T1>(pGridData_crop_deAp, pGridData_crop_d,
				Nx, Ny, kernelWidth, beta, params.gridOS);
	}
	else {
		deapodization3d<T1>(pGridData_crop_deAp, pGridData_crop_d,
				Nx, Ny, Nz, kernelWidth, beta, params.gridOS);
	}

	// Copy results from gridData_crop_d to outR_d and outI_d
	// gridData_crop_d is cufftComplex, interleaving
	// De-interleaving the data from cufftComplex to outR_d-and-outI_d

	if (Nz==1) {
		deinterleave_data2d<T1>(pGridData_crop_deAp, outR_d, outI_d, Nx, Ny);
	}
	else {
		deinterleave_data3d<T1>(pGridData_crop_deAp, outR_d, outI_d, Nx, Ny, Nz);
	}

#pragma acc exit data copyout(outR_d[0:imageNumElems],outI_d[0:imageNumElems]) delete(pGridData_crop_d,pGridData_d,pGridData, pGridData_crop_deAp)
	delete[] gridData_crop_d;
	delete[] gridData_crop_deAp;
	free(samples);
	delete[] gridData;
	delete[] gridData_d;
}

//Calculates the gridded forward fourier transform
template<typename T1>
void
computeFd_CPU_Grid(
		int numK_per_coil, const T1* __restrict kx, const T1* __restrict ky, const T1* __restrict kz,
		const T1* __restrict dR, const T1* __restrict dI, int Nx, int Ny, int Nz,
		T1 gridOS, T1* __restrict outR_d, T1* __restrict outI_d, const T1 kernelWidth, const T1 beta, const T1* LUT,
		const uword sizeLUT)
{

	/*
	 *  Based on Eqn. (5) of Beatty's gridding paper:
	 *  "Rapid Gridding Reconstruction With a Minimal Oversampling Ratio"
	 *
	 *  Note that Beatty use their kernel width to be equal twice the window
	 *  width parameter used by Jackson et al.
	 */
	/*
	T1 kernelWidth = .0;
	T1 beta = MRI_PI * std::sqrt( (gridOS - 0.5) * (gridOS - 0.5) *
								  (kernelWidth * kernelWidth*4.0) /
								  (gridOS * gridOS) - 0.8
	);i
	*/
	parameters <T1> params;
	params.sync = 0;
	params.binsize = 128;

	params.useLUT = 1;
	params.kernelWidth = kernelWidth;
	params.gridOS = gridOS;
	params.imageSize[0] = Nx;//gridSize is gridOS times larger than imageSize.
	params.imageSize[1] = Ny;
	params.imageSize[2] = Nz;
	params.gridSize[0] = CEIL(gridOS*(T1) Nx);
	params.gridSize[1] = CEIL(gridOS*(T1) Ny);
	if (params.gridSize[0]%2)//3D case, gridOS is adjusted on the z dimension:
		params.gridSize[0] += 1;//That why we need to make sure here that the xy
	if (params.gridSize[1]%2)//dimensions have even sizes.
		params.gridSize[1] += 1;
	params.gridSize[2] = (Nz==1) ? Nz : (CEIL(gridOS*(T1) Nz));// 2D or 3D
	params.numSamples = numK_per_coil;

	complex <T1>* samples = new complex<T1>[params.numSamples];


	if (samples==NULL) {
		printf("ERROR: Unable to allocate memory for input data\n");
		exit(1);
	}

	unsigned int n = params.numSamples;
	//
	for (int i = 0; i<params.numSamples; i++) {
		if (ABS(kx[i])>(Nx/(T1) 2.0) ||
				ABS(ky[i])>(Ny/(T1) 2.0) ||
				ABS(kz[i])>(Nz/(T1) 2.0)
				) {

			printf("\nError:k-space trajectory out of range [-N/2,N/2]:\n      gridding requires that k-space should be contained within the window -N/2 to N/2.\n");
			cout << "kx = " << kx[i] << " ky = " << ky[i] << " kz = " << kz[i] << " i = " << i << endl;
			exit(1);
		}
		else {

			samples[i].real((T1) 0.0);
			samples[i].imag((T1) 0.0);

		}
	}


	///*
	// grid_size in xy-axis has to be divisible-by-two:
	//       (required by the cropImageRegion)
	// grid_size in z-axis has to be divisible-by-four:
	//       (required by the function gridding_GPU_3D(.))

	if (1==Nz) {
		//round grid size (xy-axis) to the next divisible-by-two.
		gridOS = (T1) 2.0*CEIL((gridOS*(T1) Nx)/(T1) 2.0)/(T1) Nx;
	}
	else {
		//round grid size (z-axis) to the next divisible-by-four.
		gridOS = (T1) 4.0*CEIL((gridOS*(T1) Nz)/(T1) 4.0)/(T1) Nz;
	}
	//

	int gridNumElems = params.gridSize[0]*
			params.gridSize[1]*
			params.gridSize[2];
	int imageNumElems = params.imageSize[0]*
			params.imageSize[1]*
			params.imageSize[2];

	//allocate gridData
	complex <T1>* gridData = new complex<T1>[imageNumElems];
	// Have to set 'gridData' to zero.
	// Because they will be involved in accumulative operations
	// inside gridding functions.
	for (int i = 0; i<imageNumElems; i++) {
		gridData[i].real(dR[i]);
		gridData[i].imag(dI[i]);
	}
	complex <T1>* gridData_d = new complex<T1>[imageNumElems];
	complex <T1>* gridData_os_d = new complex<T1>[gridNumElems];
	complex <T1>* gridData_os = new complex<T1>[gridNumElems];
	T1* pGridData_d = reinterpret_cast<T1*>(gridData_d);
	T1* pGridData_os_d = reinterpret_cast<T1*>(gridData_os_d);
	T1* pGridData_os = reinterpret_cast<T1*>(gridData_os);
	T1* pGridData = reinterpret_cast<T1*>(gridData);

#pragma acc enter data copyin(pGridData[0:2*imageNumElems]) create(pGridData_d[0:2*imageNumElems],pGridData_os[0:2*gridNumElems],pGridData_os_d[0:2*gridNumElems])
	// deapodization
	if (Nz==1) {
		deapodization2d<T1>(pGridData_d, pGridData,
				Nx, Ny, kernelWidth, beta, params.gridOS);
	}
	else {
		deapodization3d<T1>(pGridData_d, pGridData,
				Nx, Ny, Nz, kernelWidth, beta, params.gridOS);
	}

	//zero pad

	if (Nz==1) {
		zero_pad2d<T1>(pGridData_os, pGridData_d,
				Nx, Ny, params.gridOS);
	}
	else {

		zero_pad3d<T1>(pGridData_os, pGridData_d,
				Nx, Ny, Nz, params.gridOS);

	}

	if (Nz==1) {
		fftshift2<T1>(pGridData_os_d, pGridData_os, params.gridSize[0],
				params.gridSize[1]);
	}
	else {
		fftshift3<T1>(pGridData_os_d, pGridData_os, params.gridSize[0],
				params.gridSize[1], params.gridSize[2]);
	}


	// ifftn(gridData)
#ifdef _OPENACC // We're on GPU
	// Inside this region the device data pointer will be used
	//cout << "about to reach openacc region in forward transform" << endl;

//#pragma acc data copy(pGridData_os_d[0:2*gridNumElems])
//{
#pragma acc host_data use_device(pGridData_os_d)
		{
	   // Query OpenACC for CUDA stream
			void* stream = acc_get_cuda_stream(acc_async_sync);

	   // Launch FFT on the GPU
			if (Nz==1) {
				fft2dGPU(pGridData_os_d, params.gridSize[0], params.gridSize[1], stream);
			}
			else {
				fft3dGPU(pGridData_os_d, params.gridSize[0], params.gridSize[1], params.gridSize[2], stream);
			}
		}
#else // We're on CPU
	if (Nz==1) {
		fft2dCPU(pGridData_os_d, params.gridSize[0], params.gridSize[1]);
	}
	else {
		fft3dCPU(pGridData_os_d, params.gridSize[0], params.gridSize[1], params.gridSize[2]);
	}
#endif
	// ifftshift(gridData):
	if (Nz==1) {
		ifftshift2<T1>(pGridData_os, pGridData_os_d, params.gridSize[0], params.gridSize[1]);
	}
	else {
		ifftshift3<T1>(pGridData_os, pGridData_os_d, params.gridSize[0], params.gridSize[1], params.gridSize[2]);
	}

	// Gridding with CPU - silver
	if (Nz==1) {
		gridding_Silver_2D<T1>(n, params, kx, ky, beta, samples, LUT, sizeLUT,
				gridData_os);
	}
	else {
		gridding_Silver_3D<T1>(n, params, kx, ky, kz, beta, samples, LUT, sizeLUT,
				gridData_os);
	}

	/*
	// crop the center region of the "image".
	if(Nz==1)
	{
		crop_center_region2d(gridData_crop_d, gridData,
							 params.imageSize[0], params.imageSize[1],
							 params.gridSize[0], params.gridSize[1]);
	}
	else {
		crop_center_region3d(gridData_crop_d, gridData,
							 params.imageSize[0], params.imageSize[1], params.imageSize[2],
							 params.gridSize[0], params.gridSize[1], params.gridSize[2]);
	}

	// Copy results from gridData_crop_d to outR_d and outI_d
	// gridData_crop_d is cufftComplex, interleaving
	// De-interleaving the data from cufftComplex to outR_d-and-outI_d
	if(Nz==1)
	{
		deinterleave_data2d(gridData_crop_d, outR_d, outI_d, Nx, Ny);
	}
	else
	{
		deinterleave_data3d(gridData_crop_d, outR_d, outI_d, Nx, Ny, Nz);
	}
	 */

	for (int ii = 0; ii<n; ii++) {
		outR_d[ii] = samples[ii].real();
		outI_d[ii] = samples[ii].imag();
	}
	//deallocate samples
#pragma acc exit data delete(pGridData_d[0:2*imageNumElems],pGridData_os[0:2*gridNumElems],pGridData_os_d[0:2*gridNumElems],pGridData[0:2*imageNumElems])
	delete[] samples;
	delete[] gridData;
	delete[] gridData_d;
	delete[] gridData_os;
	delete[] gridData_os_d;
	//delete[] sampleDensity;
	//delete(gridData_crop_d);
}

#endif
