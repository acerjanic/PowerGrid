//
//  gridding.hpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 4/8/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef PowerGrid_gridding_hpp
#define PowerGrid_gridding_hpp

#include "fftw3.h"

using namespace arma;


// 2D adjoint gridding on CPU
template<typename T1>
int
gridding_Gold_2D(unsigned int n, parameters<T1> params, T1 beta, ReconstructionSample<T1> *__restrict sample,
                 const T1 *LUT, const uword sizeLUT,
                 complex<T1> *__restrict gridData, T1 *__restrict sampleDensity)
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
    
    T1 kernelWidth = params.kernelWidth;
    //T1 beta = 18.5547;
    T1 gridOS = params.gridOS;
    
    unsigned int Nx = params.imageSize[0];
    unsigned int Ny = params.imageSize[1];
    //unsigned int Nz = params.imageSize[2];
    
    //Jiading GAI
    //float t0 = t[0];
    
    for (unsigned int i=0; i < n; i++)
    {
        ReconstructionSample<T1> pt = sample[i];
        
        //Jiading GAI
        //float atm = hanning_d(t[i], tau, l, t0);//a_l(t_m)
        
        shiftedKx = (gridOS)*(pt.kX+((T1)Nx)/2.0f);
        shiftedKy = (gridOS)*(pt.kY+((T1)Ny)/2.0f);
        //shiftedKz = ((float)gridOS)*(pt.kZ+((float)Nz)/2);
        
        //if(shiftedKx < 0.0f)
        //   shiftedKx = 0.0f;
        //if(shiftedKx > ((float)gridOS)*((float)Nx));
        //   shiftedKx = ((float)gridOS)*((float)Nx);
        //if(shiftedKy < 0.0f)
        //   shiftedKy = 0.0f;
        //if(shiftedKy > ((float)gridOS)*((float)Ny));
        //   shiftedKy = ((float)gridOS)*((float)Ny);
        
        
        NxL = (int)(fmax(0.0f,std::ceil(shiftedKx - kernelWidth*(gridOS)/2.0f)));
        NxH = (int)(fmin((gridOS*(T1)Nx-1.0f),std::floor(shiftedKx + kernelWidth*(gridOS)/2.0f)));
        
        NyL = (int)(fmax(0.0f,std::ceil(shiftedKy - kernelWidth*(gridOS)/2.0f)));
        NyH = (int)(fmin((gridOS*(T1)Ny-1.0f),std::floor(shiftedKy + kernelWidth*(gridOS)/2.0f)));
        
        //NzL = (int)(fmax(0.0f,ceil(shiftedKz - kernelWidth*((float)gridOS)/2)));
        //NzH = (int)(fmin((float)(gridOS*Nz-1),floor(shiftedKz + kernelWidth*((float)gridOS)/2)));


        for(nx=NxL; nx<=NxH; ++nx)
        {
            distX = fabs(shiftedKx - ((T1)nx))/(gridOS);
            if (params.useLUT) {
                kbX = kernel_value_LUT(distX, LUT, sizeLUT, kernelWidth);
                //cout << "KBX = " << kbX << endl;

            } else {
                kbX = bessi0(beta * std::sqrt(1.0 - (2.0 * distX / kernelWidth) * (2.0 * distX / kernelWidth))) /
                      kernelWidth;

            }

            if (kbX != kbX) {//if kbX = NaN
                kbX = 0;
            }
            
            for(ny=NyL; ny<=NyH; ++ny)
            {
                distY = fabs(shiftedKy - ((T1)ny))/(gridOS);
                if (params.useLUT){
                    kbY = kernel_value_LUT(distY, LUT, sizeLUT, kernelWidth);
                } else {
                    kbY = bessi0(beta * std::sqrt(1.0 - (2.0 * distY / kernelWidth) * (2.0 * distY / kernelWidth))) /
                          kernelWidth;
                }

                if (kbY != kbY) {//if kbY = NaN
                    kbY = 0;
                }

                /* kernel weighting value */
                //if (params.useLUT){
                //    w = kbX * kbY;
                //} else {

                w = kbX * kbY;
                
                /* grid data */
                idx = ny + (nx)*params.gridSize[1]/* + (nz)*gridOS*Nx*gridOS*Ny*/;
                //gridData[idx].x += (w*pt.real*atm);
                //gridData[idx].y += (w*pt.imag*atm);
                gridData[idx].real(gridData[idx].real()+w*pt.real);
                gridData[idx].imag(gridData[idx].imag()+w*pt.imag);
                /* estimate sample density */
                sampleDensity[idx] += w;
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
gridding_Gold_3D(unsigned int n, parameters<T1> params,T1 beta, ReconstructionSample<T1> *__restrict  sample,
                 const T1 *LUT, const uword sizeLUT,
                 complex<T1> *__restrict  gridData, T1 *__restrict  sampleDensity)
{
    unsigned int NxL, NxH;
    unsigned int NyL, NyH;
    unsigned int NzL, NzH;
    
    unsigned int nx;
    unsigned int ny;
    unsigned int nz;
    
    int idx;
    
    T1 w;
    
    T1 shiftedKx, shiftedKy, shiftedKz;
    T1 distX, kbX, distY, kbY, distZ, kbZ;
    
    T1 kernelWidth = params.kernelWidth;
    //T1 beta = 18.5547;
    T1 gridOS = params.gridOS;
    
    unsigned int Nx = params.imageSize[0];
    unsigned int Ny = params.imageSize[1];
    unsigned int Nz = params.imageSize[2];
    
    //Jiading GAI
    //float t0 = t[0];
    
    for (unsigned int i=0; i < n; i++)
    {
        ReconstructionSample<T1> pt = sample[i];
        
        //Jiading GAI
        //float atm = hanning_d(t[i], tau, l, t0);//a_l(t_m)
        
        shiftedKx = (gridOS)*(pt.kX+((T1)Nx)/2.0f);
        shiftedKy = (gridOS)*(pt.kY+((T1)Ny)/2.0f);
        shiftedKz = (gridOS)*(pt.kZ+((T1)Nz)/2.0f);
        
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

        NxL = (int)(fmax(0.0f,std::ceil(shiftedKx - kernelWidth*(gridOS)/2.0f)));
        NxH = (int)(fmin((gridOS*(T1)Nx-1.0f),std::floor(shiftedKx + kernelWidth*(gridOS)/2.0f)));

        NyL = (int)(fmax(0.0f,std::ceil(shiftedKy - kernelWidth*(gridOS)/2.0f)));
        NyH = (int)(fmin((gridOS*(T1)Ny-1.0f),std::floor(shiftedKy + kernelWidth*(gridOS)/2.0f)));

        NzL = (int)(fmax(0.0f,std::ceil(shiftedKz - kernelWidth*(gridOS)/2.0f)));
        NzH = (int)(fmin((gridOS*(T1)Nz-1.0f),std::floor(shiftedKz + kernelWidth*(gridOS)/2.0f)));

        for(nz=NzL; nz<=NzH; ++nz)
        {
            distZ = fabs(shiftedKz - ((T1)nz))/(gridOS);
            if (params.useLUT) {
                kbZ = kernel_value_LUT(distZ, LUT, sizeLUT, kernelWidth);
                //cout << "KBX = " << kbX << endl;

            } else {
                kbZ = bessi0(beta * std::sqrt(1.0 - (2.0 * distZ / kernelWidth) * (2.0 * distZ / kernelWidth))) /
                      kernelWidth;

            }

            if (kbZ != kbZ) {//if kbY = NaN
                kbZ = 0;
            }
            for(nx=NxL; nx<=NxH; ++nx)
            {
                distX = fabs(shiftedKx - ((T1) nx)) / (gridOS);
                if (params.useLUT) {
                    kbX = kernel_value_LUT(distX, LUT, sizeLUT, kernelWidth);
                    //cout << "KBX = " << kbX << endl;

                } else {
                    kbX = bessi0(beta * std::sqrt(1.0 - (2.0 * distX / kernelWidth) * (2.0 * distX / kernelWidth))) /
                          kernelWidth;

                }

                if (kbX != kbX) {//if kbX = NaN
                    kbX = 0;
                }

                for(ny=NyL; ny<=NyH; ++ny)
                {
                    distY = fabs(shiftedKy - ((T1) ny)) / (gridOS);
                    if (params.useLUT){
                        kbY = kernel_value_LUT(distY, LUT, sizeLUT, kernelWidth);
                    } else {
                        kbY = bessi0(
                                beta * std::sqrt(1.0 - (2.0 * distY / kernelWidth) * (2.0 * distY / kernelWidth))) /
                              kernelWidth;
                    }

                    if (kbY != kbY) {//if kbY = NaN
                        kbY = 0;
                    }

                    /* kernel weighting value */
                    //if (params.useLUT){
                    //    w = kbX * kbY;
                    //} else {

                    w = kbX * kbY * kbZ;

                    /* grid data */
                    idx = ny + (nx)*params.gridSize[1] + (nz)*params.gridSize[0]*params.gridSize[1];
                    //gridData[idx].x += (w*pt.real*atm);
                    //gridData[idx].y += (w*pt.imag*atm);
                    gridData[idx].real(gridData[idx].real()+w*pt.real);
                    gridData[idx].imag(gridData[idx].imag()+w*pt.imag);
                    
                    
                    /* estimate sample density */
                    sampleDensity[idx] += w;
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
gridding_Silver_2D(unsigned int n, parameters<T1> params,const T1  *kx, const T1 *ky, T1 beta, complex<T1> *__restrict sample,
                   const T1 *LUT, const uword sizeLUT,
                 complex<T1> * __restrict gridData, T1 *__restrict sampleDensity)
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

    T1 kernelWidth = params.kernelWidth;
    //T1 beta = 18.5547;
    T1 gridOS = params.gridOS;

    unsigned int Nx = params.imageSize[0];
    unsigned int Ny = params.imageSize[1];
    //unsigned int Nz = params.imageSize[2];

    //Jiading GAI
    //float t0 = t[0];
    #pramga acc parallel loop
    for (unsigned int i=0; i < n; i++)
    {
        complex<T1> pt = sample[i];

        //Jiading GAI
        //float atm = hanning_d(t[i], tau, l, t0);//a_l(t_m)

        shiftedKx = (gridOS)*(kx[i]+((T1)Nx)/2.0f);
        shiftedKy = (gridOS)*(ky[i]+((T1)Ny)/2.0f);

        //shiftedKz = ((float)gridOS)*(pt.kZ+((float)Nz)/2);

        //if(shiftedKx < 0.0f)
        //   shiftedKx = 0.0f;
        //if(shiftedKx > ((float)gridOS)*((float)Nx));
        //   shiftedKx = ((float)gridOS)*((float)Nx);
        //if(shiftedKy < 0.0f)
        //   shiftedKy = 0.0f;
        //if(shiftedKy > ((float)gridOS)*((float)Ny));
        //   shiftedKy = ((float)gridOS)*((float)Ny);


        NxL = (int)(fmax(0.0f,std::ceil(shiftedKx - kernelWidth*(gridOS)/2.0f)));
        NxH = (int)(fmin((gridOS*(T1)Nx-1.0f),std::floor(shiftedKx + kernelWidth*(gridOS)/2.0f)));

        NyL = (int)(fmax(0.0f,std::ceil(shiftedKy - kernelWidth*(gridOS)/2.0f)));
        NyH = (int)(fmin((gridOS*(T1)Ny-1.0f),std::floor(shiftedKy + kernelWidth*(gridOS)/2.0f)));

        //NzL = (int)(fmax(0.0f,ceil(shiftedKz - kernelWidth*((float)gridOS)/2)));
        //NzH = (int)(fmin((float)(gridOS*Nz-1),floor(shiftedKz + kernelWidth*((float)gridOS)/2)));

        for(nx=NxL; nx<=NxH; ++nx)
        {
            distX = fabs(shiftedKx - ((T1)nx))/(gridOS);
            if (params.useLUT){
                kbX = kernel_value_LUT(distX, LUT, sizeLUT, kernelWidth);
            } else {
                kbX = bessi0(beta*std::sqrt(1.0-(2.0*distX/kernelWidth)*(2.0*distX/kernelWidth)))/kernelWidth;

            }

            if (kbX != kbX) {//if kbX = NaN
                kbX = 0;
            }

            for(ny=NyL; ny<=NyH; ++ny)
            {
                distY = fabs(shiftedKy - ((T1)ny))/(gridOS);
                if (params.useLUT){
                    kbY = kernel_value_LUT(distY, LUT, sizeLUT, kernelWidth);
                } else {
                    kbY = bessi0(beta*std::sqrt(1.0-(2.0*distY/kernelWidth)*(2.0*distY/kernelWidth)))/kernelWidth;
                }

                if (kbY != kbY) {//if kbY = NaN
                    kbY = 0;
                }

                /* kernel weighting value */
                //if (params.useLUT){
                //    w = kbX * kbY;
                //} else {
                    w = kbX * kbY;
                //}
                /* grid data */
                idx = ny + (nx)*params.gridSize[1]/* + (nz)*gridOS*Nx*gridOS*Ny*/;
                //gridData[idx].x += (w*pt.real*atm);
                //gridData[idx].y += (w*pt.imag*atm);
                sample[i].real(sample[i].real()+w*gridData[idx].real());
                sample[i].imag(sample[i].imag()+w*gridData[idx].imag());
                /* estimate sample density */
                sampleDensity[i] += w;
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
gridding_Silver_3D(unsigned int n, parameters<T1> params,const T1  *kx, const T1 *ky, const T1 *kz, T1 beta, complex<T1> *__restrict  sample,
                   const T1 *LUT, const uword sizeLUT,
                 complex<T1> *__restrict  gridData, T1 *__restrict  sampleDensity)
{
    unsigned int NxL, NxH;
    unsigned int NyL, NyH;
    unsigned int NzL, NzH;

    unsigned int nx;
    unsigned int ny;
    unsigned int nz;

    int idx;

    T1 w;

    T1 shiftedKx, shiftedKy, shiftedKz;
    T1 distX, kbX, distY, kbY, distZ, kbZ;

    T1 kernelWidth = params.kernelWidth;
    //T1 beta = 18.5547;
    T1 gridOS = params.gridOS;

    unsigned int Nx = params.imageSize[0];
    unsigned int Ny = params.imageSize[1];
    unsigned int Nz = params.imageSize[2];

    //Jiading GAI
    //float t0 = t[0];

    for (unsigned int i=0; i < n; i++)
    {
        complex<T1> pt = sample[i];

        //Jiading GAI
        //float atm = hanning_d(t[i], tau, l, t0);//a_l(t_m)

        shiftedKx = (gridOS)*(kx[i]+((T1)Nx)/2.0f);
        shiftedKy = (gridOS)*(ky[i]+((T1)Ny)/2.0f);
        shiftedKz = (gridOS)*(kz[i]+((T1)Nz)/2.0f);

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


        NxL = (int)(fmax(0.0f,std::ceil(shiftedKx - kernelWidth*(gridOS)/2.0f)));
        NxH = (int)(fmin((gridOS*(T1)Nx-1.0f),std::floor(shiftedKx + kernelWidth*((float)gridOS)/2.0f)));

        NyL = (int)(fmax(0.0f,std::ceil(shiftedKy - kernelWidth*(gridOS)/2.0f)));
        NyH = (int)(fmin((gridOS*(T1)Ny-1.0f),std::floor(shiftedKy + kernelWidth*((float)gridOS)/2.0f)));

        NzL = (int)(fmax(0.0f,std::ceil(shiftedKz - kernelWidth*(gridOS)/2.0f)));
        NzH = (int)(fmin((gridOS*(T1)Nz-1.0f),std::floor(shiftedKz + kernelWidth*((float)gridOS)/2.0f)));

        for(nz=NzL; nz<=NzH; ++nz)
        {
            distZ = fabs(shiftedKz - ((T1)nz))/(gridOS);
            if (params.useLUT) {
                kbZ = kernel_value_LUT(distZ, LUT, sizeLUT, kernelWidth);
                //cout << "KBX = " << kbX << endl;

            } else {
                kbZ = bessi0(beta * std::sqrt(1.0 - (2.0 * distZ / kernelWidth) * (2.0 * distZ / kernelWidth))) /
                      kernelWidth;

            }

            if (kbZ != kbZ) {//if kbY = NaN
                kbZ = 0;
            }
            for(nx=NxL; nx<=NxH; ++nx)
            {
                distX = fabs(shiftedKx - ((T1) nx)) / (gridOS);
                if (params.useLUT) {
                    kbX = kernel_value_LUT(distX, LUT, sizeLUT, kernelWidth);
                    //cout << "KBX = " << kbX << endl;

                } else {
                    kbX = bessi0(beta * std::sqrt(1.0 - (2.0 * distX / kernelWidth) * (2.0 * distX / kernelWidth))) /
                          kernelWidth;

                }

                if (kbX != kbX) {//if kbX = NaN
                    kbX = 0;
                }

                for (ny = NyL; ny <= NyH; ++ny) {
                    distY = fabs(shiftedKy - ((T1) ny)) / (gridOS);
                    if (params.useLUT) {
                        kbY = kernel_value_LUT(distY, LUT, sizeLUT, kernelWidth);
                    } else {
                        kbY = bessi0(
                                beta * std::sqrt(1.0 - (2.0 * distY / kernelWidth) * (2.0 * distY / kernelWidth))) /
                              kernelWidth;
                    }

                    if (kbY != kbY) {//if kbY = NaN
                        kbY = 0;
                    }
                    /* kernel weighting value */
                    //if (params.useLUT){
                    //    w = kbX * kbY;
                    //} else {

                    w = kbX * kbY * kbZ;

                    /* grid data */
                    idx = ny + (nx)*params.gridSize[1] + (nz)*params.gridSize[0]*params.gridSize[1];
                    //gridData[idx].x += (w*pt.real*atm);
                    //gridData[idx].y += (w*pt.imag*atm);
                    sample[i].real(sample[i].real()+w*gridData[idx].real());
                    sample[i].imag(sample[i].imag()+w*gridData[idx].imag());


                    /* estimate sample density */
                    sampleDensity[i] += w;
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
                   int numK_per_coil, const T1  *__restrict kx, const T1 *__restrict ky, const T1 *__restrict kz,
                   const T1 *__restrict dR, const T1 *__restrict dI, int Nx, int Ny, int Nz,
                   T1 gridOS, T1 *__restrict outR_d, T1 *__restrict outI_d, const T1 kernelWidth, const T1 beta,
                   const T1 *LUT, const uword sizeLUT)
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

    parameters<T1> params;
    params.sync=0;
    params.binsize=128;

    params.useLUT = 1;
    params.kernelWidth = kernelWidth;
    params.gridOS = gridOS;
    params.imageSize[0] = Nx;//gridSize is gridOS times larger than imageSize.
    params.imageSize[1] = Ny;
    params.imageSize[2] = Nz;
    params.gridSize[0]  = std::ceil(gridOS*(T1)Nx);
    params.gridSize[1]  = std::ceil(gridOS*(T1)Ny);
    if(params.gridSize[0]%2)//3D case, gridOS is adjusted on the z dimension:
        params.gridSize[0] += 1;//That why we need to make sure here that the xy
    if(params.gridSize[1]%2)//dimensions have even sizes.
        params.gridSize[1] += 1;
    params.gridSize[2]  = (Nz==1)?Nz:(std::ceil(gridOS*(T1)Nz));// 2D or 3D
    params.numSamples = numK_per_coil;

    T1 *sampleDensity;


    ReconstructionSample<T1>* samples; //Input Data
    //allocate samples
    samples = (ReconstructionSample<T1>*) malloc(params.numSamples*sizeof(ReconstructionSample<T1>));
    
    if (samples == NULL){
        printf("ERROR: Unable to allocate memory for input data\n");
        exit(1);
    }
    
    unsigned int n =  params.numSamples;
    //
    for(int i=0; i<params.numSamples; i++){
        if( std::abs(kx[i])>(Nx/2.0) ||
                std::abs(ky[i])>(Ny/2.0) ||
                std::abs(kz[i])>(Nz/2.0)
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
            
            samples[i].sdc = 1.0f;
            //samples[i].t = t[i];
        }
    }    ///*
    // grid_size in xy-axis has to be divisible-by-two:
    //       (required by the cropImageRegion)
    // grid_size in z-axis has to be divisible-by-four:
    //       (required by the function gridding_GPU_3D(.))
    if(1==Nz) {
        //round grid size (xy-axis) to the next divisible-by-two.
        gridOS = 2.0f * std::ceil((gridOS * (T1)Nx) / 2.0f) / (T1) Nx;
    }
    else {
        //round grid size (z-axis) to the next divisible-by-four.
        gridOS = 4.0f * std::ceil((gridOS * (T1)Nz) / 4.0f) / (T1) Nz;
    }
    // */
    int gridNumElems  = params.gridSize[0] *
    params.gridSize[1] *
    params.gridSize[2] ;
    
    int imageNumElems = params.imageSize[0] *
    params.imageSize[1] *
    params.imageSize[2] ;


    //allocate gridData
    complex<T1> *gridData = new complex<T1>[gridNumElems];
    sampleDensity = new T1[gridNumElems];

    // Have to set 'gridData' and 'sampleDensity' to zero.
    // Because they will be involved in accumulative operations
    // inside gridding functions.
    for(int i=0;i<gridNumElems;i++)
    {
        gridData[i].real(0.0);
        gridData[i].imag(0.0);
        sampleDensity[i] = 0.0;
    }

    // Gridding with CPU - gold
    if(Nz==1)
    {
        gridding_Gold_2D<T1>(n, params, beta, samples, LUT, sizeLUT,
                         gridData, sampleDensity);
    }
    else
    {
        gridding_Gold_3D<T1>(n, params, beta, samples, LUT, sizeLUT,
                         gridData, sampleDensity);

    }
    //cx_vec temp(gridData,gridNumElems);
    //savemat("/shared/mrfil-data/data/PowerGridTest/64_64_16_4coils/gridData.mat","img",temp);

    complex<T1> *gridData_d = new complex<T1>[gridNumElems];
    //memcpy(gridData_d,gridData,gridNumElems*sizeof(complex<T1>));

    // ifftshift(gridData):
    if(Nz==1)
    {
        ifftshift2(gridData_d,gridData,params.gridSize[0], params.gridSize[1]);
    }
    else
    {
        ifftshift3(gridData_d,gridData,params.gridSize[0],params.gridSize[1],params.gridSize[2]);
    }
    //cx_vec temp2(gridData_d,gridNumElems);
    //savemat("/shared/mrfil-data/data/PowerGridTest/64_64_16_4coils/ifftshiftData.mat","img",temp2);

    // ifftn(gridData):
    fftw_plan plan;
    if(Nz==1)
    {
        plan = fftw_plan_dft_2d(params.gridSize[0],
                                     params.gridSize[1], (fftw_complex *)gridData_d, (fftw_complex *)gridData_d, FFTW_BACKWARD, FFTW_ESTIMATE);
    }
    else
    {
        plan = fftw_plan_dft_3d(params.gridSize[2],
                                params.gridSize[1], params.gridSize[0], (fftw_complex *)gridData_d, (fftw_complex *)gridData_d, FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    /* Inverse transform 'gridData_d' in place. */
    fftw_execute(plan);
    fftw_destroy_plan(plan);
    // fftshift(gridData):
    if(Nz==1)
    {
        fftshift2(gridData, gridData_d,params.gridSize[0],
                            params.gridSize[1]);
    }
    else
    {
        fftshift3(gridData, gridData_d, params.gridSize[0],
                            params.gridSize[1], params.gridSize[2]);
    }

    complex<T1> *gridData_crop_d = new complex<T1>[imageNumElems];
    //memcpy(gridData_crop_d,gridData_d,gridNumElems*sizeof(complex<T1>));
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
    // deapodization
    if(Nz==1)
    {
        deapodization2d(gridData_crop_d, gridData_crop_d,
                        Nx, Ny, kernelWidth, beta, params.gridOS);
    }
    else
    {
        deapodization3d(gridData_crop_d, gridData_crop_d,
                        Nx, Ny, Nz, kernelWidth, beta, params.gridOS);
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

    //deallocate samples
    free(samples);
    delete(gridData);
    delete(gridData_d);
    delete(sampleDensity);
    delete(gridData_crop_d);

}
//Calculates the gridded forward fourier transform
template<typename T1>
void
computeFd_CPU_Grid(
        int numK_per_coil, const T1  *__restrict kx, const T1 *__restrict ky, const T1 *__restrict kz,
        const T1 *__restrict dR, const T1 *__restrict dI, int Nx, int Ny, int Nz,
        T1 gridOS, T1 *__restrict outR_d, T1 *__restrict outI_d, const T1 kernelWidth, const T1 beta, const T1 *LUT,
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
    T1 kernelWidth = 4.0;
    T1 beta = MRI_PI * std::sqrt( (gridOS - 0.5) * (gridOS - 0.5) *
                                  (kernelWidth * kernelWidth*4.0) /
                                  (gridOS * gridOS) - 0.8
    );
    */
    parameters<T1> params;
    params.sync=0;
    params.binsize=128;

    params.useLUT = 1;
    params.kernelWidth = kernelWidth;
    params.gridOS = gridOS;
    params.imageSize[0] = Nx;//gridSize is gridOS times larger than imageSize.
    params.imageSize[1] = Ny;
    params.imageSize[2] = Nz;
    params.gridSize[0]  = std::ceil(gridOS*(T1)Nx);
    params.gridSize[1]  = std::ceil(gridOS*(T1)Ny);
    if(params.gridSize[0]%2)//3D case, gridOS is adjusted on the z dimension:
        params.gridSize[0] += 1;//That why we need to make sure here that the xy
    if(params.gridSize[1]%2)//dimensions have even sizes.
        params.gridSize[1] += 1;
    params.gridSize[2]  = (Nz==1)?Nz:(std::ceil(gridOS*(T1)Nz));// 2D or 3D
    params.numSamples = numK_per_coil;

    T1 *sampleDensity;


    complex<T1>* samples; //Input Data
    //allocate samples
    samples = (complex<T1>*) malloc(params.numSamples*sizeof(complex<T1>));

    if (samples == NULL){
        printf("ERROR: Unable to allocate memory for input data\n");
        exit(1);
    }

    unsigned int n =  params.numSamples;
    //
    for(int i=0; i<params.numSamples; i++){
        if( std::abs(kx[i])>(Nx/2.0) ||
            std::abs(ky[i])>(Ny/2.0) ||
            std::abs(kz[i])>(Nz/2.0)
                ) {

            printf("\nError:k-space trajectory out of range [-N/2,N/2]:\n      gridding requires that k-space should be contained within the window -N/2 to N/2.\n");
            cout << "kx = " << kx[i] << " ky = " << ky[i] << " kz = " << kz[i] << " i = " << i << endl;
            exit(1);
        }
        else {



            samples[i].real(0.0);
            samples[i].imag(0.0);

        }
    }


     ///*
    // grid_size in xy-axis has to be divisible-by-two:
    //       (required by the cropImageRegion)
    // grid_size in z-axis has to be divisible-by-four:
    //       (required by the function gridding_GPU_3D(.))

    if(1==Nz) {
        //round grid size (xy-axis) to the next divisible-by-two.
        gridOS = 2.0f * std::ceil((gridOS * (T1)Nx) / 2.0f) / (T1) Nx;
    }
    else {
        //round grid size (z-axis) to the next divisible-by-four.
        gridOS = 4.0f * std::ceil((gridOS * (T1)Nz) / 4.0f) / (T1) Nz;
    }
    //

    int gridNumElems  = params.gridSize[0] *
                        params.gridSize[1] *
                        params.gridSize[2] ;

    int imageNumElems = params.imageSize[0] *
                        params.imageSize[1] *
                        params.imageSize[2] ;


    //allocate gridData
    complex<T1> *gridData = new complex<T1>[imageNumElems];
    sampleDensity = new T1[gridNumElems];

    // Have to set 'gridData' and 'sampleDensity' to zero.
    // Because they will be involved in accumulative operations
    // inside gridding functions.
    for(int i=0;i<imageNumElems;i++)
    {
        gridData[i].real(dR[i]);
        gridData[i].imag(dI[i]);
        sampleDensity[i] = 0.0;
    }

    complex<T1> *gridData_d = new complex<T1>[imageNumElems];
    //memcpy(gridData_d,gridData,imageNumElems*sizeof(complex<T1>));

    // deapodization

    if(Nz==1)
    {
        deapodization2d(gridData_d, gridData,
                        Nx, Ny, kernelWidth, beta, params.gridOS);
    }
    else
    {
        deapodization3d(gridData_d, gridData,
                        Nx, Ny, Nz, kernelWidth, beta, params.gridOS);
    }

    complex<T1> *gridData_os = new complex<T1>[gridNumElems];

    //zero pad

    if(Nz==1)
    {
        zero_pad2d(gridData_os, gridData_d,
                        Nx, Ny, params.gridOS);
    }
    else
    {

        zero_pad3d(gridData_os, gridData_d,
                   Nx, Ny, Nz, params.gridOS);

    }

    complex<T1> *gridData_os_d = new complex<T1>[gridNumElems];
    //memcpy(gridData_os_d, gridData_os, gridNumElems*sizeof(complex<T1>));
    // fftshift(gridData):

    if(Nz==1)
    {
        fftshift2(gridData_os_d, gridData_os,params.gridSize[0],
                  params.gridSize[1]);
    }
    else
    {
        fftshift3(gridData_os_d, gridData_os, params.gridSize[0],
                  params.gridSize[1], params.gridSize[2]);
    }

    // ifftn(gridData):

    fftw_plan plan;
    if(Nz==1)
    {
        plan = fftw_plan_dft_2d(params.gridSize[0],
                                params.gridSize[1], (fftw_complex *)gridData_os_d, (fftw_complex *)gridData_os_d, FFTW_FORWARD, FFTW_ESTIMATE);
    }
    else
    {
        plan = fftw_plan_dft_3d(params.gridSize[2],
                                params.gridSize[1], params.gridSize[0], (fftw_complex *)gridData_os_d, (fftw_complex *)gridData_os_d, FFTW_FORWARD, FFTW_ESTIMATE);
    }

    /* Inverse transform 'gridData_d' in place. */
    fftw_execute(plan);
    fftw_destroy_plan(plan);

    // ifftshift(gridData):

    if(Nz==1)
    {
        ifftshift2(gridData_os, gridData_os_d, params.gridSize[0], params.gridSize[1]);
    }
    else
    {
        ifftshift3(gridData_os, gridData_os_d, params.gridSize[0], params.gridSize[1], params.gridSize[2]);
    }


    // Gridding with CPU - silver
    if(Nz==1)
    {
        gridding_Silver_2D<T1>(n, params, kx, ky, beta, samples, LUT, sizeLUT,
                               gridData_os, sampleDensity);
    }
    else
    {
        gridding_Silver_3D<T1>(n, params, kx, ky, kz, beta, samples, LUT, sizeLUT,
                               gridData_os, sampleDensity);
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

    for(int ii = 0; ii < n; ii++) {
        outR_d[ii] = samples[ii].real();
        outI_d[ii] = samples[ii].imag();
    }

    //deallocate samples
    free(samples);
    //free(LUT);
    delete(gridData);
    delete(gridData_d);
    delete(gridData_os);
    delete(gridData_os_d);
    delete(sampleDensity);
    //delete(gridData_crop_d);
}

#endif
