//
//  test_3D.h
//  PowerGrid
//
//  Created by Joe holtrop on 4/7/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef __PowerGrid__test_3D__
#define __PowerGrid__test_3D__

#include "PowerGrid.h"

using namespace arma;

template<typename T1,typename T2>
void test_3D(string dataPath)
{
    string testPath = dataPath;

    uword Nx =64;
    uword Ny =64;
    uword Nz =16;

    //Setup image space coordinates/trajectory
    Mat<T2> ix(Nx,Ny,Nz);
    Mat<T2> iy(Nx,Ny,Nz);
    Mat<T2> iz(Nx,Ny,Nz);
    Mat<T2> FM(Nx,Ny,Nz);
    
    ix.zeros();
    iy.zeros();
    iz.zeros();
    FM.zeros();
    
    //generate the image space cordinates of the voxels we want to reconstruct
    // after vectorizing ix and iy the image coordinates must match the Field and SENSe map image coordinates
    for(uword ii = 0; ii < 64; ii++) { //y
        for (uword jj = 0; jj < 64; jj++) { //x

            ix(ii,jj) = ((T2)jj-32.0)/64.0;
            iy(ii,jj) = ((T2)ii-32.0)/64.0;

        }
    }

    Col<T2> kx;
    loadmat(testPath+"kx.mat","kx",&kx);
    Col<T2> ky;
    loadmat(testPath+"ky.mat","ky",&ky);
    Col<T2> kz;
    loadmat(testPath+"kz.mat","kz",&kz);

    uword nro;
    nro = kx.n_elem;

    Col<T2> tvec;
    loadmat(testPath+"t.mat","t",&tvec);
    
    uword L = 10;
    loadmat(testPath+"FM.mat","FM",&FM);

    uword nc = 4;
    loadmat(testPath+"SMap.mat","SMap",&SMap);

    // Fourier transfrom operator
    Ggrid<T1,T2> G(nro,Nx,Ny,Nz,kx,ky,kz,vectorise(ix),vectorise(iy),vectorise(iz));

    // Field correction operation
    FieldCorrection<T1, T2, Ggrid<T1,T2>> A(G,vectorise(FM),vectorise(tvec),nro,Nx*Ny*Nz,L);

   // Sense operation
    SENSE<cx_double, FieldCorrection<T1, T2, Ggrid<T1,T2>>> S(A,SMap,nro*nc,Nx*Ny*Nz,nc);

    // Variables needed for the recon: Penalty object, num of iterations
    umat ReconMask;
    ReconMask.ones(Nx,Ny,Nz);
    
    QuadPenalty<T1>R(Nx,Ny,Nz,0,ReconMask);
    
    uword niter = 10;
    Col<T1> xinit = zeros<Col<cx_double>>(Nx,Ny,Nz); // initial estimate of x
    Mat<T1> W;
    W = eye<Mat<T1>>(A.n1,A.n1); // Should be the size of k-space data: Is it right?

    Col<T1> data;
    loadmat(testPath+"data.mat","data",&data);

    Col<T1> x_t;
    x_t = pwls_pcg1<T1, SENSE<cx_double, FieldCorrection<T1, T2, Ggrid<T1,T2>>>,QuadPenalty<T1>>(xinit, S, W, data, R, niter);

    savemat(testPath+"test_3D.mat","img",x_t);
    
}
#endif
