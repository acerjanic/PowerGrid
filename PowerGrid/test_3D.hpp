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

template<typename T1>
int test_3D(string dataPath)
{
    typedef complex<T1> CxT1;
    string testPath = dataPath;

    uword Nx =64;
    uword Ny =64;
    uword Nz =16;

    //Setup image space coordinates/trajectory
    Cube<T1> ix(Nx,Ny,Nz);
    Cube<T1> iy(Nx,Ny,Nz);
    Cube<T1> iz(Nx,Ny,Nz);
    Col<T1> FM;
    Col<CxT1> SMap;
    
    ix.zeros();
    iy.zeros();
    iz.zeros();
    
    //generate the image space cordinates of the voxels we want to reconstruct
    // after vectorizing ix and iy the image coordinates must match the Field and SENSe map image coordinates
    for(uword ii = 0; ii < Ny; ii++) { //y
        for (uword jj = 0; jj < Nx; jj++) { //x
            for (uword kk = 0; kk < Nz; kk++) { //z

                ix(ii, jj, kk) = ((T1) jj - (T1) Nx / 2.0) / ((T1) Nx);
                iy(ii, jj, kk) = ((T1) ii - (T1) Ny / 2.0) / ((T1) Ny);
                iz(ii, jj, kk) = ((T1) kk - (T1) Nz / 2.0) / ((T1) Nz);
            }
        }
    }

    Col<T1> kx;
    loadmat(testPath+"kx.mat","kx",&kx);
    Col<T1> ky;
    loadmat(testPath+"ky.mat","ky",&ky);
    Col<T1> kz;
    loadmat(testPath+"kz.mat","kz",&kz);

    uword nro;
    nro = kx.n_elem;

    Col<T1> tvec;
    loadmat(testPath+"t.mat","t",&tvec);
    
    uword L = 1;
    loadmat(testPath+"FM.mat","FM",&FM);
    FM.zeros();

    uword nc = 4;
    loadmat(testPath+"SMap.mat","SMap",&SMap);

    // Fourier transfrom operator
    cout << "Initializing Ggrid" << endl;
    Ggrid<T1> G(nro,2.0,Nx,Ny,Nz,kx,ky,kz,vectorise(ix),vectorise(iy),vectorise(iz));

    // Field correction operation
    cout << "Initializing FieldCorrection" << endl;
    FieldCorrection<T1, Ggrid<T1>> A(G,vectorise(FM),vectorise(tvec),nro,Nx*Ny*Nz,L);

   // Sense operation
    cout << "Iniitalizing SENSE" << endl;
    SENSE<T1, FieldCorrection<T1, Ggrid<T1>>> S(A,SMap,nro,Nx*Ny*Nz,nc);
    //SENSE<cx_double, Ggrid<T1,T2>> S(G,SMap,nro*nc,Nx*Ny*Nz,nc);

    // Variables needed for the recon: Penalty object, num of iterations
    ucube ReconMask(Nx,Ny,Nz);
    ReconMask.ones();

    cout << "Iniitalizing QuadPenalty" << endl;
    QuadPenalty<T1>R(Nx,Ny,Nz,0,ReconMask);
    cout << "QuadPenalty setup successfull" << endl;

    uword niter = 5;
    Col<CxT1> xinit(Nx*Ny*Nz); // initial estimate of x
    xinit.zeros();
    Col<T1> W;
    //W = eye<sp_mat<T1>>(A.n1,A.n1); // Should be the size of k-space data: Is it right?
    W=ones(nro*nc);


    cout << "loading data" << endl;
    Col<CxT1> data;
    loadmat(testPath+"data.mat","data",&data);

    Col<CxT1> x_t;
    cout << "heading into solve_pwls_pcg" << endl;
    //x_t = solve_pwls_pcg<T1, SENSE<cx_double, FieldCorrection<T1, T2, Ggrid<T1,T2>>>,QuadPenalty<T1>>(xinit, S, W, data, R, niter);
    x_t = S/data;
    savemat(testPath+"test_3D.mat","img",x_t);

    return 0;
    
}

#endif
