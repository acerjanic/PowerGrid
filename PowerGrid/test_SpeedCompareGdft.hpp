//
//  test_SpeedCompareGdft.h
//  PowerGrid
//
//  Created by Joe holtrop on 4/7/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef __PowerGrid__SpeedCompare__Gdft
#define __PowerGrid__SpeedCompare__Gdft

#include "PowerGrid.h"

using namespace arma;

template<typename T1,typename T2>
int test_SpeedCompareGdft(string dataPath, uword Nx, uword Ny, uword Nz, uword L, uword niter,uword nc) {
    string testPath = dataPath;

    //uword Nx =64;
    //uword Ny =64;
    //uword Nz =16;

    //Setup image space coordinates/trajectory
    Cube<T2> ix(Nx,Ny,Nz);
    Cube<T2> iy(Nx,Ny,Nz);
    Cube<T2> iz(Nx,Ny,Nz);
    Col<T2> FM;
    Col<T1> SMap;
    
    ix.zeros();
    iy.zeros();
    iz.zeros();
    
    //generate the image space cordinates of the voxels we want to reconstruct
    // after vectorizing ix and iy the image coordinates must match the Field and SENSe map image coordinates
    for(uword ii = 0; ii < Ny; ii++) { //y
        for (uword jj = 0; jj < Nx; jj++) { //x
            for (uword kk = 0; kk < Nz; kk++) { //z

                ix(ii, jj, kk) = ((T2) jj - (T2) Nx / 2.0) / ((T2) Nx);
                iy(ii, jj, kk) = ((T2) ii - (T2) Ny / 2.0) / ((T2) Ny);
                iz(ii, jj, kk) = ((T2) kk - (T2) Nz / 2.0) / ((T2) Nz);
            }
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
    
    //uword L = 1;
    loadmat(testPath+"FM.mat","FM",&FM);
    FM.zeros();



    // Fourier transfrom operator
    // cout << "Initializing Ggrid" << endl;
    // Ggrid<T1,T2> Gg(nro,2.0,Nx,Ny,Nz,kx,ky,kz,vectorise(ix),vectorise(iy),vectorise(iz));

    // Field correction operation
    //cout << "Initializing FieldCorrection" << endl;
    //FieldCorrection<T1, T2, Ggrid<T1,T2>> A(Gg,vectorise(FM),vectorise(tvec),nro,Nx*Ny*Nz,L);

    cout << "Initializing Gdft" << endl;
    Gdft<T1,T2> Gd(nro,Nx*Ny*Nz,kx,ky,kz,vectorise(ix),vectorise(iy),vectorise(iz),vectorise(FM),vectorise(tvec));

    //uword nc = 4;
    loadmat(testPath+"SMap.mat","SMap",&SMap);

    cout << "Iniitalizing SENSE gdft" << endl;
    SENSE<cx_double, Gdft<T1,T2>> Sd(Gd,SMap,nro,Nx*Ny*Nz,nc);

    // Sense operation
    //cout << "Iniitalizing SENSE Ggrid" << endl;
    //SENSE<cx_double, FieldCorrection<T1, T2, Ggrid<T1,T2>>> Sg(A,SMap,nro,Nx*Ny*Nz,nc);

    cout << "loading data" << endl;
    Col<T1> data;
    loadmat(testPath+"data.mat","data",&data);

    // Variables needed for the recon: Penalty object, num of iterations
    ucube ReconMask(Nx,Ny,Nz);
    ReconMask.ones();

    cout << "Iniitalizing QuadPenalty" << endl;
    QuadPenalty<T1>R(Nx,Ny,Nz,0);
    cout << "QuadPenalty setup successfull" << endl;

    //uword niter = 10;
    Col<T1> xinit(Nx*Ny*Nz); // initial estimate of x
    xinit.zeros();
    T2 W;
    //W = eye<sp_mat<T1>>(A.n1,A.n1); // Should be the size of k-space data: Is it right?
    W=1.0;

    //Col<T1> x_t;
    //cout << "heading into PWLS_pcg1" << endl;
    //x_t = pwls_pcg1<T1, SENSE<cx_double, FieldCorrection<T1, T2, Ggrid<T1,T2>>>,QuadPenalty<T1>>(xinit, S, W, data, R, niter);
    //x_t = S/data;
    //savemat(testPath+"test_3D.mat","img",x_t);
/*
    cout << "Runing adjoint with ggrid" << endl;
    Col<T1> x_Sg_adjoint;
    x_Sg_adjoint = Sg/data;
    savemat(testPath+"test_adjoint_Sggrid.mat","img",x_Sg_adjoint);

    cout << "Runing forward with ggrid" << endl;
    Col<T1> x_Sg_forward;
    x_Sg_forward = Sg*x_Sg_adjoint;
    savemat(testPath+"test_forward_Sggrid.mat","img",x_Sg_forward);

    cout << "Runing adjoint with gdft" << endl;
    Col<T1> x_Sd_adjoint;
    x_Sd_adjoint = Sd/data;
    savemat(testPath+"test_adjoint_Sdft.mat","img",x_Sd_adjoint);

    cout << "Runing forward with gdft" << endl;
    Col<T1> x_Sd_forward;
    x_Sd_forward = Sd*x_Sd_adjoint;
    savemat(testPath+"test_forward_Sdft.mat","img",x_Sd_forward);
*/
    cout << "Runing pwls with ggrid" << endl;
    Col<T1> test_pwls;
    test_pwls = pwls_pcg1<T1,SENSE<cx_double, Gdft<T1,T2>>,QuadPenalty<T1>>(xinit, Sd, W, data, R, niter);
    savemat(testPath+"test_pwls.mat","img",test_pwls);
/*
    cout << "Runing pwls with ggrid" << endl;
    Col<T1> test_pwls;
    test_pwls = pwls_pcg1<T1,  Ggrid<T1,T2>,QuadPenalty<T1>>(xinit, Gg, W, data, R, niter);
    savemat(testPath+"test_pwls.mat","img",test_pwls);
    */

    return 0;
    
}

#endif
