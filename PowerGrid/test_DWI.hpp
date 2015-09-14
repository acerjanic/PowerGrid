//
//  test_DWI.h
//  PowerGrid
//
//  Created by Joe holtrop on 8/13/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef __PowerGrid__DWI
#define __PowerGrid__DWI

#include "PowerGrid.h"

using namespace arma;

template<typename T1>
int test_DWI(string dataPath, uword Nx, uword Ny, uword Nz, uword L, uword niter,uword nc) {
    string testPath = dataPath;
    typedef complex<T1> CxT1;
    //uword Nx =64;
    //uword Ny =64;
    //uword Nz =16;

    //Setup image space coordinates/trajectory
    Cube<T1> ix(Nx,Ny,Nz);
    Cube<T1> iy(Nx,Ny,Nz);
    Cube<T1> iz(Nx,Ny,Nz);
    Col<T1> FM;
    Col<CxT1> SMap;
	Col<T1> PMap;
    
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
    
    //uword L = 1;
    loadmat(testPath+"FM.mat","FM",&FM);
    //FM.zeros();

    //uword nc = 4;
    loadmat(testPath+"SMap.mat","SMap",&SMap);

	//nonlinear motion correction for DWI
	loadmat(testPath+"PMap.mat","PMap",&PMap);
	//PMap.zeros();
	//The navigator phase is being referenced to zero
	DWICGMC<T1> S_DWI( kx, ky, kz, Nx, Ny, Nz, nc, tvec, SMap, vectorise(FM), 0-PMap);

    cout << "loading data" << endl;
    Col<CxT1> data;
    loadmat(testPath+"data.mat","data",&data);

    // Variables needed for the recon: Penalty object, num of iterations

    cout << "Iniitalizing QuadPenalty" << endl;
    //TVPenalty<T1>R(Nx,Ny,Nz,100000000.0,.000001);
	//TVPenalty<T1>R(Nx,Ny,Nz,0.0,.01);
	QuadPenalty<T1>R(Nx,Ny,Nz,0.0);
    cout << "QuadPenalty setup successfull" << endl;

    //uword niter = 10;
    Col<CxT1> xinit(Nx*Ny*Nz); // initial estimate of x
	xinit.zeros();
    Col<T1> W;
    //W = eye<sp_mat<T1>>(A.n1,A.n1); // Should be the size of k-space data: Is it right?
    //W=1.0;
	W = ones(nro*nc);

    cout << "Runing pwls with ggrid" << endl;
    Col<CxT1> test_DWI_img;
    test_DWI_img = solve_pwls_pcg<T1,  DWICGMC<T1> ,QuadPenalty<T1>>(xinit, S_DWI, W, data, R, niter);
	savemat(testPath+"test_DWICGMC.mat","img",test_DWI_img);

    return 0;
    
}

#endif
