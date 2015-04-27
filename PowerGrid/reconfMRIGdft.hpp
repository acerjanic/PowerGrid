//
//  test_3D.h
//  PowerGrid
//
//  Created by Joe holtrop on 4/7/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef __PowerGrid__reconfMRIGdft__
#define __PowerGrid__reconfMRIGdft__
#include <sstream>
#include "PowerGrid.h"
using namespace arma;

template<typename T1,typename T2>
int reconfMRIGdft(string dataPath, uword Nx, uword Ny, uword Nz, uword L, uword niter,uword nc, uword startIndex, uword endIndex ) {
    string testPath = dataPath;


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
    //loadmat(testPath+"kz.mat","kz",&kz);
   kz.zeros(kx.n_elem);
    uword nro;
    nro = kx.n_elem;

    Col<T2> tvec;
    loadmat(testPath+"t.mat","t",&tvec);
    

    // Fourier transfrom operator
  //  cout << "Initializing Ggrid" << endl;

    // Field correction operation
   // cout << "Initializing FieldCorrection" << endl;

    //cout << "Initializing Gdft" << endl;

    //uword nc = 4;

   // cout << "Iniitalizing SENSE gdft" << endl;
   // SENSE<cx_double, Gdft<T1,T2>> Sd(Gd,SMap,nro,Nx*Ny*Nz,nc);

    // Sense operation
    //cout << "Iniitalizing SENSE Ggrid" << endl;
    //SENSE<cx_double, FieldCorrection<T1, T2, Ggrid<T1,T2>>> Sg(A,SMap,nro,Nx*Ny*Nz,nc);

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
    

    cout << " Start recon data with Gdft"<< endl;
    int slices = 20;

    for(int jj = 1; jj< slices +1; jj++)   
    {
        for (int ii = startIndex; ii< endIndex +1; ii++)
        {   
   string name1;
   std::stringstream ss1;
   ss1 <<  ii;
   name1 = ss1.str();//passing the time point number

   string name2;
   std::stringstream ss2;
   ss2 <<  jj;
   name2 = ss2.str(); // passing the slice number
   
    Col<T1> data;
    loadmat(testPath + "data_" + name1 +"_" + name2 + "_data.mat","data",&data);
     
    Col<T2> FM;
    loadmat(testPath + "FM_" + name2 + ".mat","FM",&FM);
    Gdft<T1,T2> Gd(nro,Nx*Ny*Nz,kx,ky,kz,vectorise(ix),vectorise(iy),vectorise(iz),vectorise(FM),vectorise(tvec));

    cout << " data and FM loaded" << endl;
    Col<T1> SMap; 
    loadmat(testPath + "SMap" + name2 + ".mat","SMap",&SMap);
    SENSE<cx_double,  Gdft<T1,T2>> Sd(Gd,SMap,nro,Nx*Ny*Nz,nc);
  
    cout << "sense created" << endl;
    Col<T1> img;
    img = pwls_pcg1<T1, SENSE< T1 , Gdft<T1,T2>>,QuadPenalty<T1>>(xinit, Sd, W, data, R, niter);
    savemat(testPath + "img_gdft_" + name1 + "_" + name2 + ".mat","img",img);
       }
   }
    return 0;
    
}

#endif


