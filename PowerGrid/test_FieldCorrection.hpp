//
//  test_FieldCorrection.h
//  PowerGrid
//
//  Created by Joe holtrop on 4/7/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef __PowerGrid__test_FieldCorrection__
#define __PowerGrid__test_FieldCorrection__

#define N 240
#include "PowerGrid.h"

using namespace arma;

template<typename T1>
Col<T1> test_FieldCorrection(string testPath)
{
    typedef complex<T1> CxT1;
    //string testPath = "/Users/alexcerjanic/Developer/PG/Resources/";

    //Setup image space coordinates/trajectory
    Mat<T1> ix(N,N);
    Mat<T1> iy(N,N);
    Mat<T1> iz(N,N);
    Mat<T1> FM_gdft(N,N);
    
    ix.zeros();
    iy.zeros();
    iz.zeros();
    FM_gdft.zeros();
    
    //generate the image space cordinates of the voxels we want to reconstruct
    // after vectorizing ix and iy the image coordinates must match the Field and SENSe map image coordinates
    for(uword ii = 0; ii < N; ii++) { //y
        for (uword jj = 0; jj < N; jj++) { //x
            ix(ii,jj) = ((T1)jj-32.0)/(double)N;
            iy(ii,jj) = ((T1)ii-32.0)/(double)N;
        }
    }

    Col<T1> kx;
    loadmat(testPath+"kx.mat","kx",&kx);
    Col<T1> ky;
    loadmat(testPath+"ky.mat","ky",&ky);

    uword nro;
    nro = kx.n_elem;

    Col<T1> kz;
    kz.zeros(nro);

    Col<T1> tvec_gdft;
    tvec_gdft.zeros(nro);
    
    uword L = 10;
    Mat<T1> FM;
    loadmat(testPath+"FM.mat","FM",&FM);

  //FM.zeros();
   T1 tsamp = 5e-6; //sampling rate
   Col<T1> tvec;

   tvec.zeros(nro);
   for (uword ii=0; ii<nro; ii++) {
		   tvec(ii) = ii;
   }
   tvec = tvec*tsamp;

    // Forward operator
    Ggrid<T1> G(nro, 2.0, N, N, 1, kx, ky, kz, vectorise(ix), vectorise(iy), vectorise(iz));

   //cout << "min tvec = " << tvec.min() << endl;

    FieldCorrection<T1, Ggrid<T1>> A(G,vectorise(FM),vectorise(tvec),nro,N*N,L);

    // Variables needed for the recon: Penalty object, num of iterations
    umat ReconMask;
    ReconMask.ones(N,N);

    QuadPenalty<T1> R(N,N, 1, 1, ReconMask);
    
    uword niter = 10;
    Col<CxT1> xinit = zeros<Col<CxT1>>(N*N); // initial estimate of x
    Mat<T1> W;
    W = eye<Mat<T1>>(A.n1,A.n1); // Should be the size of k-space data: Is it right?

    Col<CxT1> data;
    loadmat(testPath+"data_onecoil_FM.mat","data",&data);
    /*
    Col<cx_double> TestAdjoint;
    TestAdjoint = G / data;
    savemat(testPath+"testGgridAdjoint.mat","testGgridAdjoint",TestAdjoint);
    Col<cx_double> TestForward;
    TestForward = G * TestAdjoint;
    savemat(testPath+"testGgridForward.mat","testGgridForward",TestForward);
    TestAdjoint = G / TestForward;
    savemat(testPath+"testGgridAdjoint2.mat","testGgridAdjoint2",TestAdjoint);
*/
    Col<CxT1> x_t;
    x_t = solve_pwls_pcg<T1,FieldCorrection<T1,Ggrid<T1>>,QuadPenalty<T1>>(xinit, A, W, data, R, niter);
    
    return x_t;
    
}
#endif
