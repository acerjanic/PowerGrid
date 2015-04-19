//
//  test_FieldCorrection.h
//  PowerGrid
//
//  Created by Joe holtrop on 4/7/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef __PowerGrid__test_FieldCorrection__
#define __PowerGrid__test_FieldCorrection__

#include "PowerGrid.h"

using namespace arma;

template<typename T1,typename T2>
Col<T1> test_FieldCorrection()
{
	 string testPath = "/Users/alexcerjanic/Developer/PG/Resources/";

    //Setup image space coordinates/trajectory
    Mat<T2> ix(64,64);
    Mat<T2> iy(64,64);
    Mat<T2> iz(64,64);
    Mat<T2> FM_gdft(64,64);
    
    ix.zeros();
    iy.zeros();
    iz.zeros();
    FM_gdft.zeros();
    
    //generate the image space cordinates of the voxels we want to reconstruct
    // after vectorizing ix and iy the image coordinates must match the Field and SENSe map image coordinates
    for(uword ii = 0; ii < 64; ii++) { //y
        for (uword jj = 0; jj < 64; jj++) { //x
            ix(ii,jj) = ((T2)jj-32.0)/64.0;
            iy(ii,jj) = ((T2)ii-32.0)/64.0;
        }
    }

    Col<T2> kx;
    loadmat(testPath+"kx_sp.mat","kx",&kx);
    Col<T2> ky;
    loadmat(testPath+"ky_sp.mat","ky",&ky);

    uword nro;
    nro = kx.n_elem;

    Col<T2> kz;
    kz.zeros(nro);

    Col<T2> tvec_gdft;
    tvec_gdft.zeros(nro);
    
    uword L = 10;
    Mat<T2> FM;
    loadmat(testPath+"FM.mat","FM",&FM);

  //FM.zeros();
   T2 tsamp = 5e-6; //sampling rate
   Col<T2> tvec;

   tvec.zeros(nro);
   for (uword ii=0; ii<nro; ii++) {
		   tvec(ii) = ii;
   }
   tvec = tvec*tsamp;

    // Forward operator
    Ggrid<T1, T2> G(nro, 4.0, 64, 64, 1, kx, ky, kz, vectorise(ix), vectorise(iy), vectorise(iz));

   //cout << "min tvec = " << tvec.min() << endl;

    FieldCorrection<T1, T2, Ggrid<T1,T2>> A(G,vectorise(FM),vectorise(tvec),nro,64*64,L);

    // Variables needed for the recon: Penalty object, num of iterations
    umat ReconMask;
    ReconMask.ones(64,64);

    QuadPenalty<T1> R(64, 64, 1, 0, ReconMask);
    
    uword niter = 10;
    Col<T1> xinit = zeros<Col<cx_double>>(64*64); // initial estimate of x
    Mat<T1> W;
    W = eye<Mat<T1>>(A.n1,A.n1); // Should be the size of k-space data: Is it right?

    Col<T1> data;
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
    Col<T1> x_t;
    x_t = pwls_pcg1<T1,FieldCorrection<T1,T2,Ggrid<T1,T2>>,QuadPenalty<T1>>(xinit, A, W, data, R, niter);
    
    return x_t;
    
}
#endif
