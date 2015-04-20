//
//  test_gdft.h
//  PowerGrid
//
//  Created by Alex Cerjanic on 4/7/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef __PowerGrid__test_gdft__
#define __PowerGrid__test_gdft__

#include "PowerGrid.h"

using namespace arma;
//We use the two objects to capture the complex<T2> as T1 and T2 as the real (double or float) data type. This can be simplified a bit...
template<typename T1,typename T2>
Col<T1> test_gdft(const Col<T1> d,const Col<T2> kx,const Col<T2> ky, const Col<T2> kz)
{
	 string testPath = "/shared/mrfil-data/jholtrop/repos/PowerGrid/Resources/";

    //Setup image space coordinates/trajectory
    Mat<T2> ix(64,64);
    Mat<T2> iy(64,64);
    Mat<T2> iz(64,64);
    ix.zeros();
    iy.zeros();
    iz.zeros();

    Mat<T2> FM;
    loadmat(testPath+"FM.mat","FM",&FM);
    //FM.zeros();
   T2 tsamp = 5e-6;
   Col<T2> t;
   t.zeros(kx.n_elem);
   for (uword ii=0; ii<kx.n_elem; ii++) {
		   t(ii) = ii;
   }
   t = t*tsamp;



    for(uword ii = 0; ii < 64; ii++) {
        for (uword jj = 0; jj < 64; jj++) {
            ix(ii,jj) = ((T2)jj-32.0)/64.0;
            iy(ii,jj) = ((T2)ii-32.0)/64.0;
        }
    }
    //savemat(testPath+"ix.mat","ix",ix);
    //savemat(testPath+"iy.mat","iy",iy);
    // Forward operator
    Gdft<T1,T2> G(4010,64*64,kx,ky,kz,vectorise(ix),vectorise(iy),vectorise(iz),vectorise(FM),vectorise(t));
    
    Col<T1> data;
    loadmat(testPath+"data_onecoil_FM.mat","data",&data);

    Col<cx_double> TestForward;
    TestForward = G * vectorise(conv_to<Mat<cx_double>>::from(d));
    savemat(testPath+"testGdftForward.mat","testGdftForward",TestForward);
    
    Col<cx_double> TestAdjoint;
    //TestAdjoint = G / TestForward;
    TestAdjoint = G / data;
    savemat(testPath+"testGdftAdjoint.mat","testGdftAdjoint",TestAdjoint);


    // Variables needed for the recon: Penalty object, num of iterations
    //Recon mask is used to limit the support of the object. Very useful for SENSE recons.
    umat ReconMask;
    ReconMask.ones(64,64);

    QuadPenalty<T1> R(64, 64, 1, 1, ReconMask);
    
    uword niter = 10;
    Col<T1> xinit = zeros<Col<cx_double>>(64*64); // initial estimate of x
    Mat<T1> W;
    W = eye<Mat<T1>>(G.n2,G.n2); // W is the weighting matrix. Sometimes we might want to weight the data inside of pwls,
    //  we don't really use it much at the moment.

    Col<T1> x_t;
    //We call pwls_pcg as our nonlinear solver to minmize our cost function to find the image.
    //pwls_pcg has an implict cost function of C = G'(y - G*x) + R(x), where G is the forward transform, G' the adjoint transform, R(x) (QuadPenalty for now) the regularization penalty for image roughness,
    // y is the measured data, and x the image.
    x_t = pwls_pcg1<T1,Gdft<T1,T2>,QuadPenalty<T1>>(xinit, G, W, TestForward, R, niter);
    
    return x_t;
    
}
#endif /* defined(__PowerGrid__test_gdft__) */
