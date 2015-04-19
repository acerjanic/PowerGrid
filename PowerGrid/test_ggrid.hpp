//
//  test_gdft.h
//  PowerGrid
//
//  Created by Alex Cerjanic on 4/7/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef __PowerGrid__test_ggrid__
#define __PowerGrid__test_ggrid__

#include "PowerGrid.h"

using namespace arma;
// Test function for Ggrid that uses gridding on the CPU (for now) to approximate the forward and adjoint Fourier Transforms

template<typename T1,typename T2>
Col<T1> test_ggrid(const Col<T1> d,const Col<T2> kx,const Col<T2> ky, const Col<T2> kz)
{
    string testPath = "/shared/mrfil-data/jholtrop/repos/PowerGrid/Resources/";
    cout << "Entered test_ggrid" << endl;
    //Setup image space coordinates/trajectory
    Mat<T2> ix(64,64);
    Mat<T2> iy(64,64);
    Mat<T2> iz(64,64);
    Mat<T2> FM(64,64);
    Mat<T2> t(64,64);
    ix.zeros();
    iy.zeros();
    iz.zeros();
    FM.zeros();
    t.zeros();
    for(uword ii = 0; ii < 64; ii++) {
        for (uword jj = 0; jj < 64; jj++) {
            ix(ii,jj) = ((T2)ii-32.0)/64.0;
            iy(ii,jj) = ((T2)jj-32.0)/64.0;
        }
    }
    cout << "Saving some data" << endl;
    savemat(testPath+"ix.mat","ix",ix);
    savemat(testPath+"iy.mat","iy",iy);
    // Forward operator
    Ggrid<T1,T2> G(4010,64,64,1,kx,ky,kz,vectorise(ix),vectorise(iy),vectorise(iz));
    cout << "About to run forward transform" << endl;
    Col<cx_double> TestForward;
    TestForward = G * vectorise(conv_to<Mat<cx_double>>::from(d));
    savemat(testPath+"testGgridForward.mat","testGgridForward",TestForward);
    cout << "About to run adjoint transform" << endl;

    Col<cx_double> TestAdjoint;
    TestAdjoint = G / TestForward;
    savemat(testPath+"testGgridAdjoint.mat","testGgridAdjoint",TestAdjoint);



    // Variables needed for the recon: Penalty object, num of iterations
    umat ReconMask;
    ReconMask.ones(64,64);

    QuadPenalty<T1> R(64, 64, 1, 1, ReconMask);
    cout << "Just made the penalty object" << endl;

    uword niter = 20;
    Col<T1> xinit = zeros<Col<cx_double>>(64*64); // initial estimate of x
    Mat<T1> W;
    W = eye<Mat<T1>>(G.n2,G.n2); // Should be the size of k-space data: Is it right?
    cout << "About to call pwls_pcg1 from test_ggrid" << endl;
    Col<T1> x_t;
    x_t = pwls_pcg1<T1,Ggrid<T1,T2>,QuadPenalty<T1>>(xinit, G, W, TestForward, R, niter);
    
    return x_t;
    
}
#endif /* defined(__PowerGrid__test_ggrid__) */
