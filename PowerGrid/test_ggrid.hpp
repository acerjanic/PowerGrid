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

template<typename T1>
Col<T1> test_ggrid(const Col<complex<T1>> d,const Col<T1> kx,const Col<T1> ky, const Col<T1> kz)
{
    typedef complex<T1> CxT1;
    string testPath = "/shared/mrfil-data/jholtrop/repos/PowerGrid/Resources/";
    cout << "Entered test_ggrid" << endl;
    //Setup image space coordinates/trajectory
    Mat<T1> ix(64,64);
    Mat<T1> iy(64,64);
    Mat<T1> iz(64,64);
    Mat<T1> FM(64,64);
    Mat<T1> t(64,64);
    ix.zeros();
    iy.zeros();
    iz.zeros();
    FM.zeros();
    t.zeros();
    for(uword ii = 0; ii < 64; ii++) {
        for (uword jj = 0; jj < 64; jj++) {
            ix(ii,jj) = ((T1)ii-32.0)/64.0;
            iy(ii,jj) = ((T1)jj-32.0)/64.0;
        }
    }
    cout << "Saving some data" << endl;
    savemat(testPath+"ix.mat","ix",ix);
    savemat(testPath+"iy.mat","iy",iy);
    // Forward operator
    Ggrid<T1> G(4010,64,64,1,kx,ky,kz,vectorise(ix),vectorise(iy),vectorise(iz));
    cout << "About to run forward transform" << endl;
    Col<CxT1> TestForward;
    TestForward = G * vectorise(conv_to<Mat<cx_double>>::from(d));
    savemat(testPath+"testGgridForward.mat","testGgridForward",TestForward);
    cout << "About to run adjoint transform" << endl;

    Col<CxT1> TestAdjoint;
    TestAdjoint = G / TestForward;
    savemat(testPath+"testGgridAdjoint.mat","testGgridAdjoint",TestAdjoint);



    // Variables needed for the recon: Penalty object, num of iterations
    umat ReconMask;
    ReconMask.ones(64,64);

    QuadPenalty<T1> R(64, 64, 1, 1, ReconMask);
    cout << "Just made the penalty object" << endl;

    uword niter = 20;
    Col<CxT1> xinit = zeros<Col<cx_double>>(64*64); // initial estimate of x
    Col<T1> W;
    W = ones(G.n2); // Should be the size of k-space data: Is it right?
    cout << "About to call solve_pwls_pcg from test_ggrid" << endl;
    Col<CxT1> x_t;
    x_t = solve_pwls_pcg<T1,Ggrid<T1>,QuadPenalty<T1>>(xinit, G, W, TestForward, R, niter);
    
    return x_t;
    
}
#endif /* defined(__PowerGrid__test_ggrid__) */
