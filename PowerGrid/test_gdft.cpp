//
//  test_gdft.cpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 4/7/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#include "test_gdft.h"

template<typename T1,typename T2>
T1 test_gdft(const Col<T1>& d,const Col<T2>& kx,const Col<T2>& ky, const Col<T2>& kz)
{
    // Formard operator
    Gdft<cx_double,double> G(64,64,kx,ky,kz);
    Col<cx_double> TestForward;
    TestForward = G * vectorise(conv_to<Mat<cx_double>>::from(test));
    
    // Variables needed for the recon: Penalty object, num of iterations
    umat ReconMask;
    ReconMask.ones(64,64);
    
    QuadPenalty<cx_double>R(64,64,1,ReconMask);
    
    uword niter = 2;
    Col<cx_double> xinit = zeros<Col<cx_double>>(64*64); // initial estimate of x
    Mat<cx_double> W;
    W = eye<Mat<cx_double>>(G.n1,G.n1); // Should be the size of k-space data: Is it right?
    
    Col<cx_double> x_t;
    x_t = pwls_pcg1<cx_double,Gdft<cx_double,double>,QuadPenalty<cx_double>>(xinit, G, W, TestForward, R, niter);
    
    return x_t;

}
