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

template<typename T1,typename T2>
Col<T1> test_gdft(const Col<T1> d,const Col<T2> kx,const Col<T2> ky, const Col<T2> kz)
{
    // Formard operator
    Gdft<T1,T2> G(64,64,kx,ky,kz);
    
    // Variables needed for the recon: Penalty object, num of iterations
    umat ReconMask;
    ReconMask.ones(64,64);
    
    QuadPenalty<T1>R(64,64,1,ReconMask);
    
    uword niter = 2;
    Col<T1> xinit = zeros<Col<cx_double>>(64*64); // initial estimate of x
    Mat<T1> W;
    W = eye<Mat<T1>>(G.n1,G.n1); // Should be the size of k-space data: Is it right?
    
    Col<T1> x_t;
    x_t = pwls_pcg1<T1,Gdft<T1,T2>,QuadPenalty<T1>>(xinit, G, W, vectorise(d), R, niter);
    
    return x_t;
    
}
#endif /* defined(__PowerGrid__test_gdft__) */
