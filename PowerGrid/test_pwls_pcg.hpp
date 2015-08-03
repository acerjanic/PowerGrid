//
//  test_FieldCorrection.h
//  PowerGrid
//
//  Created by Joe Holtrop on 4/13/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef PowerGrid_test_FieldCorrection_h
#define PowerGrid_test_FieldCorrection_h


template<typename T1>
Mat<std::complex<T1>> test_pwls_pcg()
{
    typedef std::complex<T1> CxT1;
    std::string testPath = "/shared/mrfil-data/jholtrop/repos/PowerGrid/Resources/";

    // Load real data
    Mat<T1> test;

    //Load our read double matrix object. (You need to match type to avoide mangling the data.)
    loadmat(testPath+"test.mat","test",&test);

    // Formard operator
    Gfft<T1> G(64,64);
    Col<CxT1> TestForward;
    TestForward = G *vectorise(conv_to<Mat<cx_double>>::from(test));

    // Variables needed for the recon: Penalty object, num of iterations
    umat ReconMask;
    ReconMask.ones(64,64);

    QuadPenalty<T1> R(64, 64, 1, 1, ReconMask);

    uword niter = 2;
    Col<T1> xinit = zeros<Col<CxT1>>(64*64); // initial estimate of xd
    double W;
    W = 1.0;
    //Mat<cx_double> W;
    //W = eye<Mat<cx_double>>(G.n1*G.n1,G.n1*G.n1); // Should be the size of k-space data: Is it right?

    Col<CxT1> x_t;
    x_t = solve_pwls_pcg<T1,Gfft<T1>,QuadPenalty<T1>>(xinit, G, W, TestForward, R, niter);

    return x_t;
}

#endif
