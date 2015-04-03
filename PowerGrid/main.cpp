//
//  main.cpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 3/12/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

// Tutorial for linear least square fitting
// Author: Pierre Moulon
// Date: 8 December 2009
// Objective:
// Fit a 2D line with to a set of points
// Specifically, find a and b in the model y = ax + b
//
// Direct application of the example in:
// http://en.wikipedia.org/wiki/Linear_least_squares#Motivational_example


#include "PowerGrid.h"
#include "../Support/CeempleMatio.h"

using namespace arma;
using namespace std;


int main(int argc, char** argv)
{
    
    Mat<cx_double> test;
    
    loadmat("/Users/alexcerjanic/Developer/PowerGrid/Resources/test.mat","test",&test);
    
    Gfft<Col<cx_double>> G(64,64);

    Col<cx_double> testForward = G*vectorise(test);
    
    Col<cx_double> testAdjoint = G/testForward;
    
    savemat("/Users/alexcerjanic/Developer/PowerGrid/Resources/testForwardReal.mat","testForwardReal", real(testForward).eval());
    savemat("/Users/alexcerjanic/Developer/PowerGrid/Resources/testForwardImag.mat","testForwardImag", imag(testForward).eval());

    savemat("/Users/alexcerjanic/Developer/PowerGrid/Resources/testAdjointReal.mat","testAdjointReal", real(testAdjoint).eval());
    savemat("/Users/alexcerjanic/Developer/PowerGrid/Resources/testAdjointImag.mat","testAdjointImag", imag(testAdjoint).eval());


    return 0;
}

