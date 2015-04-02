//
//  Gfft.cpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 3/25/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//
#ifndef __PowerGrid__Gfft__hpp
#define __PowerGrid__Gfft__hpp

#include "Gfft.h"

//Forward transform operation
cx_colvec operator*(const cx_colvec& d)
{
    
    //Return data as an armadillo object
    colvec output = d;
    
    
    //equivalent to returning col(output) in MATLAB with IRT
    return output.vectorise();
    
}

//Adjoint transform operation
cx_colvec operator/(const cx_colvec& d)
{
    //Return data as an armadillo object
    cx_colvec output = d;
    
    Mat<cx_double> d2D =
    //equivalent to returning col(output) in MATLAB with IRT
    return output.vectorise();
    
}
#endif