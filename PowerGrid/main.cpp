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


#include <iostream>
#include "armadillo"

namespace arma {
#include "op_circshift_bones.hpp"
#include "op_circshift_meat.hpp"
#include "fn_circshift.hpp"
}

using namespace arma;
using namespace std;


int main(int argc, char** argv)
{
    
    mat A = randu<mat>(4,5);
    mat B = randu<mat>(4,5);
    
    cout << A << endl;
    
    cout << circshift(A,1,0) << endl;
    
    cout << A*B.t() << endl;
    
    return 0;
}

