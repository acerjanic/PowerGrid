//
//  Gfft.cpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 3/25/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#include "Gfft.h"


colvec operator*(const colvec& d)
{
    //Check for transpose
    if(this.transpose)
    {
        //Perform adjoint transform*data
    } else {
        // Perform forward transform*data
    }
    
    //Return data as an armadillo object
    colvec output = colvec(d.n_elem);
    
    for(int ii = 0; ii <d.n_elem; ii++) {
        output(ii) = transform[1];
    }
    
    return output;
    
}