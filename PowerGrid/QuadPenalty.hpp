//
//  QuadPenalty.hpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 4/4/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef PowerGrid_QuadPenalty_hpp
#define PowerGrid_QuadPenalty_hpp

using namespace arma;

template <typename T1>
class QuadPenalty
{
public:
    QuadPenalty();
    
    //Class members
    uword Nx;
    uword Ny;
    double DeltaX;
    double DeltaY;
    double Beta;
    Col<uint_8> ReconMask;
    
    //Custom Class Constructor
    QuadPenalty(uword Nx,uword Ny,double Beta,Col<uint_8> ReconMask)
    {
        //Set Class Memebers
        this->Nx = Nx;
        this->Ny = Ny;
        this->DeltaX = 1.0/(double)Nx;
        this->DeltaX = 1.0/(double)Ny;
        this->ReconMask = ReconMask;
        
        //Create Sparse Differencing Matricies
        sp_mat D = speye<sp_mat>(Nx,Ny);
        
        
    }
    
    //Class Methods
    
};

#endif
