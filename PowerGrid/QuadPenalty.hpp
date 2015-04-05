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
    Col<uword> ReconMask;
    
    sp_mat C1;
    sp_mat C2;
    
    
    //Custom Class Constructor
    QuadPenalty(uword Nx,uword Ny,double Beta,Col<uword> ReconMask)
    {
        //Set Class Memebers
        this->Nx = Nx;
        this->Ny = Ny;
        this->DeltaX = 1.0/(double)Nx;
        this->DeltaX = 1.0/(double)Ny;
        this->ReconMask = ReconMask;
        
        //Create Sparse Differencing Matricies
        this->C1 = arma::kron(speye<sp_mat>(Nx,Nx), speye(Ny,Ny) - circshift(speye(Ny,Ny), 1,0));
        this->C2 = arma::kron(speye<sp_mat>(Nx,Nx) - circshift(speye(Nx,Nx), 1,0),speye(Ny,Ny));
       
    }
    
    //Class Methods
    double Penalty(const T1& d)
    {
        return this->Beta*(pow(C1*d,2.0)+pow(C2*d,2.0))*(this->DeltaX*this->DeltaY);
    }
    
    T1 Gradient(const T1& d)
    {
        return (C1*d+C2*d)*(this->DeltaX*this->DeltaY);
    }
    
    
    T1 Denom(const T1& d)
    {
        sp_mat C = join_vert(C1,C2);
        return trans(C*d)*(C*d);
    }
    
};

#endif
