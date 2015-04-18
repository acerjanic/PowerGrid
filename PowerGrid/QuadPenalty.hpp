//
//  QuadPenalty.hpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 4/4/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef PowerGrid_QuadPenalty_hpp
#define PowerGrid_QuadPenalty_hpp

#include <armadillo>

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
    Mat<uword> ReconMask;
    
    mat C1;
    mat C2;
    
    
    //Custom Class Constructor
    QuadPenalty(uword Nx,uword Ny,double Beta,Mat<uword> ReconMask)
    {
        //Set Class Memebers
        this->Nx = Nx;
        this->Ny = Ny;
        this->DeltaX = 1.0/(double)Nx;
        this->DeltaX = 1.0/(double)Ny;
        this->ReconMask = ReconMask;
        
        //Create Sparse Differencing Matricies
        mat DenseX;
        mat DenseY;
        
        //this->C1 = arma::kron(DenseX.eye(Nx,Nx), DenseY.eye(Ny,Ny) - circshift(DenseY.eye(Ny,Ny), 1,0));
        //this->C2 = arma::kron(DenseX.eye(Nx,Nx) - circshift(DenseX.eye(Nx,Nx), 1,0),DenseY.eye(Ny,Ny));
       
    }
    
    //Class Methods


    Col<T1>Cd(const Col<T1>& d) const
    {
       //Only 2D for the moment.
        Col<T1> out(Nx*Ny);
        uword offset;
        for(uword  ll = 0; ll < 2; ll++) {
            for( uword jj =0; jj < 2; jj++) {
                offset = ll + jj*Ny;
                for(uword ii = offset; ii < Ny*Nx; ii++) {
                    out(ii) = d(ii) - d(ii - offset);
                }
            }
        }


        return out;
    }

    double Penalty(const Col<T1>& d) const
    {   Col<T1> x = Cd(d);
        return this->Beta*(pow(this->Cd,2.0))*(this->DeltaX*this->DeltaY);
    }

    Col<T1> Gradient(const Col<T1>& d) const
    {
        return this->Cd(d)*(this->DeltaX*this->DeltaY);
    }


    Col<T1> Denom(const Col<T1>& d) const
    {
        //mat C = join_vert(C1,C2);
        //double weight = conv_to<double>::from((trans(C*d)*(C*d)));
        double weight = 1; //Testing to see how slow this line acutally is.
        return weight*ones<Col<T1>>(d.n_rows);
    }
    
};

#endif
