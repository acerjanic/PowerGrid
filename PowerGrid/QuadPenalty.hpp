//
//  QuadPenalty.hpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 4/4/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef PowerGrid_QuadPenalty_hpp
#define PowerGrid_QuadPenalty_hpp

#include "armadillo"

using namespace arma;

template <typename T1>
class QuadPenalty
{
public:
    QuadPenalty();

    //Class members
    uword Nx;
    uword Ny;
    uword Nz;
    double DeltaX;
    double DeltaY;
    double DeltaZ;
    double Beta;
    Mat<uword> ReconMask;

    mat C1;
    mat C2;


    //Custom Class Constructor
    QuadPenalty(uword Nx,uword Ny,uword Nz, double Beta,Mat<uword> ReconMask)
    {
        //Set Class Memebers
        this->Nx = Nx;
        this->Ny = Ny;
        this->Nz = Nz;
        this->DeltaX = 1.0/(double)Nx;
        this->DeltaY = 1.0/(double)Ny;
        this->DeltaZ = 1.0/(double)Nz;

        this->ReconMask = ReconMask;

    }

    //Class Methods


    Col<T1>Cd(const Col<T1>& d) const
    {

        Col<T1> out(Nx*Ny*Nz*8);
        uword offset;
        for(uword  ll = 0; ll < 2; ll++) {
            for( uword jj =0; jj < 2; jj++) {
                for (uword kk = 0; kk < 2; kk++) {
                    offset = ll + jj*Nx + kk*Nx*Ny;
                    for(uword ii = offset; ii < Ny*Nx*Nz; ii++) {
                        out(ii+((ll+2*jj+4*kk)*(Nx*Ny*Nz))) = d(ii) - d(ii - offset);
                    }
                }
            }
        }


        return out;
    }

    Col<T1>Ct(const Col<T1>& d) const
    {

        Col<T1> out(Nx*Ny*Nz);
        uword offset;
        for(uword  ll = 0; ll < 2; ll++) {
            for( uword jj =0; jj < 2; jj++) {
                for (uword kk = 0; kk < 2; kk++) {
                    offset = ll + jj*Nx + kk*Nx*Ny;
                    for(uword ii = offset; ii < Ny*Nx*Nz; ii++) {
                        out(ii) += d(ii+((ll+2*jj+4*kk)*(Nx*Ny*Nz)));
                        out(ii-offset) -= d(ii+((ll+2*jj+4*kk)*(Nx*Ny*Nz)));
                    }
                }
            }
        }


        return out;
    }


    double Penalty(const Col<T1>& d) const
    {   
        return Beta*(pow(Cd(d),2.0))*(DeltaX*DeltaY*DeltaZ);
    }

    Col<T1> Gradient(const Col<T1>& d) const
    {
       // THIS IS NOT RIGHT: NEEDS HELP!
        return Beta*Ct(Cd(d))*(DeltaX*DeltaY*DeltaZ);
    }


    Col<T1> Denom(const Col<T1>& d) const
    {
        //mat C = join_vert(C1,C2);
        //double weight = conv_to<double>::from((trans(C*d)*(C*d)));
        double weight = 1.0/d.n_rows; //Testing to see how slow this line acutally is.
        return weight*ones<Col<T1>>(d.n_rows);
    }

};

#endif
