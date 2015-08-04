//
//  TVPenalty.hpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 6/8/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef PowerGrid_TVPenalty_hpp
#define PowerGrid_TVPenalty_hpp

//#include <armadillo>

using namespace arma;

template <typename T1>
class TVPenalty: public Robject<T1>
{
typedef complex<T1> CxT1;
public:
    TVPenalty();

    // It was declared as type Mat<uword> and the 3D type was a cube. We need to vectorize it before it is passed to QuadPenalty.
    //Custom Class Constructor
    TVPenalty(uword nx,uword ny,uword nz, T1 beta, T1 delta)
    {
        //Set Class Members
        this->Nx = nx;
        this->Ny = ny;
        this->Nz = nz;
        this->DeltaX = 1.0/(double)nx;
        this->DeltaY = 1.0/(double)ny;
        this->DeltaZ = 1.0/(double)nz;
        this->Beta = beta;
        this->Delta = delta;

    }

    //Class Methods

    Col<CxT1> wpot(const Col<CxT1>& d) const
    {   Col<T1> temp = abs(d / this->Delta);
        Col<CxT1> out = 1.0/sqrt(1.0 + conv_to<Col<CxT1>>::from(temp % temp));
        return out;
    }

    Col<CxT1> pot(const Col<CxT1>& d) const
    {
        Col <CxT1> out = conv_to<Col<CxT1>>::from(this->Delta * this->Delta * (sqrt(1.0 + (abs(d / this->Delta) % abs(d / this->Delta))) - 1.0));
        return out;
    }



private:
    T1 Delta;
};

#endif //PowerGrid_TVPenalty_hpp
