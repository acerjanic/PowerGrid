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
public:
    TVPenalty();

    // It was declared as type Mat<uword> and the 3D type was a cube. We need to vectorize it before it is passed to QuadPenalty.
    //Custom Class Constructor
    TVPenalty(uword nx,uword ny,uword nz, double beta, double delta)
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

    Col<T1> wpot(const Col<T1>& d) const
    {   Col<double> temp = abs(d / this->Delta);
        Col<T1> out = 1.0/sqrt(1.0 + conv_to<Col<T1>>::from(temp % temp));
        return out;
    }

    Col<T1> pot(const Col<T1>& d) const
    {
        Col <T1> out = this->delta * this->delta * (sqrt(1.0 + (abs(d / this->Delta) % abs(d / this->Delta))) - 1.0);
        return out;
    }



private:
    double Delta;
};

#endif //PowerGrid_TVPenalty_hpp
