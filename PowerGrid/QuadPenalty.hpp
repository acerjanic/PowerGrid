//
//  QuadPenalty.hpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 4/4/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef PowerGrid_QuadPenalty_hpp
#define PowerGrid_QuadPenalty_hpp

//#include <armadillo>

using namespace arma;

template <typename T1>
class QuadPenalty: public Robject<T1>
{
public:
    QuadPenalty();

    // It was declared as type Mat<uword> and the 3D type was a cube. We need to vectorize it before it is passed to QuadPenalty.
    //Custom Class Constructor
    QuadPenalty(uword nx,uword ny,uword nz, double beta)
    {
        //Set Class Memebers
        this->Nx = nx;
        this->Ny = ny;
        this->Nz = nz;
        this->DeltaX = 1.0/(double)nx;
        this->DeltaY = 1.0/(double)ny;
        this->DeltaZ = 1.0/(double)nz;
        this->Beta = beta;

    }

    //Class Methods

    Col<T1> wpot(const Col<T1>& d) const
    {
        return ones < Col < T1 >> (d.n_rows);
    }

    Col<T1> pot(const Col<T1>& d) const
    {
        return d % d / 2.0;
    }


};

#endif
