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
class QuadPenalty: protected Robject<T1>
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

    double Penalty(const Col<T1>& x) const
    {
        Col <T1> d = zeros< Col < T1 >>(x.n_rows);
        double penal = 0;
        uword nd = 0;
        if(this->Nz == 1) {
            nd = 2;
        } else {
            nd = 3;
            cout << "Setting dimension to 3 in reg." << endl;
        }

        for(uword ii = 0; ii < nd; ii++) {
            d = this->Cd(x,ii);
            d = this->pot(d);

            penal = penal + abs(sum(d));
        }

        return (this->DeltaX*this->DeltaY*this->DeltaZ)*this->Beta*penal;
    }

    Col<T1> Gradient(const Col<T1>& x) const
    {
        Col <T1> g = zeros< Col < T1 >>(x.n_rows);
        Col <T1> d = zeros< Col < T1 >>(x.n_rows);
        uword nd = 0;
        if(this->Nz == 1) {
            nd = 2;
        } else {
            nd = 3;
            cout << "Setting dimension to 3 in reg." << endl;

        }

        for(uword ii = 0; ii < nd; ii++) {
            d = this->Cd(x,ii);
            d = this->wpot(d);
            d = this->Ctd(d,ii);
            g = g + d;
        }
        return this->Beta*g*(this->DeltaX*this->DeltaY*this->DeltaZ);
    }


    double Denom(const Col<T1>& ddir) const
    {

        Col <T1> d = zeros< Col < T1 >>(ddir.n_rows);
        double penal = 0;
        uword nd = 0;
        if(this->Nz == 1) {
            nd = 2;
        } else {
            nd = 3;
        }

        for(uword ii = 0; ii < nd; ii++) {
            d = this->Cd(ddir,ii);
            d = this->wpot(d);

            penal = penal + abs(sum(d));
        }

        return (this->DeltaX*this->DeltaY*this->DeltaZ)*this->Beta*penal;
    }

};

#endif
