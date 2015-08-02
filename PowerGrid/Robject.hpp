//
// Created by Alex Cerjanic on 7/7/15.
//

#ifndef PowerGrid_Robject_hpp
#define PowerGrid_Robject_hpp

using namespace arma;

template <typename T1>
class Robject
{
public:
    // Robject();
    Robject() {};
    //Class members
    uword Nx;
    uword Ny;
    uword Nz;
    double DeltaX;
    double DeltaY;
    double DeltaZ;
    double Beta;


    // It was declared as type Mat<uword> and the 3D type was a cube. We need to vectorize it before it is passed to QuadPenalty.
    //Custom Class Constructor
    Robject(uword nx,uword ny,uword nz, double beta)
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

    //Class Methods - Declared virtual so they can be implemented in the base classes. Also they are virtual so that if you try to call Robject, things crash rather than give un results.
    virtual Col<T1> wpot(const Col<T1>& d) const;

    virtual Col<T1> pot(const Col<T1>& d) const;

    Col<T1>Cd(const Col<T1>& d, uword dim) const
    {
        /*
        Col<T1> out(Nx*Ny*Nz);
        uword offset;
        for(uword  ll = 0; ll < 2; ll++) {
            for(uword jj =0; jj < 2; jj++) {
                for (uword kk = 0; kk < 2; kk++) {
                    offset = ll + jj*Ny + kk*Nx*Ny;
                    for(uword ii = offset; ii < Ny*Nx*Nz; ii++) {
                        out(ii+((ll+2*jj+4*kk)*(Nx*Ny*Nz))) = d(ii) - d(ii - offset);
                    }
                }
            }
        }
        */
        Col<T1> out(Nx*Ny*Nz);

        uword ll,jj,kk;
        switch (dim)
        {
            case (uword)0:
                ll = 1;
                jj = 0;
                kk = 0;
                break;
            case (uword)1:
                ll = 0;
                jj = 1;
                kk = 0;
                break;
            case (uword)2:
                ll = 0;
                jj = 0;
                kk = 1;
                break;
            default:
                cout << "Warning regularization along dimension greater than 3! Undefined case!" << endl;
        }
        uword offset = ll + jj*Ny + kk*Nx*Ny;
        for(uword ii = offset; ii < Ny*Nx*Nz; ii++) {
            out(ii) = d(ii) - d(ii - offset);
        }

        return out;
    }

    Col<T1>Ctd(const Col<T1>& d, uword dim) const
    {
        /*
        Col<T1> out(Nx*Ny*Nz);
        uword offset;
        for(uword  ll = 0; ll < 2; ll++) {
            for( uword jj =0; jj < 2; jj++) {
                for (uword kk = 0; kk < 2; kk++) {
                    offset = ll + jj*Ny + kk*Nx*Ny;
                    for(uword ii = offset; ii < Ny*Nx*Nz; ii++) {
                        out(ii) += d(ii+((ll+2*jj+4*kk)*(Nx*Ny*Nz)));
                        out(ii-offset) -= d(ii+((ll+2*jj+4*kk)*(Nx*Ny*Nz)));
                    }
                }
            }
        }
        */
        Col<T1> out(Nx*Ny*Nz);

        uword ll,jj,kk;
        switch (dim)
        {
            case (uword)0:
                ll = 1;
                jj = 0;
                kk = 0;
                break;
            case (uword)1:
                ll = 0;
                jj = 1;
                kk = 0;
                break;
            case (uword)2:
                ll = 0;
                jj = 0;
                kk = 1;
                break;
            default:
                cout << "Warning regularization along dimension greater than 3! Undefined case!" << endl;
        }
        uword offset = ll + jj*Ny + kk*Nx*Ny;
        for(uword ii = offset; ii < Ny*Nx*Nz; ii++) {
            out(ii-offset) = d(ii - offset) - d(ii);
        }

        return out;
    }
public:
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

    T1 Denom(const Col<T1>& ddir,const Col<T1>& x) const
    {

        Col <T1> d = zeros< Col < T1 >>(ddir.n_rows);
        Col <T1> Cx = zeros< Col < T1 >>(ddir.n_rows);
        T1 penal = 0;
        T1 temp;
        uword nd = 0;
        if(this->Nz == 1) {
            nd = 2;
        } else {
            nd = 3;
        }

        for(uword ii = 0; ii < nd; ii++) {
            d = this->Cd(ddir,ii);
            Cx = this->wpot(this->Cd(x,ii));
            Cx = Cx % this->Ctd(ddir,ii);
            temp = as_scalar(d*Cx.t());
            penal += temp;
        }

        return (this->DeltaX*this->DeltaY*this->DeltaZ)*this->Beta*penal;
    }



};

#endif //PowerGrid_Robject_hpp
