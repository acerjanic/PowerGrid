/*
   (C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
   All rights reserved.

   See LICENSE.txt for the University of Illinois/NCSA Open Source license.

   Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
 */

/*****************************************************************************

    File Name   [Robject.cpp]

    Synopsis    [Base class of implementing quadratic and non-quadratic
                    regularization.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

*****************************************************************************/

#include "Robject.h"
#include <chrono>  // for high_resolution_clock

// template <typename T1> Robject<T1>::Robject(){};

// It was declared as type Mat<uword> and the 3D type was a cube. We need to
// vectorize it before it is passed to QuadPenalty.
// Custom Class Constructor
template <typename T1>
Robject<T1>::Robject(uword nx, uword ny, uword nz, T1 beta, uword dims2penalize) {
        // Set Class Memebers
        this->Nx = nx;
        this->Ny = ny;
        this->Nz = nz;
        this->Beta = beta;
        this->Dims2Penalize = dims2penalize;
}

// Class Methods - Declared virtual so they can be implemented in the base
// classes. Also they are virtual so that if you try to call Robject, things
// crash rather than give un results.

template <typename T1>
Col<complex<T1>> Robject<T1>::Cd(const Col<complex<T1> > &d, uword dim) const {

        Col<complex<T1> > out(Nx * Ny * Nz);
        out.zeros();
        uword ll, jj, kk;
        switch (dim) {
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
                cout << "Warning regularization along dimension greater than 3! "
                        "Undefined case!"
                     << endl;
        }
        uword offset = ll + jj * Ny + kk * Nx * Ny;
        for (uword ii = offset; ii < Ny * Nx * Nz; ii++) {
                out(ii) = d(ii) - d(ii - offset);
        }

        return out;
}

template <typename T1>
Col<complex<T1>> Robject<T1>::Ctd(const Col<complex<T1> > &d, uword dim) const {
        Col<complex<T1> > out(Nx * Ny * Nz);
        out.zeros();
        uword ll, jj, kk;
        switch (dim) {
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
                cout << "Warning regularization along dimension greater than 3! "
                        "Undefined case!"
                     << endl;
        }
        
        uword n = Ny * Nx * Nz;
        uword offset = ll + jj * Ny + kk * Nx * Ny;

        for (uword i = 0; i < offset; i++) {
                out(i) = -d(i + offset);
        }
        for (uword i = offset; i < n - offset; i++) {
                out(i) = d(i) - d(i + offset);
        }
        for (uword i = n - offset; i < n; i++) {
                out(i) = d(i);
        }

        return out;
}

template <typename T1>
T1 Robject<T1>::Penalty(const Col<complex<T1> > &x) const {
        RANGE()
        Col<complex<T1> > d = zeros<Col<complex<T1> > >(x.n_rows);
        T1 penal = 0;
        uword nd = 0;
        if ((this->Nz == 1) || (this->Dims2Penalize == 2)) {
                nd = 2;
        } else {
                nd = 3;
                // cout << "Setting dimension to 3 in reg." << endl;
        }

        for (uword ii = 0; ii < nd; ii++) {
                d = this->Cd(x, ii);
                d = this->pot(d);
                penal = penal + abs(sum(d));
        }

        return this->Beta * penal;
}
template <typename T1>
Col<complex<T1>> Robject<T1>::Gradient(const Col<complex<T1> > &x) const {
        RANGE()
        Col<complex<T1> > g = zeros<Col<complex<T1> > >(x.n_rows);
        Col<complex<T1> > d = zeros<Col<complex<T1> > >(x.n_rows);
        uword nd = 0;
        if ((this->Nz == 1) || (this->Dims2Penalize == 2)){
                nd = 2;
        } else {
                nd = 3;
                // cout << "Setting dimension to 3 in reg." << endl;
        }

        for (uword ii = 0; ii < nd; ii++) {
                d = this->Cd(x, ii);
                d = this->dpot(d);
                d = this->Ctd(d, ii);
                g += d;

        }

        return this->Beta * g;
}

template <typename T1>
complex<T1> Robject<T1>::Denom(const Col<complex<T1> > &ddir,
                               const Col<complex<T1> > &x) const {
        RANGE()

        Col<complex<T1>> Cdir = zeros<Col<complex<T1>>>(ddir.n_rows);
        Col<complex<T1>> Cx = zeros<Col<complex<T1>>>(ddir.n_rows);
        complex<T1> penal = 0;
        complex<T1> temp;
        complex<T1> cxBeta(this->Beta, 0);
        uword nd = 0;
        if ((this->Nz == 1) || (this->Dims2Penalize == 2)){
                nd = 2;
        } else {
                nd = 3;
        }
        for (uword ii = 0; ii < nd; ii++) {
                Cdir = this->Cd(ddir, ii);
                Cx = this->wpot(this->Cd(x, ii));
                Cx %= Cdir;
                temp = as_scalar(cdot(Cdir,Cx));
                penal += temp;

        }

        return penal * cxBeta;
}

template <typename T1>
pgCol<pgComplex<T1>> Robject<T1>::Cd(const pgCol<pgComplex<T1>>& d, arma::uword dim) const {
    pgCol<pgComplex<T1>> out(Nx * Ny * Nz);
    out.zeros();
    arma::uword ll, jj, kk;
    switch (dim) {
    case 0: ll = 1; jj = 0; kk = 0; break;
    case 1: ll = 0; jj = 1; kk = 0; break;
    case 2: ll = 0; jj = 0; kk = 1; break;
    default:
        std::cout << "Warning: regularization along dimension > 3 undefined!" << std::endl;
        return out;
    }
    arma::uword offset = ll + jj * Ny + kk * Nx * Ny;
    for (arma::uword ii = offset; ii < Ny * Nx * Nz; ii++) {
        out.at(ii) = d.at(ii) - d.at(ii - offset);
    }
    return out;
}

template <typename T1>
pgCol<pgComplex<T1>> Robject<T1>::Ctd(const pgCol<pgComplex<T1>>& d, arma::uword dim) const {
    pgCol<pgComplex<T1>> out(Nx * Ny * Nz);
    out.zeros();
    arma::uword ll, jj, kk;
    switch (dim) {
    case 0: ll = 1; jj = 0; kk = 0; break;
    case 1: ll = 0; jj = 1; kk = 0; break;
    case 2: ll = 0; jj = 0; kk = 1; break;
    default:
        std::cout << "Warning: regularization along dimension > 3 undefined!" << std::endl;
        return out;
    }
    arma::uword n = Ny * Nx * Nz;
    arma::uword offset = ll + jj * Ny + kk * Nx * Ny;
    for (arma::uword i = 0; i < offset; i++) {
        out.at(i) = pgComplex<T1>(T1(0), T1(0)) - d.at(i + offset);
    }
    for (arma::uword i = offset; i < n - offset; i++) {
        out.at(i) = d.at(i) - d.at(i + offset);
    }
    for (arma::uword i = n - offset; i < n; i++) {
        out.at(i) = d.at(i);
    }
    return out;
}

template <typename T1>
T1 Robject<T1>::Penalty(const pgCol<pgComplex<T1>>& x) const {
    RANGE()
    T1 penal = T1(0);
    arma::uword nd = ((this->Nz == 1) || (this->Dims2Penalize == 2)) ? 2 : 3;
    for (arma::uword ii = 0; ii < nd; ii++) {
        pgCol<pgComplex<T1>> d = this->Cd(x, ii);
        d = this->pot(d);
        pgComplex<T1> s = sum(d);
        penal += std::abs(s.real());
    }
    return this->Beta * penal;
}

template <typename T1>
pgCol<pgComplex<T1>> Robject<T1>::Gradient(const pgCol<pgComplex<T1>>& x) const {
    RANGE()
    pgCol<pgComplex<T1>> g(x.n_elem);
    g.zeros();
    arma::uword nd = ((this->Nz == 1) || (this->Dims2Penalize == 2)) ? 2 : 3;
    for (arma::uword ii = 0; ii < nd; ii++) {
        pgCol<pgComplex<T1>> d = this->Cd(x, ii);
        d = this->dpot(d);
        d = this->Ctd(d, ii);
        g += d;
    }
    g %= pgComplex<T1>(this->Beta, T1(0));
    return g;
}

template <typename T1>
pgComplex<T1> Robject<T1>::Denom(const pgCol<pgComplex<T1>>& ddir,
                                  const pgCol<pgComplex<T1>>& x) const {
    RANGE()
    pgComplex<T1> penal(T1(0), T1(0));
    pgComplex<T1> cxBeta(this->Beta, T1(0));
    arma::uword nd = ((this->Nz == 1) || (this->Dims2Penalize == 2)) ? 2 : 3;
    for (arma::uword ii = 0; ii < nd; ii++) {
        pgCol<pgComplex<T1>> Cdir = this->Cd(ddir, ii);
        pgCol<pgComplex<T1>> Cx = this->wpot(this->Cd(x, ii));
        Cx %= Cdir;
        pgComplex<T1> temp = cdot(Cdir, Cx);
        penal += temp;
    }
    return penal * cxBeta;
}

// Explicit Instantiation
template class Robject<float>;
template class Robject<double>;
