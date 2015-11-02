/*
 * solve_art.hpp
 *
 *  Created on: 2015-10-30
 *      Author: Alex Cerjanic
 */

#ifndef POWERGRID_SOLVE_ART_HPP_
#define POWERGRID_SOLVE_ART_HPP_

#include <cstdlib>

using namespace arma;

template<typename T1, typename Tobj>
Col <complex<T1>> solve_art(const Col <complex<T1>> &xInitial, Tobj const &A, Col <complex<T1>> const &yi,
                            uword niter) {
    typedef complex <T1> CxT1;
    Col <CxT1> dataEst;
    Col <CxT1> imgError;
    uword N = arma::floor(std::sqrt(xInitial.n_row));
    Col <CxT1> img = (1 / (2 * N)) * (A / xInitial);

    for (unsigned int ii = 0; ii < niter; ii++) {
        dataEst = (1 / (2 * N)) * (A * img);
        imgError = (1 / (2 * N)) * (A / (yi - dataEst));
        img = img + imgError;
    }
    return img;
}


#endif /* POWERGRID_SOLVE_ART_HPP_ */
