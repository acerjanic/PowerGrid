/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [solve_grad_desc.hpp]

    Synopsis    [Object implementing a gradient descent solver without line
                    search.]

    Description [This object is only a demonstration of how to use the class
                    objects in PowerGrid. Not recommended for use in practical
                    reconstructions.]

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

/// @file solve_grad_desc.hpp
/// @brief Gradient descent solver (demonstration only — not for production use).

#ifndef POWERGRID_SOLVE_GRAD_DESC_HPP_
#define POWERGRID_SOLVE_GRAD_DESC_HPP_

#include <cstdlib>

using namespace arma;

/// @brief Solve an image reconstruction problem via gradient descent (no line search).
///
/// This is a demonstration solver showing how to use the PowerGrid encoding
/// operators. It is **not recommended** for practical reconstructions; use
/// solve_pwls_pcg instead.
///
/// @tparam T1    Floating-point precision type (`float` or `double`).
/// @tparam Tobj  Encoding operator type (must support `A * x` and `A / y`).
/// @param xInitial  Initial image estimate.
/// @param A         Encoding operator.
/// @param yi        Measured k-space data.
/// @param niter     Number of gradient descent iterations.
/// @returns         Image estimate after @p niter iterations.
template<typename T1, typename Tobj>
Col <complex<T1>> solve_grad_desc(const Col <complex<T1>> &xInitial, Tobj const &A, Col <complex<T1>> const &yi,
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


#endif /* POWERGRID_SOLVE_GRAD_DESC_HPP_ */
