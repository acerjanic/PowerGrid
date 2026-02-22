/// @file fftshift.hpp
/// @brief Circular FFT-shift utilities for 1-D and 2-D Armadillo arrays.

//
//  fftshift.hpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 4/2/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef __PowerGrid__fftshift__hpp
#define __PowerGrid__fftshift__hpp
#include "armadillo"


using namespace arma;
using namespace std;

/// @brief Circular shift of a 1-D or 2-D array along a single dimension.
///
/// Shifts @p X by floor(size/2) positions along the given dimension, implementing
/// the standard FFT-shift for centring the zero-frequency component.
///
/// @tparam T1   Armadillo array type (e.g., `Col<cx_float>`, `Mat<float>`).
/// @param X     Input array (modified in-place via circshift copy).
/// @param dim   Dimension to shift: 0 = rows, 1 = columns.
/// @returns     Shifted copy of @p X.
template<typename T1>
arma_inline
T1 fftshift(T1& X, uword dim);

/// @brief 2-D FFT-shift: shift both rows and columns by floor(size/2).
///
/// Equivalent to calling `fftshift(X, 0)` then `fftshift(X, 1)`.
///
/// @tparam T1   Armadillo array type.
/// @param X     Input 2-D array.
/// @returns     2-D shifted copy of @p X.
template<typename T1>
arma_inline
T1 fftshift(T1 X);

template<typename T1>
arma_inline
T1 fftshift(
            T1& X,
            uword dim
            ) {
    uword size = 0;
    if (dim == 0) {
        size = X.n_rows;
    } else if(dim == 1) {
        size = X.n_cols;
    } else if(dim > 1) {
        arma_extra_debug_sigprint();
        arma_debug_print( "fftshift(): trying to shift dimension greater than 2");
    }

    T1 out = circshift(X,dim,std::floor(size/2));
    return out;
}

template<typename T1>
arma_inline
T1 fftshift(
            T1 X
            ) {

    T1 out = circshift(circshift(X,0,std::floor(X.n_rows/2)),1,std::floor(X.n_cols/2));
    return out;
}



#endif /* defined(__PowerGrid__fftshift__hpp) */
