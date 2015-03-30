//
//  fftshift.hpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 3/12/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef PowerGrid_fftshift_hpp
#define PowerGrid_fftshift_hpp


//  Create an armadillo compatible circshift function similar to MATLAB's

template<typename T1>
inline const Gen<T1, op_circshift>
circshift(const Base<typename T1::elem_type,T1>& X,
const uword n_dim, const uword shift)
{
    arma_extra_debug_sigprint();
    arma_debug_check ( (n_dim > 3), "circshfit(): trying to shift dimension greater than 3");
    
    
    return Gen<

}



#endif
