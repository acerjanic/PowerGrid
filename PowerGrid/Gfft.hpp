//
//  Gfft.hpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 3/25/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//
#ifndef __PowerGrid__Gfft__hpp
#define __PowerGrid__Gfft__hpp
//We want to operate on many types of variables (assume of type Col<type>)
template<typename T1>
class Gfft {
    
public:
    //Default Class Constructor and Destructor
    Gfft();
    //Class Constructor
    Gfft(uword a, uword b)
    {
        n1 = a;
        n2 = b;
    }
    
    //Class variables go here.
    uword n1 = 0; //Data size (n1*n1)
    uword n2 = 0; //Image size (n2,n2)

    //Overloaded methods for forward and adjoint transform
    //Forward transform operation
    Col<T1> operator*(const Col<T1>& d) const
    {
        uword stx = (this->n2-this->n1)/2;

        //Create 2D array object for use with the fft
        Mat<T1> d2D = reshape(d,this->n2,this->n2);
        
        d2D = fftshift(fft2(fftshift(d2D)));
        Mat<T1> d2Dtrimmed = d2D(span(stx,stx+this->n1-1),span(stx,stx+this->n1-1));
        //equivalent to returning col(output) in MATLAB with IRT
        return vectorise(d2Dtrimmed);
        
    }
    
    //Adjoint transform operation
    Col<T1> operator/(const Col<T1>& d) const
    {
        uword stx = (this->n2-this->n1)/2;

        //Create 2D array object for use with the fft
        Mat<T1> d2D = reshape(d,this->n1,this->n1);
        Mat<T1> d2Dtrimmed = zeros<Mat<T1>>(this->n2,this->n2);

        d2Dtrimmed = d2D(span(stx,stx+this->n1-1),span(stx,stx+this->n1-1));

        //d2Dtrimmed = this->n2*this->n2*fftshift(fft2(fftshift(d2Dtrimmed)));
        //NOTE: Armadillo's ffts appear to be normalized with the 1/N^2 in the 2D case already build into ifft2.
        //We don't need the n2*n2 term anymore, we can match the MATLAB Gfft behavior without it.
        
        d2Dtrimmed = fftshift(ifft2(fftshift(d2Dtrimmed)));
        //equivalent to returning col(output) in MATLAB with IRT
        return vectorise(d2Dtrimmed);
        
    }

};
#endif
