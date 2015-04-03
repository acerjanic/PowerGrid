//
//  Gfft.hpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 3/25/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//
#ifndef __PowerGrid__Gfft__hpp
#define __PowerGrid__Gfft__hpp

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
    
    uword n1 = 0; //Data size (n1*n1)
    uword n2 = 0; //Image size (n2,n2)

    //Overloaded methods for forward and adjoint transform
    //Forward transform operation
    T1 operator*(const T1& d)
    {
        uword stx = (this->n2-this->n1)/2;
        //Create 2D array object for use with the fft
        Mat<cx_double> d2D = reshape(d,this->n2,this->n2);
        
        d2D = fftshift(fft2(fftshift(d2D)));
        Mat<cx_double> d2Dtrimmed = d2D(span(stx,stx+this->n1-1),span(stx,stx+this->n1-1));
        //equivalent to returning col(output) in MATLAB with IRT
        return vectorise(d2Dtrimmed);
        
    }
    
    //Adjoint transform operation
    T1 operator/(const T1& d)
    {
        uword stx = (this->n2-this->n1)/2;

        //Create 2D array object for use with the fft
        Mat<cx_double> d2D = reshape(d,this->n1,this->n1);
        Mat<cx_double> d2Dtrimmed = d2D(span(stx,stx+this->n1-1),span(stx,stx+this->n1-1));

        d2Dtrimmed = this->n2*this->n2*ifft2(d2Dtrimmed);
        //equivalent to returning col(output) in MATLAB with IRT
        return vectorise(d2Dtrimmed);
        
    }

};
#endif