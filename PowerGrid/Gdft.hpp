//
//  Gdft.hpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 4/4/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef PowerGrid_Gdft_hpp
#define PowerGrid_Gdft_hpp
template<typename T1> //This is of type complex<double> or complex<float>, or any other type like float or single
class Gdft {
    
public:
    //Default Class Constructor and Destructor
    Gdft();
    //Class Constructor
    Gdft(uword a, uword b, Col<complex<T1>> &k1, Col<complex<T1>> &k2, Col<complex<T1>> &k3) //Change these argumenst as you need to setup the object
    {
        n1 = a;
        n2 = b;
        kx = k1;
        ky = k2;
        kz = k3;
    }
    
    //Class variables go here. Change as necessary
    uword n1 = 0;
    uword n2 = 0;
    
    Col<T1> kx;
    Col<T1> ky;
    Col<T1> kz;
    
    //Overloaded methods for forward and adjoint transform
    //Forward transform operation
    Col<complex<T1>> operator*(const Col<complex<T1>>& d) const//Don't change these arguments
    {
        //This is just specifying size assuming things are the same size, change as necessary
        uword dataLength = size(d,1);
        Col<T1> realData = real(d);
        Col<T1> imagData = imag(d);
        //Now we grab the data out of armadillo with the memptr() function
        //This returns a pointer of the type of the elements of the array/vector/matrix/cube (3d matrix)
        //Armadillo uses column major like MATLAB and Fortran, but different from 2D C++ arrays which are row major.
        T1  *realDataPtr = realData.memptr();
        T1  *imagDataPtr = imagData.memptr();
        
        Col<T1> realXformedData;
        Col<T1> imagXformedData;
        realXformedData.copy_size(realData);
        imagXformedData.copy_size(realData);
        realXformedData.zeros();
        imagXformedData.zeros();
        
        Col<T1> kx_r = real(this->kx);
        Col<T1> ky_r = real(this->ky);
        Col<T1> kz_r = real(this->kz);
        Col<T1> kx_i = imag(this->kx);
        Col<T1> ky_r = imag(this->ky);
        Col<T1> kz_r = imag(this->kz);
        
        T1  *realXformedDataPtr = realXformedData.memptr();
        T1  *imagXformedDataPtr = realXformedData.memptr();
        //Process data here, like calling a brute force transform, dft...
        // I assume you create the pointers to the arrays where the transformed data will be stored
        // realXformedDataPtr and imagXformedDataPtr and they are of type float*
        ftCpu<T1>(realXformedDataPtr,imagXformedDataPtr,
                  realDataPtr, imagDataPtr, kx_r.memptr(),
                  ky_r.memptr(), kz_r.memptr(),
                  kx_i.memptr(), ky_i.memptr(), kz_i.memptr(),
                  this->n1, this->n2
                  );
        
        //To return data, we need to put our data back into Armadillo objects
        //We are telling the object how long it is because it will copy the data back into managed memory
        Col<T1> realXformedData(realDataXformPtr, dataLength);
        Col<T1> imagXformedData(imagDataXformPtr, dataLength);
        
        //We can free the realDataXformPtr and imagDataXformPtr at this point and Armadillo will manage armadillo object memory as things change size or go out of scope and need to be destroyed
        
        Col<cx_float> XformedData(dataLength);
        XformedData.set_real(realXformedData);
        XformedData.set_imag(imagXformedData);
        
        return conv_to<T1>::from(XformedData); //Return a vector of type T1
        
    }
    
    //Adjoint transform operation
    Col<T1> operator/(const Col<T1>& d) const
    {
        
        uword dataLength = size(d,1);
        
        Col<T1> realData = real(d);
        Col<T1> imagData = imag(d);
        
        float *realDataPtr = realData.memptr();
        float *imagDataPtr = imagData.memptr();
        
        Col<T1> realXformedData;
        Col<T1> imagXformedData;
        realXformedData.copy_size(realData);
        imagXformedData.copy_size(realData);
        realXformedData.zeros();
        imagXformedData.zeros();
        
        Col<T1> kx_r = real(this->kx);
        Col<T1> ky_r = real(this->ky);
        Col<T1> kz_r = real(this->kz);
        Col<T1> kx_i = imag(this->kx);
        Col<T1> ky_r = imag(this->ky);
        Col<T1> kz_r = imag(this->kz);
        
        T1  *realXformedDataPtr = realXformedData.memptr();
        T1  *imagXformedDataPtr = realXformedData.memptr();
        //Process data here, like calling a brute force transform, dft...
        // I assume you create the pointers to the arrays where the transformed data will be stored
        // realXformedDataPtr and imagXformedDataPtr and they are of type float*
        iftCpu<T1>(realXformedDataPtr,imagXformedDataPtr,
                  realDataPtr, imagDataPtr, kx_r.memptr(),
                  ky_r.memptr(), kz_r.memptr(),
                  kx_i.memptr(), ky_i.memptr(), kz_i.memptr(),
                  this->n1, this->n2
                  );
        
        Col<float> realXformedData(realDataXformPtr, dataLength);
        Col<float> imagXformedData(imagDataXformPtr, dataLength);
        
        //We can free the realDataXformPtr and imagDataXformPtr at this point and Armadillo will manage armadillo object memory as things change size or go out of scope and need to be destroyed
        
        Col<cx_float> XformedData(dataLength);
        XformedData.set_real(realXformedData);
        XformedData.set_imag(imagXformedData);
        
        return conv_to<T1>::from(XformedData); //Return a vector of type T1
        
    }
    
};

#endif
