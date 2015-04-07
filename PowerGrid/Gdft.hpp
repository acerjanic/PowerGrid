//
//  Gdft.hpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 4/4/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef PowerGrid_Gdft_hpp
#define PowerGrid_Gdft_hpp
template<typename T1,typename T2> //This is of type complex<double> or complex<float>, or any other type like float or single
class Gdft {
    
public:
    //Default Class Constructor and Destructor
    Gdft();
    //Class Constructor
    Gdft(uword a, uword b, const Col<T2> &k1, const Col<T2> &k2, const Col<T2> &k3) //Change these argumenst as you need to setup the object
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
    
    Col<T2> kx;
    Col<T2> ky;
    Col<T2> kz;
    
    //Overloaded methods for forward and adjoint transform
    //Forward transform operation
    Col<T1> operator*(const Col<T1>& d) const//Don't change these arguments
    {
        //This is just specifying size assuming things are the same size, change as necessary
        uword dataLength = d.n_rows;
        Col<T2> realData = real(d);
        Col<T2> imagData = imag(d);
        //Now we grab the data out of armadillo with the memptr() function
        //This returns a pointer of the type of the elements of the array/vector/matrix/cube (3d matrix)
        //Armadillo uses column major like MATLAB and Fortran, but different from 2D C++ arrays which are row major.
        T2* realDataPtr = realData.memptr();
        T2* imagDataPtr = imagData.memptr();
        
        Col<T2> realXformedData;
        Col<T2> imagXformedData;
        realXformedData.copy_size(realData);
        imagXformedData.copy_size(realData);
        realXformedData.zeros();
        imagXformedData.zeros();
        
        Col<T2> kx_r = real(this->kx);
        Col<T2> ky_r = real(this->ky);
        Col<T2> kz_r = real(this->kz);
        Col<T2> kx_i = imag(this->kx);
        Col<T2> ky_i = imag(this->ky);
        Col<T2> kz_i = imag(this->kz);
        
        T2* realXformedDataPtr = realXformedData.memptr();
        T2* imagXformedDataPtr = imagXformedData.memptr();
        //Process data here, like calling a brute force transform, dft...
        // I assume you create the pointers to the arrays where the transformed data will be stored
        // realXformedDataPtr and imagXformedDataPtr and they are of type float*
        ftCpu<T2>(realXformedDataPtr,imagXformedDataPtr,
                  realDataPtr, imagDataPtr, kx_r.memptr(),
                  ky_r.memptr(), kz_r.memptr(),
                  kx_i.memptr(), ky_i.memptr(), kz_i.memptr(),
                  this->n1, this->n2
                  );
        
        //To return data, we need to put our data back into Armadillo objects
        //We are telling the object how long it is because it will copy the data back into managed memory
        //realXformedData(realXformedDataPtr, dataLength);
        //imagXformedData(imagXformedDataPtr, dataLength);
        
        //We can free the realDataXformPtr and imagDataXformPtr at this point and Armadillo will manage armadillo object memory as things change size or go out of scope and need to be destroyed
        
        Col<T1> XformedData(dataLength);
        XformedData.set_real(realXformedData);
        XformedData.set_imag(imagXformedData);
        
        return conv_to<Col<T1>>::from(XformedData); //Return a vector of type T1
        
    }
    
    //Adjoint transform operation
    Col<T1> operator/(const Col<T1>& d) const
    {
        
        uword dataLength = size(d,1);
        
        Col<T2> realData = real(d);
        Col<T2> imagData = imag(d);
        
        T2* realDataPtr = realData.memptr();
        T2* imagDataPtr = imagData.memptr();
        
        Col<T2> realXformedData;
        Col<T2> imagXformedData;
        realXformedData.copy_size(realData);
        imagXformedData.copy_size(realData);
        realXformedData.zeros();
        imagXformedData.zeros();
        
        Col<T2> kx_r = real(this->kx);
        Col<T2> ky_r = real(this->ky);
        Col<T2> kz_r = real(this->kz);
        Col<T2> kx_i = imag(this->kx);
        Col<T2> ky_i = imag(this->ky);
        Col<T2> kz_i = imag(this->kz);
        
        T2* realXformedDataPtr = realXformedData.memptr();
        T2* imagXformedDataPtr = imagXformedData.memptr();
        //Process data here, like calling a brute force transform, dft...
        // I assume you create the pointers to the arrays where the transformed data will be stored
        // realXformedDataPtr and imagXformedDataPtr and they are of type float*
        iftCpu<T2>(realXformedDataPtr,imagXformedDataPtr,
                  realDataPtr, imagDataPtr, kx_r.memptr(),
                  ky_r.memptr(), kz_r.memptr(),
                  kx_i.memptr(), ky_i.memptr(), kz_i.memptr(),
                  this->n1, this->n2
                  );
        
        //realXformedData(realXformedDataPtr, dataLength);
        //imagXformedData(imagXformedDataPtr, dataLength);
        
        //We can free the realDataXformPtr and imagDataXformPtr at this point and Armadillo will manage armadillo object memory as things change size or go out of scope and need to be destroyed
        
        Col<T1> XformedData(dataLength);
        XformedData.set_real(realXformedData);
        XformedData.set_imag(imagXformedData);
        
        return conv_to<Col<T1>>::from(XformedData); //Return a vector of type T1
        
    }
    
};

#endif
