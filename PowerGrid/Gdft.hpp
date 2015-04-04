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
    Gdft(uword a, uword b) //Change these argumenst as you need to setup the object
    {
        n1 = a;
        n2 = b;
    }
    
    //Class variables go here. Change as necessary
    uword n1 = 0;
    uword n2 = 0;
    
    //Overloaded methods for forward and adjoint transform
    //Forward transform operation
    Col<T1> operator*(const Col<T1>& d) //Don't change these arguments
    {
        //This is just specifying size assuming things are the same size, change as necessary
        uword dataLength = size(d,1);
        Col<float> realData = real(d);
        Col<float> imagData = imag(d);
        
        //Now we grab the data out of armadillo with the memptr() function
        //This returns a pointer of the type of the elements of the array/vector/matrix/cube (3d matrix)
        //Armadillo uses column major like MATLAB and Fortran, but different from 2D C++ arrays which are row major.
        float  *realDataPtr = realData.memptr();
        float  *imagDataPtr = imagData.memptr();
        
        //Process data here, like calling a brute force transform, dft...
        // I assume you create the pointers to the arrays where the transformed data will be stored
        // realXformedDataPtr and imagXformedDataPtr and they are of type float*
        
        
        //To return data, we need to put our data back into Armadillo objects
        //We are telling the object how long it is because it will copy the data back into managed memory
        Col<float> realXformedData(realDataXformPtr, dataLength);
        Col<float> imagXformedData(imagDataXformPtr, dataLength);
        
        //We can free the realDataXformPtr and imagDataXformPtr at this point and Armadillo will manage armadillo object memory as things change size or go out of scope and need to be destroyed
        
        Col<cx_float> XformedData(dataLength);
        XformedData.set_real(realXformedData);
        XformedData.set_imag(imagXformedData);
        
        return conv_to<T1>::from(XformedData); //Return a vector of type T1
        
    }
    
    //Adjoint transform operation
    T1 operator/(const T1& d)
    {
        
        uword dataLength = size(d,1);

        Col<float> realData = real(d);
        Col<float> imagData = imag(d);
        
        float *realDataPtr = realData.memptr();
        float *imagDataPtr = imagData.memptr();
        
        //Process data here
        // I assume you create the pointers to the arrays where the transformed data will be stored
        // realXformedDataPtr and imagXformedDataPtr and they are of type float*
        
        //To return data, we need to put our data back into Armadillo objects
        //We are telling the object how long it is because it will copy the data back into managed memory
        //Later on we can talk about how to keep the data in already allocated memory, but it is more complicated and error prone.
        
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
