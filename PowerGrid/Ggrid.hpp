//
//  Ggrid.hpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 4/8/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#ifndef PowerGrid_Ggrid_h
#define PowerGrid_Ggrid_h
template<typename T1,typename T2> //This is of type complex<double> or complex<float>, or any other type like float or single
class Ggrid {
    
public:
    
    //Default Class Constructor and Destructor
    Ggrid();

    //~Ggrid();
    //Class Constructor
    Ggrid(uword dataLength, T2 gridos, uword nx, uword ny, uword nz, const Col <T2> &k1, const Col <T2> &k2,
          const Col <T2> &k3, const Col <T2> &i1, const Col <T2> &i2,
          const Col <T2> &i3) //Change these argumenst as you need to setup the object
    {
        cout << "Entering the class constructor for Ggrid" << endl;
        n1 = nx*ny*nz;
        n2 = dataLength;
        Nx = nx;
        Ny = ny;
        Nz = nz;
        ix = i1;
        iy = i2;
        iz = i3;
        kx = k1;
        ky = k2;
        kz = k3;
        gridOS = gridos;
        //Set Beta
        kernelWidth = 4.0;
        beta = MRI_PI * std::sqrt((gridOS - 0.5) * (gridOS - 0.5) *
                                  (kernelWidth * kernelWidth * 4.0) /
                                  (gridOS * gridOS) - 0.8);

        //Deal with the LUT
        //Generating Look-Up Table
        cout << "Calculating look up table" << endl;
        calculateLUT(beta, kernelWidth, LUT, sizeLUT);
        #pragma acc enter data copyin(LUT[0:sizeLUT])

    }

    //Class destructor to free LUT
    ~Ggrid() {
        if (LUT) {
            #pragma acc exit data delete(LUT)
            free(LUT);
        }
    }
    
    //Class variables go here. Change as necessary
    uword n1 = 0;
    uword n2 = 0;
    uword Nx = 0;
    uword Ny = 0;
    uword Nz = 0;
    
    Col<T2> kx; //k-space coordinates
    Col<T2> ky;
    Col<T2> kz;
    Col<T2> ix; //image space coordinates
    Col<T2> iy;
    Col<T2> iz;

    T2 gridOS; //grid oversampling
    T2 *LUT = 0; // Lookup table for the gridding operations
    uword sizeLUT = 0;
    T2 beta; //beta factor for gridding not the same as beta in regularization!
    T2 kernelWidth; //Kaiser Bessel Kernel Support
    
    //Overloaded methods for forward and adjoint transform
    //Forward transform operation using gridding
    Col<T1> operator*(const Col<T1>& d) const//Don't change these arguments
    {
        cout << "Entering forward operator overload in Ggrid." << endl;
        //This is just specifying size assuming things are the same size, change as necessary
        //uword dataLength = d.n_rows;
        /* TODO: Fix this. It eats up more than 32GB of memory!
        Col<T2> FM(ix.n_rows*iy.n_rows);
        Col<T2> t(n2);
        FM.zeros();
        t.zeros();
        */
        Col<T2> realData = real(d);
        Col<T2> imagData = imag(d);
        //Now we grab the data out of armadillo with the memptr() function
        //This returns a pointer of the type of the elements of the array/vector/matrix/cube (3d matrix)
        //Armadillo uses column major like MATLAB and Fortran, but different from 2D C++ arrays which are row major.
        T2* realDataPtr = realData.memptr();
        T2* imagDataPtr = imagData.memptr();
        
        Col<T2> realXformedData(this->n2);
        Col<T2> imagXformedData(this->n2);
        //realXformedData.zeros();
        //imagXformedData.zeros();
        
        T2* realXformedDataPtr = realXformedData.memptr();
        T2* imagXformedDataPtr = imagXformedData.memptr();
        //Process data here, like calling a brute force transform, dft...
        // I assume you create the pointers to the arrays where the transformed data will be stored
        // realXformedDataPtr and imagXformedDataPtr and they are of type float*
        /*
        ftCpu<T2>(realXformedDataPtr,imagXformedDataPtr,
                  realDataPtr, imagDataPtr, kx.memptr(),
                  ky.memptr(), kz.memptr(),
                  ix.memptr(), iy.memptr(), iz.memptr(),
                  FM.memptr(), t.memptr(),
                  this->n2, this->n1
        );
        */
        //T2 gridOS = 2.0;
        computeFd_CPU_Grid<T2>(n2, kx.memptr(),  ky.memptr(),  kz.memptr(),
                               realDataPtr, imagDataPtr, Nx, Ny, Nz,
                               gridOS, realXformedDataPtr, imagXformedDataPtr, kernelWidth, beta, LUT, sizeLUT);


        //To return data, we need to put our data back into Armadillo objects
        //We are telling the object how long it is because it will copy the data back into managed memory
        //realXformedData(realXformedDataPtr, dataLength);
        //imagXformedData(imagXformedDataPtr, dataLength);
        
        //We can free the realDataXformPtr and imagDataXformPtr at this point and Armadillo will manage armadillo object memory as things change size or go out of scope and need to be destroyed
        
        Col<T1> XformedData(this->n2);
        XformedData.set_real(realXformedData);
        XformedData.set_imag(imagXformedData);

        return conv_to<Col<T1>>::from(XformedData); //Return a vector of type T1
        
    }
    
    //Adjoint transform operation
    Col<T1> operator/(const Col<T1>& d) const
    {
        
        uword dataLength = n2;
        /* TODO: Fix this. It eats up more than 32GB of memory!
        Col<T2> FM(ix.n_rows*iy.n_rows);
        Col<T2> t(n2);
        FM.zeros();
        t.zeros();
        */
        Col<T2> realData = real(d);
        Col<T2> imagData = imag(d);
        
        T2* realDataPtr = realData.memptr();
        T2* imagDataPtr = imagData.memptr();
        
        Col<T2> realXformedData(n1);
        Col<T2> imagXformedData(n1);

        //realXformedData.zeros();
        //imagXformedData.zeros();
        
        T2* realXformedDataPtr = realXformedData.memptr();
        T2* imagXformedDataPtr = imagXformedData.memptr();
        //Process data here, like calling a brute force transform, dft...
        // I assume you create the pointers to the arrays where the transformed data will be stored
        // realXformedDataPtr and imagXformedDataPtr and they are of type float*

        //T2 gridOS = 2.0;

        computeFH_CPU_Grid<T2>(dataLength, kx.memptr(),  ky.memptr(),  kz.memptr(),
                           realDataPtr, imagDataPtr, Nx, Ny, Nz,
                               gridOS, realXformedDataPtr, imagXformedDataPtr, kernelWidth, beta, LUT, sizeLUT);
        /*
        iftCpu<T2>(realXformedDataPtr,imagXformedDataPtr,
                   realDataPtr, imagDataPtr, kx.memptr(),
                   ky.memptr(), kz.memptr(),
                   ix.memptr(), iy.memptr(), iz.memptr(),
                   FM.memptr(), t.memptr(),
                   this->n2, this->n1
                   );
        */
        //realXformedData(realXformedDataPtr, dataLength);
        //imagXformedData(imagXformedDataPtr, dataLength);
        
        //We can free the realDataXformPtr and imagDataXformPtr at this point and Armadillo will manage armadillo object memory as things change size or go out of scope and need to be destroyed
        
        Col<T1> XformedData(n1);
        XformedData.set_real(realXformedData);
        XformedData.set_imag(imagXformedData);
        //savemat("/shared/mrfil-data/data/PowerGridTest/64_64_16_4coils/ggrid.mat","img",XformedData);

        return conv_to<Col<T1>>::from(XformedData); //Return a vector of type T1
        
    }
    
};

#endif
