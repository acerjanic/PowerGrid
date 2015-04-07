//
//  main.cpp
//  PowerGrid
//
//  Created by Alex Cerjanic on 3/12/15.
//  Copyright (c) 2015 MRFIL. All rights reserved.
//

#include "PowerGrid.h" //Project headers.
#include "../Support/CeempleMatio.h" //Headers for using savemat and loadmat

using namespace arma; //Armdillo stuff is in the arma namespace
using namespace std; //complex type comes from the STL


int main(int argc, char** argv)
{
    //Create an empty unsized matrix of type real double
    Mat<double> test;
    
    //Mat<cx_double> testComplex;
    //Load our read double matrix object. (You need to match type to avoide mangling the data.)
    loadmat("/Users/alexcerjanic/Developer/PG/Resources/test.mat","test",&test);
    
    //loadmat("/Users/alexcerjanic/Developer/PG/Resources/testForwardTest.mat","testForward",&testComplex);
    //Create our Gfft object for type complex double and N1=64, N2=64
    Gfft<cx_double> G(64,64);

    //Perform the forward transformation and store the result in a colum vector of type complex double
    //Note the conv_to<type to convert to>::from(data) command to convert our real double phantom to type complex double
    Col<cx_double> testForward = G*vectorise(conv_to<Mat<cx_double>>::from(test));
    
    //Now perform the adjoint transform
    Col<cx_double> testAdjoint = G/testForward;
    
    //Testing the my fftshift function. Armadillo delays evaluations to batch them. You can force the evaluation with data.eval()
    //If we just passed the variable rather than the variable wrapped in a function, we wouldn't need the eval.
    //savemat("/Users/alexcerjanic/Developer/PG/Resources/testForwardOut.mat","testForwardOut",test);

    savemat("/Users/alexcerjanic/Developer/PG/Resources/out.mat","fftShiftOut", fftshift(test).eval());
    
    //Writing the transforms to disk for evaluation in MATLAB. Note that we can't read or write complex natively (yet), so lets
    //seperate real and imaginary first. Also, we can't put multiple variables in a single file yet. These are TODOs.
    savemat("/Users/alexcerjanic/Developer/PG/Resources/testForwardReal.mat","testForwardReal", real(testForward).eval());
    savemat("/Users/alexcerjanic/Developer/PG/Resources/testForwardImag.mat","testForwardImag", imag(testForward).eval());
    savemat("/Users/alexcerjanic/Developer/PG/Resources/testForwardTest.mat","testForward", testForward);
    savemat("/Users/alexcerjanic/Developer/PG/Resources/testAdjointTest.mat","testAdjoint", testAdjoint);

    savemat("/Users/alexcerjanic/Developer/PG/Resources/testAdjointReal.mat","testAdjointReal", real(testAdjoint).eval());
    savemat("/Users/alexcerjanic/Developer/PG/Resources/testAdjointImag.mat","testAdjointImag", imag(testAdjoint).eval());
    
    //Test SENSE

    //Create an empty unsized matrix of type real double
    Mat<cx_double> SMap;

    //Load our read double matrix object. (You need to match type to avoide mangling the data.)
    loadmat("/Users/alexcerjanic/Developer/PG/Resources/SMap.mat","SMap",&SMap);

    //create our SENSE object, coils =2
    SENSE<cx_double, Gfft<cx_double>> S(G,vectorise(SMap),4096*2,4096,2);

    //Perform the forward transformation and store the result in a colum vector of type complex double
     //Note the conv_to<type to convert to>::from(data) command to convert our real double phantom to type complex double
     Col<cx_double> testForward_SENSE = S*vectorise(conv_to<Mat<cx_double>>::from(test));

     //Now perform the adjoint transform
     Col<cx_double> testAdjoint_SENSE = S/testForward_SENSE;

     savemat("/Users/alexcerjanic/Developer/PG/Resources/testForwardTest_SENSE.mat","testForward_SENSE", testForward_SENSE);
     savemat("/Users/alexcerjanic/Developer/PG/Resources/testAdjointTest_SENSE.mat","testAdjoint_SENSE", testAdjoint_SENSE);
     //Test PWLS_PCG1
    
     Mat<cx_double> testPWLS = test_pwls_pcg();
     savemat("/Users/alexcerjanic/Developer/PG/Resources/testPWLS.mat","testPWLS", testPWLS);

    return 0;
}

