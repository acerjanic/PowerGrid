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
    string testPath = "/Users/alexcerjanic/Developer/PG/Resources/";

    //Mat<cx_double> testComplex;
    //Load our read double matrix object. (You need to match type to avoide mangling the data.)
    loadmat(testPath+"test.mat","test",&test);
    
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

    savemat("/shared/mrfil-data/acerja2/repos/PowerGrid/Resources/out.mat","fftShiftOut", fftshift(test).eval());
    
    //Writing the transforms to disk for evaluation in MATLAB. Note that we can't read or write complex natively (yet), so lets
    //seperate real and imaginary first. Also, we can't put multiple variables in a single file yet. These are TODOs.
    savemat(testPath+"testForwardReal.mat","testForwardReal", real(testForward).eval());
    savemat(testPath+"testForwardImag.mat","testForwardImag", imag(testForward).eval());
    savemat(testPath+"testForwardTest.mat","testForward", testForward);
    savemat(testPath+"testAdjointTest.mat","testAdjoint", testAdjoint);

    savemat(testPath+"testAdjointReal.mat","testAdjointReal", real(testAdjoint).eval());
    savemat(testPath+"testAdjointImag.mat","testAdjointImag", imag(testAdjoint).eval());
    
    //Test SENSE

    //Create an empty unsized matrix of type real double
    //Mat<double> SMap;

    //Load our read double matrix object. (You need to match type to avoide mangling the data.)
    //loadmat("/Users/alexcerjanic/Developer/PG/Resources/SMap.mat","SMap",&SMap);

    //create our SENSE object, coils =2
    //SENSE<cx_double, Gfft<cx_double>> S(G,vectorise(conv_to<Col<cx_double>>::from(SMap)),4096*2,4096,2);

    //Perform the forward transformation and store the result in a colum vector of type complex double
     //Note the conv_to<type to convert to>::from(data) command to convert our real double phantom to type complex double
     //Col<cx_double> testForward_SENSE = S*vectorise(conv_to<Mat<cx_double>>::from(test));

     //Now perform the adjoint transform
     //Col<cx_double> testAdjoint_SENSE = S/testForward_SENSE;

     //savemat("/Users/alexcerjanic/Developer/PG/Resources/testForwardTest_SENSE.mat","testForward_SENSE", testForward_SENSE);
     //savemat("/Users/alexcerjanic/Developer/PG/Resources/testAdjointTest_SENSE.mat","testAdjoint_SENSE", testAdjoint_SENSE);
     //Test PWLS_PCG1
    
     //Mat<cx_double> testPWLS = test_pwls_pcg();
     //savemat("/Users/alexcerjanic/Developer/PG/Resources/testPWLS.mat","testPWLS", testPWLS);
    
    Mat<cx_double> cxData;
    Col<double> kx;
    Col<double> ky;
    Col<double> kz;
    loadmat(testPath+"kx_sp.mat","kx",&kx);
    loadmat(testPath+"ky_sp.mat","ky",&ky);
    //loadmat("/shared/mrfil-data/jholtrop/repos/PowerGrid/Resources/data_gdft.mat","data_gdft",&cxData);

    //loadmat("/Users/alexcerjanic/Developer/PG/Resources/testGdftTraj.mat","kz",&kz);
    kz.copy_size(ky);
    kz.zeros();
    //Col<cx_double> out = test_gdft<cx_double,double>(conv_to<Col<cx_double>>::from(test),kx,ky,kz);
    //savemat(testPath+"testGdft.mat","testGdft",out);

   Col<cx_double> out_FieldCorrected = test_FieldCorrection<cx_double,double>();
   savemat(testPath+"test_FieldCorrection.mat","test_FieldCorrection",out_FieldCorrected);

    return 0;
}

