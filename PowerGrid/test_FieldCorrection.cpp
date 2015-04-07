#include "PowerGrid.h"

using namespace arma;

// Define code for testing
void test_fieldCorrection()
{
  //Create an empty unsized matrix of type real double
  Mat<double> test;

  //Mat<cx_double> testComplex;
  //Load our read double matrix object. (You need to match type to avoide mangling the data.)
  loadmat("C:\\Ceemple\\user\\PowerGrid\\Resources\\test.mat","test",&test);

  //loadmat("/Users/alexcerjanic/Developer/PowerGrid/Resources/testForwardTest.mat","testForward",&testComplex);
  //Create our Gfft object for type complex double and N1=64, N2=64
  Gfft<Col<cx_double>> G(64,64);

  //Test SENSE

  //Create an empty unsized matrix of type real double
  Col<cx_double> SMap;

  //Load our read double matrix object. (You need to match type to avoide mangling the data.)
  loadmat("C:\\Ceemple\\user\\PowerGrid\\Resources\\SMap.mat","SMap",&SMap);

  //create our SENSE object, coils =2
  SENSE<cx_double, Gfft<cx_double>> S(G,SMap,4096*2,4096,2);

  //Perform the forward transformation and store the result in a colum vector of type complex double
   //Note the conv_to<type to convert to>::from(data) command to convert our real double phantom to type complex double
   Col<cx_double> testForward_SENSE = S*vectorise(conv_to<Mat<cx_double>>::from(test));

   //Now perform the adjoint transform
   Col<cx_double> testAdjoint_SENSE = S/testForward_SENSE;

   savemat("C:\\Ceemple\\user\\PowerGrid\\Resources\\testForwardTest_SENSE.mat","testForward_SENSE", testForward_SENSE);
   savemat("C:\\Ceemple\\user\\PowerGrid\\Resources\\testAdjointTest_SENSE.mat","testAdjoint_SENSE", testAdjoint_SENSE);

   //test time segmentation

   Col<cx_double> FMap;
   Col<double> t_vec;

   loadmat("C:\\Ceemple\\user\\PowerGrid\\Resources\\FMap.mat","FMap",&FMap);
   loadmat("C:\\Ceemple\\user\\PowerGrid\\Resources\\t_vec.mat","t_vec",&t_vec);

   uword L = 20;

   FieldCorrection<cx_double, double, Gfft<cx_double>> A(G,FMap,t_vec,4096,4096,L);

   //Perform the forward transformation and store the result in a colum vector of type complex double
    //Note the conv_to<type to convert to>::from(data) command to convert our real double phantom to type complex double
    Col<cx_double> testForward_FC = A*vectorise(conv_to<Mat<cx_double>>::from(test));

    //Now perform the adjoint transform
    Col<cx_double> testAdjoint_FC = A/testForward_FC;

    savemat("C:\\Ceemple\\user\\PowerGrid\\Resources\\testForwardTest_FC.mat","testForward_FC", testForward_FC);
    savemat("C:\\Ceemple\\user\\PowerGrid\\Resources\\testAdjointTest_FC.mat","testAdjoint_FC", testAdjoint_FC);



}

