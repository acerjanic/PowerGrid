#include "PowerGrid.h"

using namespace arma;

// Define code for testing
Mat<cx_double> test_pwls_pcg()
{

   // Load real data
   Mat<double> test;

   //Load our read double matrix object. (You need to match type to avoide mangling the data.)
   loadmat("C:\\Ceemple\\user\\PowerGrid\\Resources\\test.mat","test",&test);

   // Formard operator
   Gfft<cx_double> G(64,64);
   Col<cx_double> TestForwars;
   TestForward = G *vectorise(conv_to<Mat<cx_double>>::from(test));

   // Variables needed for the recon: Penalty object, num of iterations 
   umat ReconMask;
   ReconMask.ones(64,64);

   QuadPenalty<cx_double>R(64,64,1,ReconMask);

   uword niter = 2;
   Col<cx_double> xinit = zeros<Col<cx_double>>(64*64); // initial estimate of x
   Mat<cx_double> W;
   W = eye<Mat<cx_double>>(G.n1,G.n1); // Should be the size of k-space data: Is it right?

   Col<cx_double> x_t;
   x_t = pwls_pcg1<cx_double,Gfft<cx_double>,QuadPenalty<cx_double>>(xinit, G, W, TestForward, R, niter);
	
   return x_t;
}

