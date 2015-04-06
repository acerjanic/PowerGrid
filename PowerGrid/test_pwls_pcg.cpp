#include "PowerGrid.h"

using namespace arma;

// Define code for testing
Mat<cx_double> test_pwls_pcg()
{
   // Simulate the k-space data
    Mat<cx_double> x = ones<Mat<cx_double>>(64,64);

   Gfft<Col<cx_double>> G(64,64);

   Col<cx_double> y;
   y = G *vectorise(x);

   // Variables needed for the recon: Penalty object, num of iterations 
   umat ReconMask;
   ReconMask.ones(64,64);
   QuadPenalty<Col<cx_double>>R(64,64,1,ReconMask);
   uword niter = 20;
   Col<cx_double> xinit = zeros<Col<cx_double>>(64*64); // initial estimate of x
   Col<cx_double> W;
   W << 1;
    Col<cx_double> x_t;
    x_t = pwls_pcg1<Col<cx_double>,Gfft<Col<cx_double>>,QuadPenalty<Col<cx_double>>>(xinit, G, W, y, R, niter);
	
   return x_t;
}

