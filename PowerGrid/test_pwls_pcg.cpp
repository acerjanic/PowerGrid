#include "CeempleArmadillo.h"
#include <iostream>

#include "PowerGrid.h" //Project headers.

// Define code for testing
int main()
{
   // Simulate the k-space data
   Mat<cx_double> x = ones(64,64);

   Gfft<Col<cx_double>> G(64,64);

   col<cx_double> y;
   y = G *vectorise(x);

   // Variables needed for the recon: Penalty object, num of iterations 
   
   QuadPenalty<col<cx_double>>R;
   unword niter = 20;	
   col<cx_double> xinit = zeros(64*64,1); // initial estimate of x
   col<cx_double> W = 1;
   x_t = pwls_pcg1(xinit, G, W, y, R, niter);
	
   return 0;
}

