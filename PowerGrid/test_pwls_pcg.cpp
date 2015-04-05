#include "CeempleArmadillo.h"
#include <iostream>

// question 
// no preconditinned gradient
// no stepper choice
// function Robject
// complexe reall????
// .t needs to be transpose conjugate

// Function prototype
void pwls_pcg(mat x, mat A, mat W, vec yi, Robj R, double niter);

// Define code for testing
int main()
{
   return 0;
}

// Function 

void pwls_pcg(cx_mat& x, fast_mr const& A, mat const& W, cx_mat const& yi, Robj const& R, double niter)
{

 // Initialize projection

  cx_mat Ax = A*x;

  for (unsigned int ii(0); ii < niter; +ii)
  {
    // Compute negative gradient
    cx_mat ngrad = A.t() * (W * (yi - Ax));
    cx_mat pgrad = R.cgrad(R,x);
    ngrad = ngrad - pgrad;

    // Direction
    mat newinprod = real(dot_double(conj(ngrad), ngrad));

    if (ii == 0) {
      cx_mat ddir = ngrad;
      double gamma = 0.0;
    }
    else {
      if (abs(oldinprod) < 1e10) { gamma = 0.0;}
      else { gamma = newinprod / oldinprod; }

      ddir = ngrad + gamma * ddir;
    }

    cx_mat oldgrad = ngrad;
    mat oldinprod = newinprod;

   // Check if descent direction
    if (real(dot_double(ddir), ngrad) <0) {
     cout << " Warning descent direction not negative" << endl;
     return ;
    }

   // Step size in search direction

    cx_mat Adir = A * ddir;

    mat dAWAd = real(dot_double(conj(Adir),W*Adir));
    mat dAWr = real(Adir.t() *(W*(yi -Ax)));
    double step = 0.0;

    for (unsigned int j(0); j<2; ++j)
    {
      double pdenom = dot_double(abs(ddir)%abs(ddir) , R.denom(R,x + step*ddir));
      double denom = dAWAd + pednom;

      if (abs(denom) < 1e-10) { double step = 0.0;}
      else {  }

      pgrad = R.cgrad(R,x+step*ddir);
      pdot = real(dot_double(conj(ddir),pgrad));
      step = step - (-dAWr + step * dAWAd +pdot) / denom;
     }

    // Check downhill direction
    if (step < 0) {
      cout <<"Warning downhill"<<endl;
    }

    // Update
    Ax = Ax + step * Adir;
    x = x + step * ddir;

  }

}

cx_mat dot_double(cx_mat A, cx_mat B)
{
  return sum(A%B);
}

