/*
 * pwls_pcg1.hpp
 *
 *  Created on: Apr 5, 2015
 *      Author: Giangchau
 */

#ifndef POWERGRID_PWLS_PCG1_HPP_
#define POWERGRID_PWLS_PCG1_HPP_

using namespace arma;

template<typename T1, typename Tobj, typename Robj>
T1 pwls_pcg(T1& x, Tobj const& A,T1 const& W, T1 const& yi, Robj const& R, uword niter)
{

 // Initialize projection

  T1 Ax = A*x;

  for (unsigned int ii(0); ii < niter; +ii)
  {
    // Compute negative gradient
    T1 ngrad = A / (W * (yi - Ax));
    T1 pgrad = R.Gradient(R,x);
    ngrad = ngrad - pgrad;

    // Direction
    T1 newinprod = real(dot_double(conj(ngrad), ngrad));

    if (ii == 0) {
      T1 ddir = ngrad;
      double gamma = 0.0;
    }
    else {
      if (abs(oldinprod) < 1e10) { gamma = 0.0;}
      else { gamma = newinprod / oldinprod; }

      ddir = ngrad + gamma * ddir;
    }

    T1 oldgrad = ngrad;
    T1 oldinprod = newinprod;

   // Check if descent direction
    if (real(dot_double(ddir), ngrad) <0) {
     cout << " Warning descent direction not negative" << endl;
     return ;
    }

   // Step size in search direction

    T1 Adir = A * ddir;

    T1 dAWAd = real(dot_double(conj(Adir),W*Adir));
    T1 dAWr = real(Adir.t() *(W*(yi -Ax)));
    double step = 0.0;

    for (unsigned int j(0); j<2; ++j)
    {
      double pdenom = dot_double(abs(ddir)%abs(ddir) , R.denom(R,x + step*ddir));
      double denom = dAWAd + pednom;

      pgrad = R.Gradient(R,x+step*ddir);
      pdot = real(dot_double(conj(ddir),pgrad));
      step = step - (-dAWr + step * dAWAd +pdot) / denom;
     }

    // Check downhill direction
    if (step < 0) {
      cout <<"Warning downhill?"<<endl;
    }

    // Update
    Ax = Ax + step * Adir;
    x = x + step * ddir;

  }

}


T1 dot_double(T1 A, T1 B)
{
  return sum(A%B);
}



#endif /* POWERGRID_PWLS_PCG1_HPP_ */
