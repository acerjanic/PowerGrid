/*
 * pwls_pcg1.hpp
 *
 *  Created on: Apr 5, 2015
 *      Author: Giangchau
 */

#ifndef POWERGRID_PWLS_PCG1_HPP_
#define POWERGRID_PWLS_PCG1_HPP_

using namespace arma;

template<typename T1,typename T2>
complex<double> dot_double(T1 A, T2 B)
{
    complex<double> sumReturn = sum(A % B);
    return sumReturn;
}

template<typename T1>
double norm_grad(T1 g,T1 yi,Mat<cx_double> W)
{
    double normGrad = conv_to<double>::from(norm(g) / real(trans(yi) * (W * yi)));
    return normGrad;
}

template<typename T1, typename Tobj, typename Robj>
Col<T1> pwls_pcg1(Col<T1>& x, Tobj const& A,Mat<T1> const& W, Col<T1> const& yi, Robj const& R, uword niter)
{
    
    // Initialize projection
    
    Col<T1> Ax = A*x;
    double oldinprod = 0;
    double gamma = 0.0;
    Col<T1> ddir;
    Col<T1> Adir;
    double dAWAd;
    Col<T1> dAWr;
    double pdenom;
    double denom;
    
    Col<T1> ngrad;
    Col<T1> pgrad;
    double pdot;
    Col<T1> cngrad;
    Col<T1> WAdir;
    Col<T1> proj;
    Col<T1> stepIntermediate;
    double step;
    T1 rdenom;
    for (unsigned int ii = 0; ii < niter; ii++)
    {
        // Compute negative gradient
        ngrad = A / (W * (yi - Ax));
        if (norm_grad(ngrad,yi,W) < 1e-10) {
            cout << "Terminating early due to zero gradient." << endl;
            return x;
        }
        pgrad = R.Gradient(x);
        ngrad = ngrad - pgrad;
        // Direction
        cngrad = conj(ngrad);
        double newinprod = real(dot_double(cngrad, ngrad));
        
        if (ii == 0) {
            ddir = ngrad;
            
        }
        else {
            if (std::abs(oldinprod) < 1e10) {
                gamma = 0.0;
            }
            else {
                gamma = newinprod / oldinprod;
            }
            
            ddir = ngrad + gamma * ddir;
        }
        
        Col<T1> oldgrad = ngrad;
        oldinprod = newinprod;
        
        // Check if descent direction
        if (real(dot_double(ddir, ngrad)) < 0) {
            cout << " Warning descent direction not negative" << endl;
            return x;
        }
        
        // Step size in search direction
        
        Adir = A * ddir;
        WAdir = W * Adir;
        dAWAd = real(dot_double(conj(Adir).eval(),WAdir));
        proj = Adir.t() * (W * (yi -Ax));
        dAWr = conv_to<double>::from(real(proj).eval());
        step = 0.0;
        
        for (unsigned int j = 0; j < 2; j++)
        {
            pdenom = real(dot_double(pow(abs(ddir),2.0).eval(), R.Denom(x + step*ddir)));
            denom = dAWAd + pdenom;
            
            if (std::abs(denom) < 1e-10 || std::abs(denom) > 1e25) {
              if (real(dot_double(ngrad,ngrad)) < 1e-10) {
                cout << " Found exact solution" << endl;
                  return x;
              }
              else {
                cout << "inf denom" << endl;
                return x;
              }
            }

            pgrad = R.Gradient(x+step*ddir);
            pdot = real(dot_double(conj(ddir),pgrad));

            stepIntermediate = (-dAWr + step * dAWAd + pdot) / denom;
            step = step - conv_to<double>::from(stepIntermediate.eval());
        }
        
        // Check downhill direction
        if (step < 0) {
            cout <<"Warning downhill?"<<endl;
        }
        
        // Update
        Ax = Ax + step * Adir;
        x = x + (step * ddir);
        cout << "Iteration Complete = " << ii << endl;
    }
    return x;
}





#endif /* POWERGRID_PWLS_PCG1_HPP_ */
