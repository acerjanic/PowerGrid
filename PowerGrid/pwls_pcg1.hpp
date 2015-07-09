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
inline
complex<double> dot_double(const T1 &A, const T2  &B)
{
    complex<double> sumReturn = sum(A % B);
    return sumReturn;
}

template<typename T1>
inline
double norm_grad(const Col<T1> &g,const Col<T1> &yi,const double &W)
{
    double normGrad = conv_to<double>::from(norm(g) / real(trans(yi) * (W * yi)));
    return normGrad;

}

template<typename T1, typename Tobj, typename Robj>
Col<T1> pwls_pcg1(const Col<T1> &xInitial, Tobj const& A,double const& W, Col<T1> const& yi, Robj const& R, uword niter)
{
    
    // Initialize projection
    cout << "Entering pwls_pcg1" << endl;
    Col<T1> Ax = A*xInitial;
    Col<T1> x = xInitial;
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
    double newinprod;
    cout << "Entering pwls_pcg1 iteration loop" << endl;
    for (unsigned int ii = 0; ii < niter; ii++)
    {
        // Compute negative gradient
        cout << "About to calculate the gradient of the cost function" << endl;
        ngrad = A / (W * (yi - Ax));
        if (norm_grad<T1>(ngrad,yi,W) < 1e-10) {
            cout << "Terminating early due to zero gradient." << endl;
            return x;
        }
        pgrad = R.Gradient(x);
        ngrad = ngrad - pgrad;
        // Direction
        cngrad = conj(ngrad);
        newinprod = real(dot_double(cngrad, ngrad));
        
        if (ii == 0) {
            ddir = ngrad;
            
        }
        else {
            if (std::abs(oldinprod) < 1e-10) {
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
        if (real(dot_double(conj(ddir), ngrad)) < 0) {
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
        
        for (unsigned int j = 0; j < 3; j++)
        {
            //pdenom = real(dot_double(pow(abs(ddir),2.0).eval(), R.Denom(x + step*ddir)));
            //pdenom = Cdir' * (R.wpot(R.wt, Cdir) .* Cdir); Original MATLAB code from pwls_pcg1.m
            pdenom = R.Denom(ddir);
            cout << " pdenom = " << pdenom << endl;
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
            savemat("/home/vagrant/test_reg.mat","pgrad",pgrad);

            pdot = real(dot_double(conj(ddir),pgrad));

            stepIntermediate = (-dAWr + step * dAWAd + pdot) / denom;
            step = step - conv_to<double>::from(stepIntermediate.eval());
        }
        
        // Check downhill direction
        if (step < 0) {
            cout << "Warning downhill?" << endl;
        }
        
        // Update
        Ax = Ax + step * Adir;
        x = x + (step * ddir);
        cout << "Iteration Complete = " << ii << endl;
    }
    return x;
}





#endif /* POWERGRID_PWLS_PCG1_HPP_ */
