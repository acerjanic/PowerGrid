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

template<typename T1, typename Tobj, typename Robj>
T1 pwls_pcg1(T1& x, Tobj const& A,T1 const& W, T1 const& yi, Robj const& R, uword niter)
{
    
    // Initialize projection
    
    T1 Ax = A*x;
    double oldinprod = 0;
    double gamma = 0.0;
    T1 ddir;
    T1 Adir;
    double dAWAd;
    T1 dAWr;
    double pdenom;
    double denom;
    
    T1 ngrad;
    T1 pgrad;
    double pdot;
    T1 cngrad;
    T1 WAdir;
    T1 proj;
    T1 stepIntermediate;
    double step;
    for (unsigned int ii = 0; ii < niter; ii++)
    {
        // Compute negative gradient
        ngrad = A / (W * (yi - Ax));
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
        
        T1 oldgrad = ngrad;
        oldinprod = newinprod;
        
        // Check if descent direction
        if (real(dot_double(ddir, ngrad)) <0) {
            cout << " Warning descent direction not negative" << endl;
            return x;
        }
        
        // Step size in search direction
        
        Adir = A * ddir;
        WAdir = W*Adir;
        dAWAd = real(dot_double(conj(Adir).eval(),WAdir));
        proj = Adir.t() *(W*(yi -Ax));
        dAWr = conv_to<double>::from(real(proj).eval());
        step = 0.0;
        
        for (unsigned int j =0; j<2; j++)
        {
            
            pdenom = real(dot_double(pow(abs(ddir),2.0).eval() , R.Denom(x + step*ddir)));
            denom = dAWAd + pdenom;
            
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
        x = x + step * ddir;
        
    }
    return x;
}





#endif /* POWERGRID_PWLS_PCG1_HPP_ */
