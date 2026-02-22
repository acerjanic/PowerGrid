/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [solve_pwls_pcg.hpp]

    Synopsis    [Object implementing a preconditioned weighted least squares
                    solver using preconditioned conjugate gradient.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/19/2016]

 *****************************************************************************/

#ifndef POWERGRID_SOLVE_PWLS_PCG_HPP_
#define POWERGRID_SOLVE_PWLS_PCG_HPP_

#include <cstdlib>
#include <chrono>
#include "Core/PGLog.hpp"

#ifdef METAL_COMPUTE
#include "Core/pgCol.hpp"
#endif

using namespace arma;

template <typename T1>
inline complex<T1> dot_double(const Col<complex<T1>>& A,
                              const Col<complex<T1>>& B) {
  complex<T1> sumReturn = accu(A % B);
  return sumReturn;
}

template <typename T1>
inline T1 norm_grad(const Col<complex<T1>> &g, const Col<complex<T1>> &yi,
                    const Col<T1> &W) {
  T1 normGrad = conv_to<T1>::from(norm(g) / real(trans(yi) * (W % yi)));
  return normGrad;
}

template <typename T1, typename Tobj, typename Robj>
Col<complex<T1>> solve_pwls_pcg(const Col<complex<T1>> &xInitial, Tobj const &A,
                                Col<T1> const &W, Col<complex<T1>> const &yi,
                                Robj const &R, uword niter) {
  typedef complex<T1> CxT1;

#ifdef METAL_COMPUTE
  if constexpr (std::is_same<T1, float>::value) {
    // Metal path: pgCol with GPU-dispatched vector algebra
    PG_INFO("Starting PCG solver (Metal path): {} iterations requested", niter);

    // Convert inputs to pgCol
    pgCol<pgComplex<T1>> x_pg(xInitial);
    pgCol<pgComplex<T1>> yi_pg(yi);
    pgCol<T1> W_pg(W);

    // Initial forward projection (pgCol overload — no arma conversion)
    pgCol<pgComplex<T1>> Ax_pg = A * x_pg;
    if (Ax_pg.has_nan())
      PG_WARN("Ax has NaN after initial forward projection");

    pgComplex<T1> oldinprod(0, 0);
    pgComplex<T1> gamma(0, 0);
    pgCol<pgComplex<T1>> ddir_pg;
    pgCol<pgComplex<T1>> Adir_pg;
    pgComplex<T1> newinprod;

    PG_DEBUG("Entering PCG iteration loop (Metal)");
    for (unsigned int ii = 0; ii < niter; ii++) {
      // Compute negative gradient: ngrad = A' * (W .* (yi - Ax))
      pgCol<pgComplex<T1>> residual = yi_pg - Ax_pg;        // Metal cvec_sub
      pgCol<pgComplex<T1>> Wresidual = W_pg % residual;     // Metal rvec_cmul
      pgCol<pgComplex<T1>> ngrad_pg = A / Wresidual; // pgCol overload

      if (ngrad_pg.has_nan())
        PG_WARN("ngrad has NaN at iteration {}", ii);

      if (norm_grad<T1>(ngrad_pg.getArma(), yi, W) < 1e-10) {
        PG_INFO("Terminating early: zero gradient at iteration {}", ii);
        return x_pg.getArma();
      }

      // Subtract regularizer gradient
      pgCol<pgComplex<T1>> rgrad_pg(R.Gradient(x_pg.getArma())); // arma boundary
      ngrad_pg -= rgrad_pg;                                  // Metal cvec_sub

      // Direction (conjugate gradient update)
      newinprod = cdot(ngrad_pg, ngrad_pg);                  // Metal cvec_cdot
      newinprod.imag(0); // force real

      if (ii == 0) {
        ddir_pg = ngrad_pg;
      } else {
        if (abs(oldinprod) < 1e-10) {
          gamma = pgComplex<T1>(0, 0);
        } else {
          gamma = newinprod / oldinprod;
        }
        // ddir = ngrad + gamma * ddir
        pgCol<pgComplex<T1>> scaled = ddir_pg % gamma;      // Metal cvec_mul_scalar
        ddir_pg = ngrad_pg + scaled;                         // Metal cvec_add
      }

      oldinprod = newinprod;

      // Check descent direction
      pgComplex<T1> descCheck = cdot(ddir_pg, ngrad_pg);     // Metal cvec_cdot
      if (descCheck.real() < 0) {
        PG_WARN("Descent direction is not downhill at iteration {} -- stopping", ii);
        return x_pg.getArma();
      }

      // Step size in search direction
      Adir_pg = A * ddir_pg; // pgCol overload
      if (Adir_pg.has_nan())
        PG_WARN("NaN in Adir at iteration {}", ii);

      pgCol<pgComplex<T1>> WAdir_pg = W_pg % Adir_pg;       // Metal rvec_cmul
      pgComplex<T1> dAWAd = cdot(Adir_pg, WAdir_pg);        // Metal cvec_cdot
      T1 dAWAd_re = dAWAd.real();

      pgCol<pgComplex<T1>> Wresid2 = W_pg % (yi_pg - Ax_pg); // Metal rvec_cmul
      pgComplex<T1> dAWr_cx = cdot(Adir_pg, Wresid2);        // Metal cvec_cdot
      T1 dAWr_re = dAWr_cx.real();

      pgComplex<T1> step_pg(0, 0);

      for (unsigned int j = 0; j < 3; j++) {
        // Compute x + step * ddir for regularizer
        pgCol<pgComplex<T1>> xstep = x_pg + (ddir_pg % step_pg); // Metal ops

        CxT1 pdenom = R.Denom(ddir_pg.getArma(), xstep.getArma()); // arma boundary
        T1 denom_re = dAWAd_re + std::real(pdenom);
        T1 denom_im = std::imag(pdenom);

        if (std::abs(denom_re) < 1e-20 || std::isinf(denom_re) || std::isnan(denom_re)) {
          T1 n = norm(ngrad_pg);
          if (n == 0) {
            PG_INFO("Found exact solution");
            return x_pg.getArma();
          } else {
            PG_WARN("inf denom (denom_re={})", denom_re);
            return x_pg.getArma();
          }
        }

        Col<CxT1> pgrad_arma = R.Gradient(xstep.getArma());  // arma boundary
        pgCol<pgComplex<T1>> pgrad_pg(pgrad_arma);
        pgComplex<T1> pdot = cdot(ddir_pg, pgrad_pg);         // Metal cvec_cdot
        T1 pdot_re = pdot.real();

        T1 step_update = (-dAWr_re + step_pg.real() * dAWAd_re + pdot_re) / denom_re;
        step_pg.real(step_pg.real() - step_update);
      }

      // Check downhill direction
      if (step_pg.real() < 0) {
        PG_WARN("Step size is negative (downhill?) at iteration {}", ii);
      }

      // Update: Ax += step * Adir, x += step * ddir
      Ax_pg += Adir_pg % step_pg;                            // Metal cvec_mul_scalar + add
      x_pg += ddir_pg % step_pg;                             // Metal cvec_mul_scalar + add

      T1 errNorm = norm(yi_pg - Ax_pg);
      PG_DEBUG("PCG iteration {}/{}: error norm = {}", ii + 1, niter, errNorm);
    }
    return x_pg.getArma();
  }
#endif

  // Armadillo path (double, or non-Metal builds)
  PG_INFO("Starting PCG solver: {} iterations requested", niter);
  Col<CxT1> Ax = A * xInitial;
  if (Ax.has_nan())
    PG_WARN("Ax has NaN after initial forward projection");
  Col<CxT1> x = xInitial;
  CxT1 oldinprod = 0;
  CxT1 gamma = 0.0;
  Col<CxT1> ddir;
  Col<CxT1> Adir;
  CxT1 dAWAd;
  CxT1 dAWr;
  CxT1 pdenom;
  CxT1 denom;

  Col<CxT1> ngrad;
  Col<CxT1> pgrad;
  CxT1 pdot;
  Col<CxT1> WAdir;
  Col<CxT1> stepIntermediate;
  CxT1 step;
  CxT1 newinprod;

  PG_DEBUG("Entering PCG iteration loop");
  for (unsigned int ii = 0; ii < niter; ii++) {
    // Compute negative gradient

    ngrad = A / (W % (yi - Ax));
    if(ngrad.has_nan())
      PG_WARN("ngrad has NaN at iteration {}", ii);


    if (norm_grad<T1>(ngrad, yi, W) < 1e-10) {
      PG_INFO("Terminating early: zero gradient at iteration {}", ii);
      return x;
    }
    ngrad -= R.Gradient(x);

    // Direction
    newinprod = real(cdot(ngrad,ngrad));
    if (ii == 0) {
      ddir = ngrad;

    } else {
      if (std::abs(oldinprod) < 1e-10) {
        gamma = 0.0;
      } else {
        gamma = newinprod / oldinprod;
      }

      ddir = ngrad + gamma * ddir;
    }

    Col<CxT1> oldgrad = ngrad;
    oldinprod = newinprod;

    // Check if descent direction
    if (real(cdot(ddir, ngrad)) < 0) {
      PG_WARN("Descent direction is not downhill at iteration {} -- stopping", ii);
      return x;
    }

    // Step size in search direction
    Adir = A * ddir;
    if (Adir.has_nan())
      PG_WARN("NaN in Adir at iteration {}", ii);

    WAdir = W % Adir;
    dAWAd = as_scalar(real(cdot(Adir, WAdir)));
    dAWr  = as_scalar(real(Adir.t() * (W % (yi - Ax))));

    step = 0.0;

    for (unsigned int j = 0; j < 3; j++) {
      pdenom = R.Denom(ddir, x + step * ddir);
      denom = dAWAd + pdenom;
      if( denom != denom)
        PG_WARN("denom has NaN at iteration {}", ii);

      if (std::abs(denom) < 1e-20 || std::isinf(std::abs(denom)) || std::isnan(std::abs(denom))) {
        if (norm(ngrad, 2) == 0) {
          PG_INFO("Found exact solution");
          return x;
        } else {
          PG_WARN("inf denom (denom={})", std::abs(denom));
          return x;
        }
      }

      pgrad = R.Gradient(x + step * ddir);
      pdot = real(cdot(ddir, pgrad));

      stepIntermediate = (-dAWr + step * dAWAd + pdot) / denom;
      step -= as_scalar(stepIntermediate);
    }

    if (as_scalar(real(step)) < 0) {
      PG_WARN("Step size is negative (downhill?) at iteration {}", ii);
    }

    // Update
    Ax += step * Adir;
    x += (step * ddir);
    PG_DEBUG("PCG iteration {}/{}: error norm = {}", ii + 1, niter, norm(yi - Ax, 2));

  }
  return x;
}

#endif /* POWERGRID_SOLVE_PWLS_PCG_HPP_ */
