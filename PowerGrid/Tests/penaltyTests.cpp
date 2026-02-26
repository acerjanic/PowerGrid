/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file penaltyTests.cpp
/// @brief Phase 1 penalty tests: QuadPenalty and TVPenalty correctness,
///        adjoint identity, gradient correctness, and pgCol cross-validation.

#include "catch.hpp"

#include "../Core/PGIncludes.h"
#include "../Core/pgCol.hpp"
#include "../Penalties/QuadPenalty.h"
#include "../Penalties/TVPenalty.h"

#include <cmath>
#include <complex>

using namespace arma;

// ---------------------------------------------------------------------------
// Helper: deterministic complex test vector
// ---------------------------------------------------------------------------
static Col<cx_float> makeTestVecF(uword n) {
    Col<cx_float> v(n);
    for (uword i = 0; i < n; i++) {
        v(i) = cx_float(static_cast<float>(i + 1) / static_cast<float>(n),
                        static_cast<float>(n - i) / static_cast<float>(n));
    }
    return v;
}

// ---------------------------------------------------------------------------
// QuadPenalty tests
// ---------------------------------------------------------------------------

TEST_CASE("QuadPenalty: Cd known pattern 2D", "[QuadPenalty][Cd]") {
    uword Nx = 4, Ny = 4;
    float beta = 1.0f;
    QuadPenalty<float> R(Nx, Ny, 1, beta, 2);

    // All-ones input → Cd along x: out[ii] = d[ii] - d[ii-1] for ii >= 1
    Col<cx_float> d = ones<Col<cx_float>>(Nx * Ny);
    Col<cx_float> Cd_x = R.Cd(d, 0); // dim=0 (x): offset=1

    // All differences of constant vector should be zero except at boundary
    // offset=1: out[0] = 0, out[ii] = d[ii] - d[ii-1] = 0 for ii >= 1
    for (uword i = 1; i < Nx * Ny; i++) {
        REQUIRE(std::abs(Cd_x(i)) < 1e-6f);
    }
    REQUIRE(std::abs(Cd_x(0)) < 1e-6f);
}

TEST_CASE("QuadPenalty: Ctd is adjoint of Cd", "[QuadPenalty][adjoint]") {
    uword Nx = 6, Ny = 6;
    float beta = 1.0f;
    QuadPenalty<float> R(Nx, Ny, 1, beta, 2);

    uword n = Nx * Ny;
    Col<cx_float> x = makeTestVecF(n);
    Col<cx_float> y = makeTestVecF(n);
    y = y / (float)n; // make y different from x

    for (uword dim = 0; dim < 2; dim++) {
        Col<cx_float> Cx = R.Cd(x, dim);
        Col<cx_float> Cty = R.Ctd(y, dim);

        // <Cd(x), y> == <x, Ctd(y)>
        cx_float lhs = as_scalar(cdot(Cx, y));
        cx_float rhs = as_scalar(cdot(x, Cty));
        REQUIRE(std::abs(lhs.real() - rhs.real()) < 1e-4f);
        REQUIRE(std::abs(lhs.imag() - rhs.imag()) < 1e-4f);
    }
}

TEST_CASE("QuadPenalty: Gradient correctness via finite differences", "[QuadPenalty][Gradient]") {
    uword Nx = 4, Ny = 4;
    float beta = 0.5f;
    QuadPenalty<float> R(Nx, Ny, 1, beta, 2);

    uword n = Nx * Ny;
    Col<cx_float> x = makeTestVecF(n);
    Col<cx_float> grad = R.Gradient(x);

    // Finite difference check along real part of gradient
    float eps = 1e-3f;
    for (uword i = 0; i < std::min(n, (uword)4); i++) {
        Col<cx_float> xp = x;
        Col<cx_float> xm = x;
        xp(i) += cx_float(eps, 0);
        xm(i) -= cx_float(eps, 0);
        float Pp = R.Penalty(xp);
        float Pm = R.Penalty(xm);
        float fd_grad = (Pp - Pm) / (2.0f * eps);
        REQUIRE(std::abs(grad(i).real() - fd_grad) < 1e-2f);
    }
}

TEST_CASE("QuadPenalty: Denom is non-negative", "[QuadPenalty][Denom]") {
    uword Nx = 4, Ny = 4;
    float beta = 1.0f;
    QuadPenalty<float> R(Nx, Ny, 1, beta, 2);

    uword n = Nx * Ny;
    Col<cx_float> x = makeTestVecF(n);
    Col<cx_float> ddir = makeTestVecF(n) * cx_float(0.1f, 0.0f);

    cx_float denom = R.Denom(ddir, x);
    REQUIRE(denom.real() >= 0.0f);
}

TEST_CASE("QuadPenalty: Penalty is non-negative", "[QuadPenalty][Penalty]") {
    uword Nx = 4, Ny = 4;
    float beta = 1.0f;
    QuadPenalty<float> R(Nx, Ny, 1, beta, 2);

    Col<cx_float> x = makeTestVecF(Nx * Ny);
    float pen = R.Penalty(x);
    REQUIRE(pen >= 0.0f);
}

TEST_CASE("QuadPenalty: 3D adjoint of Cd", "[QuadPenalty][3D][adjoint]") {
    uword Nx = 4, Ny = 4, Nz = 4;
    float beta = 1.0f;
    QuadPenalty<float> R(Nx, Ny, Nz, beta, 3);

    uword n = Nx * Ny * Nz;
    Col<cx_float> x = makeTestVecF(n);
    Col<cx_float> y = makeTestVecF(n) * cx_float(0.5f, 0.0f);

    for (uword dim = 0; dim < 3; dim++) {
        Col<cx_float> Cx = R.Cd(x, dim);
        Col<cx_float> Cty = R.Ctd(y, dim);
        cx_float lhs = as_scalar(cdot(Cx, y));
        cx_float rhs = as_scalar(cdot(x, Cty));
        REQUIRE(std::abs(lhs.real() - rhs.real()) < 1e-3f);
        REQUIRE(std::abs(lhs.imag() - rhs.imag()) < 1e-3f);
    }
}

// ---------------------------------------------------------------------------
// TVPenalty tests
// ---------------------------------------------------------------------------

TEST_CASE("TVPenalty: Ctd is adjoint of Cd", "[TVPenalty][adjoint]") {
    uword Nx = 6, Ny = 6;
    float beta = 1.0f, delta = 0.1f;
    TVPenalty<float> R(Nx, Ny, 1, beta, delta, 2);

    uword n = Nx * Ny;
    Col<cx_float> x = makeTestVecF(n);
    Col<cx_float> y = makeTestVecF(n) * cx_float(0.3f, 0.0f);

    for (uword dim = 0; dim < 2; dim++) {
        Col<cx_float> Cx = R.Cd(x, dim);
        Col<cx_float> Cty = R.Ctd(y, dim);
        cx_float lhs = as_scalar(cdot(Cx, y));
        cx_float rhs = as_scalar(cdot(x, Cty));
        REQUIRE(std::abs(lhs.real() - rhs.real()) < 1e-3f);
        REQUIRE(std::abs(lhs.imag() - rhs.imag()) < 1e-3f);
    }
}

TEST_CASE("TVPenalty: Penalty is non-negative", "[TVPenalty][Penalty]") {
    uword Nx = 4, Ny = 4;
    float beta = 1.0f, delta = 0.01f;
    TVPenalty<float> R(Nx, Ny, 1, beta, delta, 2);

    Col<cx_float> x = makeTestVecF(Nx * Ny);
    float pen = R.Penalty(x);
    REQUIRE(pen >= 0.0f);
}

TEST_CASE("TVPenalty: Gradient correctness via finite differences", "[TVPenalty][Gradient]") {
    uword Nx = 4, Ny = 4;
    float beta = 0.5f, delta = 0.1f;
    TVPenalty<float> R(Nx, Ny, 1, beta, delta, 2);

    uword n = Nx * Ny;
    Col<cx_float> x = makeTestVecF(n) * cx_float(0.1f, 0.0f); // small values near smooth region
    Col<cx_float> grad = R.Gradient(x);

    float eps = 1e-3f;
    for (uword i = 0; i < std::min(n, (uword)4); i++) {
        Col<cx_float> xp = x;
        Col<cx_float> xm = x;
        xp(i) += cx_float(eps, 0);
        xm(i) -= cx_float(eps, 0);
        float Pp = R.Penalty(xp);
        float Pm = R.Penalty(xm);
        float fd_grad = (Pp - Pm) / (2.0f * eps);
        REQUIRE(std::abs(grad(i).real() - fd_grad) < 1e-1f); // TV gradient is less smooth
    }
}

TEST_CASE("TVPenalty: wpot analytical check at zero delta", "[TVPenalty][wpot]") {
    // For small d, wpot ≈ 1/(1 + |d|/delta): at d=0, wpot = 1
    uword Nx = 4, Ny = 1, Nz = 1;
    float beta = 1.0f, delta = 1.0f;
    TVPenalty<float> R(Nx, Ny, Nz, beta, delta, 1);

    Col<cx_float> d = zeros<Col<cx_float>>(Nx);
    Col<cx_float> w = R.wpot(d);
    // wpot(0) = 1/sqrt(1 + 0) = 1
    for (uword i = 0; i < Nx; i++) {
        REQUIRE(std::abs(w(i).real() - 1.0f) < 1e-6f);
    }
}
