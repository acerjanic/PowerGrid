#include "catch.hpp"

#include "../Core/PGIncludes.h"
#include "../Operators/Gdft.h"
#include "../Operators/Gnufft.h"
#include "../Operators/SENSE.h"
#include "../Solvers/QuadPenalty.h"
#include "../Solvers/solve_pwls_pcg.hpp"

#include <cmath>
#include <complex>

using namespace arma;

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

// Generate Cartesian image-space and k-space coordinates for an Nx x Ny grid.
//
// Image coords: ix = jj/Nx, iy = ii/Ny  (uncentered, range [0, (N-1)/N])
//   This gives G^H*G = n1*I (standard DFT orthogonality), so PCG converges
//   in one step when the system is perfectly conditioned.
//
// Gdft kspace:   kx = jj (integers 0..Nx-1), ky = ii (integers 0..Ny-1)
//   Phase = 2*pi * kx * ix = 2*pi * jj * mm / Nx — standard DFT.
//
// Gnufft kspace: centered integers [-Nx/2, Nx/2-1] as Gnufft expects.
//   The gridding shift maps kx+Nx/2 → FFT index jj, producing
//   the same DFT as Gdft when ix = m/Nx.
static void makeCartesianCoords(uword Nx, uword Ny,
    Col<float>& ix, Col<float>& iy, Col<float>& iz,
    Col<float>& kxGdft, Col<float>& kyGdft, Col<float>& kzGdft,
    Col<float>& kxGnufft, Col<float>& kyGnufft, Col<float>& kzGnufft)
{
    uword n = Nx * Ny;
    ix.set_size(n); iy.set_size(n); iz.set_size(n);
    kxGdft.set_size(n); kyGdft.set_size(n); kzGdft.set_size(n);
    kxGnufft.set_size(n); kyGnufft.set_size(n); kzGnufft.set_size(n);
    iz.zeros(); kzGdft.zeros(); kzGnufft.zeros();

    // Column-major: y (ii) varies fastest, x (jj) next — matches vectorise(Cube)
    for (uword jj = 0; jj < Nx; jj++) {
        for (uword ii = 0; ii < Ny; ii++) {
            uword idx = ii + jj * Ny;
            // Image space: uncentered positions in [0, (N-1)/N]
            ix(idx) = (float)jj / (float)Nx;
            iy(idx) = (float)ii / (float)Ny;
            // Gdft: uncentered integer k-space [0, N-1]
            // Phase: kxGdft * ix = jj * mm/Nx → standard DFT, G^H*G = n1*I
            kxGdft(idx) = (float)jj;
            kyGdft(idx) = (float)ii;
            // Gnufft: centered integer k-space [-N/2, N/2-1]
            // gridding shift: kxGnufft + Nx/2 = jj → same DFT frequency as Gdft
            kxGnufft(idx) = (float)jj - (float)Nx * 0.5f;
            kyGnufft(idx) = (float)ii - (float)Ny * 0.5f;
        }
    }
}

// Deterministic complex test vector: re[i] = i+1, im[i] = n-i
static Col<cx_float> makeTestVec(uword n) {
    Col<cx_float> v(n);
    for (uword i = 0; i < n; i++) {
        v(i) = cx_float(static_cast<float>(i + 1),
                        static_cast<float>(n - i));
    }
    return v;
}

// ---------------------------------------------------------------------------
// Test 1: Gdft adjoint property  <G*x, y> == <x, G^H*y>
// ---------------------------------------------------------------------------

TEST_CASE("Gdft: adjoint property <G*x, y> == <x, G^H*y>", "[Gdft adjoint]") {

    uword Nx = 8, Ny = 8;
    uword n1 = Nx * Ny;   // image size
    uword n2 = Nx * Ny;   // kspace samples (Cartesian, fully sampled)

    Col<float> ix, iy, iz, kxG, kyG, kzG, kxN, kyN, kzN;
    makeCartesianCoords(Nx, Ny, ix, iy, iz, kxG, kyG, kzG, kxN, kyN, kzN);

    Col<float> FM(n1, fill::zeros);
    Col<float> t(n2, fill::zeros);

    Gdft<float> G(n1, n2, kxG, kyG, kzG, ix, iy, iz, FM, t);

    Col<cx_float> x = makeTestVec(n1);
    Col<cx_float> y = makeTestVec(n2);
    // Offset y so it's distinct from x
    for (uword i = 0; i < n2; i++) {
        y(i) = cx_float(static_cast<float>(2 * i + 1),
                        -static_cast<float>(i));
    }

    SECTION("Adjoint property holds: |<G*x,y> - <x,G^H*y>| / |<G*x,y>| < 1e-3") {
        Col<cx_float> Gx  = G * x;
        Col<cx_float> GHy = G / y;

        cx_float lhs = cdot(Gx, y);
        cx_float rhs = cdot(x, GHy);

        float relError = std::abs(lhs - rhs) / std::abs(lhs);
        REQUIRE(relError < 1e-3f);
    }
}

// ---------------------------------------------------------------------------
// Test 2: Gnufft adjoint property
// ---------------------------------------------------------------------------

TEST_CASE("Gnufft: adjoint property <G*x, y> == <x, G^H*y>", "[Gnufft adjoint]") {

    uword Nx = 8, Ny = 8;
    uword n1 = Nx * Ny;
    uword n2 = Nx * Ny;
    float gridOS = 2.0f;

    Col<float> ix, iy, iz, kxG, kyG, kzG, kxN, kyN, kzN;
    makeCartesianCoords(Nx, Ny, ix, iy, iz, kxG, kyG, kzG, kxN, kyN, kzN);

    Gnufft<float> G(n2, gridOS, Nx, Ny, 1, kxN, kyN, kzN, ix, iy, iz);

    Col<cx_float> x = makeTestVec(n1);
    Col<cx_float> y(n2);
    for (uword i = 0; i < n2; i++) {
        y(i) = cx_float(static_cast<float>(2 * i + 1),
                        -static_cast<float>(i));
    }

    SECTION("Adjoint property holds: |<G*x,y> - <x,G^H*y>| / |<G*x,y>| < 1e-3") {
        Col<cx_float> Gx  = G * x;
        Col<cx_float> GHy = G / y;

        cx_float lhs = cdot(Gx, y);
        cx_float rhs = cdot(x, GHy);

        float relError = std::abs(lhs - rhs) / std::abs(lhs);
        REQUIRE(relError < 1e-3f);
    }
}

// ---------------------------------------------------------------------------
// Test 3: SENSE<Gdft> adjoint property
// ---------------------------------------------------------------------------

TEST_CASE("SENSE<Gdft>: adjoint property <S*x, y> == <x, S^H*y>", "[SENSE adjoint]") {

    uword Nx = 8, Ny = 8;
    uword n1 = Nx * Ny;  // image size
    uword n2 = Nx * Ny;  // kspace samples per coil
    uword nc = 2;        // number of coils

    Col<float> ix, iy, iz, kxG, kyG, kzG, kxN, kyN, kzN;
    makeCartesianCoords(Nx, Ny, ix, iy, iz, kxG, kyG, kzG, kxN, kyN, kzN);

    Col<float> FM(n1, fill::zeros);
    Col<float> t(n2, fill::zeros);
    Gdft<float> Gd(n1, n2, kxG, kyG, kzG, ix, iy, iz, FM, t);

    // Constant sensitivity maps: coil 0 = all 1+0i, coil 1 = all 0+1i
    Col<cx_float> senseMap(n1 * nc);
    senseMap.subvec(0,    n1 - 1).fill(cx_float(1.0f, 0.0f));
    senseMap.subvec(n1, 2*n1 - 1).fill(cx_float(0.0f, 1.0f));

    SENSE<float, Gdft<float>> S(Gd, senseMap, n2, n1, nc);

    Col<cx_float> x = makeTestVec(n1);
    Col<cx_float> y(n2 * nc);
    for (uword i = 0; i < n2 * nc; i++) {
        y(i) = cx_float(static_cast<float>(2 * i + 1),
                        -static_cast<float>(i));
    }

    SECTION("Adjoint property holds: |<S*x,y> - <x,S^H*y>| / |<S*x,y>| < 1e-3") {
        Col<cx_float> Sx  = S * x;
        Col<cx_float> SHy = S / y;

        cx_float lhs = cdot(Sx, y);
        cx_float rhs = cdot(x, SHy);

        float relError = std::abs(lhs - rhs) / std::abs(lhs);
        REQUIRE(relError < 1e-3f);
    }
}

// ---------------------------------------------------------------------------
// Test 4: PCG solver convergence on noiseless Cartesian Gdft problem
// ---------------------------------------------------------------------------

TEST_CASE("PCG solver converges on noiseless Cartesian Gdft problem", "[PCG convergence]") {

    uword Nx = 4, Ny = 4;
    uword n1 = Nx * Ny;  // image size = 16
    uword n2 = Nx * Ny;  // kspace samples = 16

    Col<float> ix, iy, iz, kxG, kyG, kzG, kxN, kyN, kzN;
    makeCartesianCoords(Nx, Ny, ix, iy, iz, kxG, kyG, kzG, kxN, kyN, kzN);

    Col<float> FM(n1, fill::zeros);
    Col<float> t(n2, fill::zeros);
    Gdft<float> G(n1, n2, kxG, kyG, kzG, ix, iy, iz, FM, t);

    // Known true image: unit impulse at index 0 (DC component)
    Col<cx_float> x_true(n1, fill::zeros);
    x_true(0) = cx_float(1.0f, 0.0f);

    // Noiseless k-space data
    Col<cx_float> yi = G * x_true;

    // All-ones data weights (unweighted least squares)
    Col<float> W(n2, fill::ones);

    // Zero regularization so PCG converges to the least-squares solution
    QuadPenalty<float> R(Nx, Ny, 1, 0.0f);

    // Starting from zero, run PCG for up to 40 iterations
    Col<cx_float> x0(n1, fill::zeros);

    SECTION("PCG recovers x_true from noiseless data with relative error < 1e-2") {
        Col<cx_float> xhat = solve_pwls_pcg<float>(x0, G, W, yi, R, 40);
        float relError = norm(xhat - x_true) / norm(x_true);
        REQUIRE(relError < 1e-2f);
    }
}
