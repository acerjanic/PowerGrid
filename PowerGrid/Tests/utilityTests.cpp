#include "catch.hpp"

#include "../Core/PGIncludes.h"
#include "../Gridding/griddingSupport.h"
#include "../FFT/fftCPU.h"
#include "../Solvers/QuadPenalty.h"

#include <complex>
#include <cmath>
#include <vector>

using namespace arma;

// ===================================================================
// bessi0: Modified Bessel Function of the First Kind I_0(x)
// ===================================================================

TEST_CASE("bessi0: Modified Bessel I_0 function", "[bessi0]") {

    SECTION("bessi0(0) == 1.0") {
        REQUIRE(bessi0(0.0f) == Approx(1.0f));
    }

    SECTION("bessi0(1) matches tabulated value 1.2661") {
        // Known value: I_0(1) = 1.2660658...
        REQUIRE(bessi0(1.0f) == Approx(1.2660658f).epsilon(1e-4f));
    }

    SECTION("bessi0(2) matches tabulated value 2.2796") {
        // Known value: I_0(2) = 2.2795853...
        REQUIRE(bessi0(2.0f) == Approx(2.2795853f).epsilon(1e-4f));
    }

    SECTION("bessi0 is symmetric: bessi0(-x) == bessi0(x)") {
        REQUIRE(bessi0(-1.5f) == Approx(bessi0(1.5f)));
        REQUIRE(bessi0(-3.0f) == Approx(bessi0(3.0f)));
    }

    SECTION("bessi0 is monotone non-decreasing for x >= 0") {
        for (float x = 0.0f; x < 5.0f; x += 0.5f) {
            REQUIRE(bessi0(x + 0.5f) >= bessi0(x));
        }
    }

    // Double precision
    SECTION("bessi0<double>(0) == 1.0") {
        REQUIRE(bessi0(0.0) == Approx(1.0));
    }

    SECTION("bessi0<double>(1) matches tabulated value") {
        // The polynomial uses float-precision coefficients so ~1e-6 relative accuracy
        REQUIRE(bessi0(1.0) == Approx(1.2660658777520082).epsilon(1e-5));
    }
}

// ===================================================================
// Kaiser-Bessel Lookup Table
// ===================================================================

TEST_CASE("Kaiser-Bessel LUT construction and accuracy", "[kbLUT]") {

    float width = 4.0f;
    // beta value corresponding to width=4 (standard gridding parameter)
    float beta = MRI_PI * std::sqrt((width * width / 4.0f) * (4.0f / 1.25f - 0.8f));
    float* LUT = nullptr;
    uword sizeLUT = 0;

    calculateLUT(beta, width, LUT, sizeLUT);

    REQUIRE(LUT != nullptr);
    REQUIRE(sizeLUT > 0);

    SECTION("LUT peak at dist=0 is positive") {
        float peak = kernel_value_LUT(0.0f, LUT, sizeLUT, width);
        REQUIRE(peak > 0.0f);
    }

    SECTION("LUT value decreases from center to edge") {
        float v0 = kernel_value_LUT(0.0f, LUT, sizeLUT, width);
        float v1 = kernel_value_LUT(1.0f, LUT, sizeLUT, width);
        REQUIRE(v0 > v1);
    }

    SECTION("LUT value at half-width (support edge) is near zero") {
        float vEdge = kernel_value_LUT(width / 2.0f, LUT, sizeLUT, width);
        REQUIRE(vEdge == Approx(0.0f).margin(0.05f));
    }

    free(LUT);
}

// ===================================================================
// fftshift2 round-trip (even dimensions)
// ===================================================================

TEST_CASE("fftshift2 round-trip for even dimensions", "[fftshift]") {

    int xdim = 4, ydim = 4;
    // fftshift2 / circshift2 operate on interleaved complex data:
    // each "pixel" is 2 floats (re, im), so total floats = 2*xdim*ydim
    int N = xdim * ydim;
    int Nc = 2 * N;  // interleaved-complex array size

    std::vector<float> data(Nc), shifted(Nc), recovered(Nc);

    // Fill with distinct values
    for (int i = 0; i < Nc; i++) {
        data[i] = static_cast<float>(i + 1);
    }

    SECTION("Two fftshift2 applications recover original (even dims)") {
        fftshift2(shifted.data(), data.data(), xdim, ydim);
        fftshift2(recovered.data(), shifted.data(), xdim, ydim);
        for (int i = 0; i < Nc; i++) {
            REQUIRE(recovered[i] == Approx(data[i]));
        }
    }

    SECTION("fftshift2 moves DC (index 0) to center") {
        // For a 4x4 complex array, DC at (x=0,y=0) should move to (x=2,y=2)
        // Real part of pixel (x,y) is at index 2*(x + y*xdim)
        std::vector<float> impulse(Nc, 0.0f);
        impulse[0] = 1.0f;   // real part of (x=0, y=0)
        std::vector<float> out(Nc, 0.0f);
        fftshift2(out.data(), impulse.data(), xdim, ydim);
        // After shift: (0,0) moves to (2,2); real part index = 2*(2 + 2*4) = 20
        REQUIRE(out[20] == Approx(1.0f));
    }
}

TEST_CASE("fftshift3 round-trip for even dimensions", "[fftshift]") {

    int xdim = 4, ydim = 4, zdim = 2;
    // Interleaved-complex: 2 floats per voxel
    int N = xdim * ydim * zdim;
    int Nc = 2 * N;

    std::vector<float> data(Nc), shifted(Nc), recovered(Nc);

    for (int i = 0; i < Nc; i++) {
        data[i] = static_cast<float>(i + 1);
    }

    SECTION("Two fftshift3 applications recover original (even dims)") {
        fftshift3(shifted.data(), data.data(), xdim, ydim, zdim);
        fftshift3(recovered.data(), shifted.data(), xdim, ydim, zdim);
        for (int i = 0; i < Nc; i++) {
            REQUIRE(recovered[i] == Approx(data[i]));
        }
    }
}

// ===================================================================
// fftCPU: 2D FFT round-trip via forward + inverse
// ===================================================================

TEST_CASE("fft2dCPU forward/inverse round-trip", "[fftCPU]") {

    uword Nx = 8, Ny = 8;
    uword count = 2 * Nx * Ny;  // interleaved complex: [re0, im0, re1, im1, ...]

    std::vector<float> data(count, 0.0f);
    std::vector<float> original(count, 0.0f);

    // Set a single real impulse at origin: re[0,0] = 1
    data[0] = 1.0f;
    original[0] = 1.0f;

    SECTION("Forward then inverse recovers original signal (with 1/N^2 normalization)") {
        fft2dCPU<float>(data.data(), Nx, Ny);
        ifft2dCPU<float>(data.data(), Nx, Ny);

        float scale = 1.0f / static_cast<float>(Nx * Ny);
        for (uword i = 0; i < count; i++) {
            REQUIRE(data[i] * scale == Approx(original[i]).margin(1e-5f));
        }
    }

    SECTION("FFT of impulse at origin produces flat (all-ones real) spectrum") {
        // Reset
        std::fill(data.begin(), data.end(), 0.0f);
        data[0] = 1.0f;

        fft2dCPU<float>(data.data(), Nx, Ny);

        // All real parts should be 1.0, all imaginary parts 0.0
        for (uword k = 0; k < Nx * Ny; k++) {
            REQUIRE(data[2 * k]     == Approx(1.0f).margin(1e-5f));  // real
            REQUIRE(data[2 * k + 1] == Approx(0.0f).margin(1e-5f));  // imag
        }
    }
}

// ===================================================================
// Difference Operator Adjoint: <Cd*x, y> == <x, Ctd*y>
// ===================================================================

TEST_CASE("Robject: Cd and Ctd satisfy the adjoint property", "[Robject adjoint]") {

    uword Nx = 4, Ny = 4, Nz = 1;
    float beta = 0.01f;
    QuadPenalty<float> R(Nx, Ny, Nz, beta);

    uword n = Nx * Ny * Nz;

    // Deterministic test vectors (avoid random seed dependency)
    Col<cx_float> x(n), y(n);
    for (uword i = 0; i < n; i++) {
        x(i) = cx_float(static_cast<float>(i + 1),
                         static_cast<float>(n - i));
        y(i) = cx_float(static_cast<float>(2 * i + 1),
                         -static_cast<float>(i));
    }

    for (uword dim = 0; dim < 2; dim++) {
        std::string label = "dim=" + std::to_string(dim);
        SECTION("<Cd(x), y> == <x, Ctd(y)> for " + label) {
            Col<cx_float> Cdx  = R.Cd(x, dim);
            Col<cx_float> Ctdy = R.Ctd(y, dim);

            cx_float lhs = dot(Cdx, y);
            cx_float rhs = dot(x, Ctdy);

            REQUIRE(lhs.real() == Approx(rhs.real()).epsilon(1e-4f));
            REQUIRE(lhs.imag() == Approx(rhs.imag()).epsilon(1e-4f));
        }
    }
}

TEST_CASE("QuadPenalty: Gradient consistency via finite differences", "[Robject adjoint]") {

    uword Nx = 4, Ny = 4, Nz = 1;
    float beta = 1.0f;
    QuadPenalty<float> R(Nx, Ny, Nz, beta);

    uword n = Nx * Ny * Nz;

    // Simple test image
    Col<cx_float> x(n);
    for (uword i = 0; i < n; i++) {
        x(i) = cx_float(static_cast<float>(i % 4), 0.0f);
    }

    SECTION("Finite difference approximation matches Gradient() for one component") {
        float eps = 1e-3f;
        Col<cx_float> grad = R.Gradient(x);

        // Check first element: (Penalty(x + eps*e0) - Penalty(x - eps*e0)) / (2*eps) == grad(0).real()
        Col<cx_float> xPlus  = x;  xPlus(0)  += cx_float(eps, 0.0f);
        Col<cx_float> xMinus = x;  xMinus(0) -= cx_float(eps, 0.0f);

        float fd = (R.Penalty(xPlus) - R.Penalty(xMinus)) / (2.0f * eps);
        REQUIRE(fd == Approx(grad(0).real()).epsilon(1e-2f));
    }
}
