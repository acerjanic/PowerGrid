#include "catch.hpp"

#ifdef _OPENACC
#include "accel.h"
#endif

#include "../Core/PGIncludes.h"
#include "../Core/pgCol.hpp"
#include "../Core/pgComplex.hpp"
#include <cmath>
#include <cstdlib>
#include <iostream>

#ifdef METAL_COMPUTE
#include "../Metal/MetalVectorOps.h"
#include "../Metal/MetalVectorOps_dispatch.hpp"
#endif

// Helper: fill pgCol<float> with pseudo-random values in [-1, 1]
static void fillRandom(pgCol<float>& v, unsigned seed) {
    std::srand(seed);
    for (arma::uword i = 0; i < v.n_elem; i++) {
        v.at(i) = (float)std::rand() / RAND_MAX * 2.0f - 1.0f;
    }
}

// Helper: fill pgCol<pgComplex<float>> with pseudo-random complex values
static void fillRandomCplx(pgCol<pgComplex<float>>& v, unsigned seed) {
    std::srand(seed);
    for (arma::uword i = 0; i < v.n_elem; i++) {
        float re = (float)std::rand() / RAND_MAX * 2.0f - 1.0f;
        float im = (float)std::rand() / RAND_MAX * 2.0f - 1.0f;
        v.at(i) = pgComplex<float>(re, im);
    }
}

// ======================================================================
// Low-level Metal API tests (call Metal functions directly)
// ======================================================================

#ifdef METAL_COMPUTE
TEST_CASE("MetalVectorOps: context singleton", "[metal_vecops]") {
    MetalVectorContext* ctx = metal_vector_get_context();
    REQUIRE(ctx != nullptr);
    // Multiple calls return same pointer
    REQUIRE(metal_vector_get_context() == ctx);
}

TEST_CASE("MetalVectorOps: real vec add/sub/mul/div", "[metal_vecops]") {
    const size_t N = 16384;
    MetalVectorContext* ctx = metal_vector_get_context();
    REQUIRE(ctx != nullptr);

    std::vector<float> A(N), B(N), C(N);
    std::srand(42);
    for (size_t i = 0; i < N; i++) {
        A[i] = (float)std::rand() / RAND_MAX * 2.0f - 1.0f;
        B[i] = (float)std::rand() / RAND_MAX * 0.9f + 0.1f; // avoid div-by-zero
    }

    // Add
    metal_vec_add(ctx, A.data(), B.data(), C.data(), N);
    for (size_t i = 0; i < N; i++) {
        REQUIRE(C[i] == Approx(A[i] + B[i]).epsilon(1e-5f));
    }

    // Sub
    metal_vec_sub(ctx, A.data(), B.data(), C.data(), N);
    for (size_t i = 0; i < N; i++) {
        REQUIRE(C[i] == Approx(A[i] - B[i]).epsilon(1e-5f));
    }

    // Mul
    metal_vec_mul(ctx, A.data(), B.data(), C.data(), N);
    for (size_t i = 0; i < N; i++) {
        REQUIRE(C[i] == Approx(A[i] * B[i]).epsilon(1e-5f));
    }

    // Div
    metal_vec_div(ctx, A.data(), B.data(), C.data(), N);
    for (size_t i = 0; i < N; i++) {
        REQUIRE(C[i] == Approx(A[i] / B[i]).epsilon(1e-5f));
    }
}

TEST_CASE("MetalVectorOps: real scalar ops", "[metal_vecops]") {
    const size_t N = 8192;
    MetalVectorContext* ctx = metal_vector_get_context();
    REQUIRE(ctx != nullptr);

    std::vector<float> A(N), C(N);
    std::srand(123);
    for (size_t i = 0; i < N; i++) {
        A[i] = (float)std::rand() / RAND_MAX * 4.0f - 2.0f;
    }

    float s = 3.14f;

    metal_vec_add_scalar(ctx, A.data(), s, C.data(), N);
    for (size_t i = 0; i < N; i++) {
        REQUIRE(C[i] == Approx(A[i] + s).epsilon(1e-5f));
    }

    metal_vec_mul_scalar(ctx, A.data(), s, C.data(), N);
    for (size_t i = 0; i < N; i++) {
        REQUIRE(C[i] == Approx(A[i] * s).epsilon(1e-5f));
    }
}

TEST_CASE("MetalVectorOps: complex vec add/sub/mul/div", "[metal_vecops]") {
    const size_t N = 8192; // complex elements
    MetalVectorContext* ctx = metal_vector_get_context();
    REQUIRE(ctx != nullptr);

    // Interleaved re/im
    std::vector<float> A(2*N), B(2*N), C(2*N);
    std::srand(99);
    for (size_t i = 0; i < 2*N; i++) {
        A[i] = (float)std::rand() / RAND_MAX * 2.0f - 1.0f;
        B[i] = (float)std::rand() / RAND_MAX * 0.9f + 0.1f;
    }

    // Complex add
    metal_cvec_add(ctx, A.data(), B.data(), C.data(), N);
    for (size_t i = 0; i < 2*N; i++) {
        REQUIRE(C[i] == Approx(A[i] + B[i]).epsilon(1e-5f));
    }

    // Complex sub
    metal_cvec_sub(ctx, A.data(), B.data(), C.data(), N);
    for (size_t i = 0; i < 2*N; i++) {
        REQUIRE(C[i] == Approx(A[i] - B[i]).epsilon(1e-5f));
    }

    // Complex mul: (ar+i*ai)(br+i*bi) = (ar*br - ai*bi) + i(ar*bi + ai*br)
    metal_cvec_mul(ctx, A.data(), B.data(), C.data(), N);
    for (size_t i = 0; i < N; i++) {
        float ar = A[2*i], ai = A[2*i+1];
        float br = B[2*i], bi = B[2*i+1];
        REQUIRE(C[2*i]   == Approx(ar*br - ai*bi).epsilon(1e-5f));
        REQUIRE(C[2*i+1] == Approx(ar*bi + ai*br).epsilon(1e-5f));
    }

    // Complex div
    metal_cvec_div(ctx, A.data(), B.data(), C.data(), N);
    for (size_t i = 0; i < N; i++) {
        float ar = A[2*i], ai = A[2*i+1];
        float br = B[2*i], bi = B[2*i+1];
        float denom = br*br + bi*bi;
        REQUIRE(C[2*i]   == Approx((ar*br + ai*bi) / denom).epsilon(1e-4f));
        REQUIRE(C[2*i+1] == Approx((ai*br - ar*bi) / denom).epsilon(1e-4f));
    }
}

TEST_CASE("MetalVectorOps: rvec_cmul", "[metal_vecops]") {
    const size_t N = 8192;
    MetalVectorContext* ctx = metal_vector_get_context();
    REQUIRE(ctx != nullptr);

    std::vector<float> W(N), X(2*N), C(2*N);
    std::srand(77);
    for (size_t i = 0; i < N; i++) {
        W[i] = (float)std::rand() / RAND_MAX * 2.0f - 1.0f;
    }
    for (size_t i = 0; i < 2*N; i++) {
        X[i] = (float)std::rand() / RAND_MAX * 2.0f - 1.0f;
    }

    metal_rvec_cmul(ctx, W.data(), X.data(), C.data(), N);
    for (size_t i = 0; i < N; i++) {
        REQUIRE(C[2*i]   == Approx(W[i] * X[2*i]).epsilon(1e-5f));
        REQUIRE(C[2*i+1] == Approx(W[i] * X[2*i+1]).epsilon(1e-5f));
    }
}

TEST_CASE("MetalVectorOps: cvec_axpy", "[metal_vecops]") {
    const size_t N = 8192;
    MetalVectorContext* ctx = metal_vector_get_context();
    REQUIRE(ctx != nullptr);

    std::vector<float> A(2*N), B(2*N), C(2*N);
    std::srand(55);
    for (size_t i = 0; i < 2*N; i++) {
        A[i] = (float)std::rand() / RAND_MAX * 2.0f - 1.0f;
        B[i] = (float)std::rand() / RAND_MAX * 2.0f - 1.0f;
    }

    float alphaRe = 0.5f, alphaIm = -0.3f;
    metal_cvec_axpy(ctx, A.data(), B.data(), C.data(), alphaRe, alphaIm, N);
    for (size_t i = 0; i < N; i++) {
        float br = B[2*i], bi = B[2*i+1];
        float pr = alphaRe * br - alphaIm * bi;
        float pi = alphaRe * bi + alphaIm * br;
        REQUIRE(C[2*i]   == Approx(A[2*i]   + pr).margin(1e-5f));
        REQUIRE(C[2*i+1] == Approx(A[2*i+1] + pi).margin(1e-5f));
    }
}

TEST_CASE("MetalVectorOps: cvec_mul_scalar", "[metal_vecops]") {
    const size_t N = 8192;
    MetalVectorContext* ctx = metal_vector_get_context();
    REQUIRE(ctx != nullptr);

    std::vector<float> A(2*N), C(2*N);
    std::srand(33);
    for (size_t i = 0; i < 2*N; i++) {
        A[i] = (float)std::rand() / RAND_MAX * 2.0f - 1.0f;
    }

    float aRe = 2.0f, aIm = -1.5f;
    metal_cvec_mul_scalar(ctx, A.data(), aRe, aIm, C.data(), N);
    for (size_t i = 0; i < N; i++) {
        float ar = A[2*i], ai = A[2*i+1];
        REQUIRE(C[2*i]   == Approx(aRe*ar - aIm*ai).margin(1e-5f));
        REQUIRE(C[2*i+1] == Approx(aRe*ai + aIm*ar).margin(1e-5f));
    }
}

TEST_CASE("MetalVectorOps: reductions", "[metal_vecops]") {
    const size_t N = 65536;
    MetalVectorContext* ctx = metal_vector_get_context();
    REQUIRE(ctx != nullptr);

    SECTION("vec_sum") {
        std::vector<float> A(N);
        float expected = 0.0f;
        std::srand(11);
        for (size_t i = 0; i < N; i++) {
            A[i] = (float)std::rand() / RAND_MAX * 2.0f - 1.0f;
            expected += A[i];
        }

        float result = metal_vec_sum(ctx, A.data(), N);
        REQUIRE(result == Approx(expected).epsilon(1e-3f));
    }

    SECTION("cvec_cdot") {
        std::vector<float> A(2*N), B(2*N);
        float expectRe = 0.0f, expectIm = 0.0f;
        std::srand(22);
        for (size_t i = 0; i < 2*N; i++) {
            A[i] = (float)std::rand() / RAND_MAX * 2.0f - 1.0f;
            B[i] = (float)std::rand() / RAND_MAX * 2.0f - 1.0f;
        }
        for (size_t i = 0; i < N; i++) {
            float ar = A[2*i], ai = A[2*i+1];
            float br = B[2*i], bi = B[2*i+1];
            // conj(A)*B = (ar-i*ai)(br+i*bi)
            expectRe += ar*br + ai*bi;
            expectIm += ar*bi - ai*br;
        }

        float outRe, outIm;
        metal_cvec_cdot(ctx, A.data(), B.data(), &outRe, &outIm, N);
        REQUIRE(outRe == Approx(expectRe).epsilon(1e-3f));
        REQUIRE(outIm == Approx(expectIm).epsilon(1e-3f));
    }

    SECTION("cvec_norm2sq") {
        std::vector<float> A(2*N);
        float expected = 0.0f;
        std::srand(33);
        for (size_t i = 0; i < 2*N; i++) {
            A[i] = (float)std::rand() / RAND_MAX * 2.0f - 1.0f;
            expected += A[i] * A[i];
        }

        float result = metal_cvec_norm2sq(ctx, A.data(), N);
        REQUIRE(result == Approx(expected).epsilon(1e-3f));
    }
}
#endif // METAL_COMPUTE

// ======================================================================
// pgCol operator tests with Metal dispatch (run regardless of METAL_COMPUTE
// since they fall through to CPU on non-Metal builds)
// ======================================================================

TEST_CASE("pgCol<float>: Metal-dispatched operators", "[pgCol_metal]") {
    // Use sizes above and below the Metal threshold (4096) to test both paths
    const arma::uword N = 16384;

    pgCol<float> A(N), B(N), C(N);
    fillRandom(A, 100);
    fillRandom(B, 200);

    SECTION("operator+ (pgCol + pgCol)") {
        C = A + B;
        for (arma::uword i = 0; i < N; i++) {
            REQUIRE(C.at(i) == Approx(A.at(i) + B.at(i)).epsilon(1e-5f));
        }
    }

    SECTION("operator- (pgCol - pgCol)") {
        C = A - B;
        for (arma::uword i = 0; i < N; i++) {
            REQUIRE(C.at(i) == Approx(A.at(i) - B.at(i)).epsilon(1e-5f));
        }
    }

    SECTION("operator% (pgCol % pgCol = element-wise mul)") {
        C = A % B;
        for (arma::uword i = 0; i < N; i++) {
            REQUIRE(C.at(i) == Approx(A.at(i) * B.at(i)).epsilon(1e-5f));
        }
    }

    SECTION("operator/ (pgCol / pgCol)") {
        // Make B positive to avoid div-by-zero
        for (arma::uword i = 0; i < N; i++) {
            B.at(i) = std::abs(B.at(i)) + 0.1f;
        }
        C = A / B;
        for (arma::uword i = 0; i < N; i++) {
            REQUIRE(C.at(i) == Approx(A.at(i) / B.at(i)).epsilon(1e-5f));
        }
    }

    SECTION("scalar add/sub/mul") {
        float s = 2.5f;
        C = A + s;
        for (arma::uword i = 0; i < N; i++) {
            REQUIRE(C.at(i) == Approx(A.at(i) + s).epsilon(1e-5f));
        }
        C = A - s;
        for (arma::uword i = 0; i < N; i++) {
            REQUIRE(C.at(i) == Approx(A.at(i) - s).epsilon(1e-5f));
        }
        C = A % s;
        for (arma::uword i = 0; i < N; i++) {
            REQUIRE(C.at(i) == Approx(A.at(i) * s).epsilon(1e-5f));
        }
    }

    SECTION("sum()") {
        float expected = 0.0f;
        for (arma::uword i = 0; i < N; i++) expected += A.at(i);
        float result = sum(A);
        REQUIRE(result == Approx(expected).epsilon(1e-3f));
    }

    SECTION("compound assignment +=, -=, %=") {
        pgCol<float> orig(N);
        for (arma::uword i = 0; i < N; i++) orig.at(i) = A.at(i);

        A += B;
        for (arma::uword i = 0; i < N; i++) {
            REQUIRE(A.at(i) == Approx(orig.at(i) + B.at(i)).epsilon(1e-5f));
        }
    }
}

TEST_CASE("pgCol<pgComplex<float>>: Metal-dispatched operators", "[pgCol_metal]") {
    const arma::uword N = 16384;

    pgCol<pgComplex<float>> A(N), B(N), C(N);
    fillRandomCplx(A, 300);
    fillRandomCplx(B, 400);

    SECTION("operator+ (complex)") {
        C = A + B;
        for (arma::uword i = 0; i < N; i++) {
            REQUIRE(C.at(i).real() == Approx(A.at(i).real() + B.at(i).real()).epsilon(1e-5f));
            REQUIRE(C.at(i).imag() == Approx(A.at(i).imag() + B.at(i).imag()).epsilon(1e-5f));
        }
    }

    SECTION("operator- (complex)") {
        C = A - B;
        for (arma::uword i = 0; i < N; i++) {
            REQUIRE(C.at(i).real() == Approx(A.at(i).real() - B.at(i).real()).epsilon(1e-5f));
            REQUIRE(C.at(i).imag() == Approx(A.at(i).imag() - B.at(i).imag()).epsilon(1e-5f));
        }
    }

    SECTION("operator% (complex element-wise mul)") {
        C = A % B;
        for (arma::uword i = 0; i < N; i++) {
            pgComplex<float> expected = A.at(i) * B.at(i);
            REQUIRE(C.at(i).real() == Approx(expected.real()).epsilon(1e-4f));
            REQUIRE(C.at(i).imag() == Approx(expected.imag()).epsilon(1e-4f));
        }
    }

    SECTION("operator/ (complex element-wise div)") {
        // Ensure B has nonzero magnitude
        for (arma::uword i = 0; i < N; i++) {
            if (std::abs(B.at(i).real()) + std::abs(B.at(i).imag()) < 0.1f) {
                B.at(i) = pgComplex<float>(1.0f, 0.5f);
            }
        }
        C = A / B;
        for (arma::uword i = 0; i < N; i++) {
            pgComplex<float> expected = A.at(i) / B.at(i);
            REQUIRE(C.at(i).real() == Approx(expected.real()).epsilon(1e-4f));
            REQUIRE(C.at(i).imag() == Approx(expected.imag()).epsilon(1e-4f));
        }
    }

    SECTION("scalar multiply (complex)") {
        pgComplex<float> s(2.0f, -1.5f);
        C = A % s;
        for (arma::uword i = 0; i < N; i++) {
            pgComplex<float> expected = A.at(i) * s;
            REQUIRE(C.at(i).real() == Approx(expected.real()).margin(1e-5f));
            REQUIRE(C.at(i).imag() == Approx(expected.imag()).margin(1e-5f));
        }
    }

    SECTION("compound assignment += (complex)") {
        pgCol<pgComplex<float>> orig(N);
        for (arma::uword i = 0; i < N; i++) orig.at(i) = A.at(i);

        A += B;
        for (arma::uword i = 0; i < N; i++) {
            REQUIRE(A.at(i).real() == Approx(orig.at(i).real() + B.at(i).real()).epsilon(1e-5f));
            REQUIRE(A.at(i).imag() == Approx(orig.at(i).imag() + B.at(i).imag()).epsilon(1e-5f));
        }
    }
}

TEST_CASE("pgCol: cdot and norm (Metal-dispatched)", "[pgCol_metal]") {
    const arma::uword N = 16384;

    pgCol<pgComplex<float>> A(N), B(N);
    fillRandomCplx(A, 500);
    fillRandomCplx(B, 600);

    SECTION("cdot vs CPU reference") {
        // CPU reference
        float refRe = 0.0f, refIm = 0.0f;
        for (arma::uword i = 0; i < N; i++) {
            float ar = A.at(i).real(), ai = A.at(i).imag();
            float br = B.at(i).real(), bi = B.at(i).imag();
            refRe += ar*br + ai*bi;
            refIm += ar*bi - ai*br;
        }

        pgComplex<float> result = cdot(A, B);
        REQUIRE(result.real() == Approx(refRe).epsilon(1e-3f));
        REQUIRE(result.imag() == Approx(refIm).epsilon(1e-3f));
    }

    SECTION("cdot Hermitian symmetry: cdot(A,B) == conj(cdot(B,A))") {
        pgComplex<float> ab = cdot(A, B);
        pgComplex<float> ba = cdot(B, A);
        REQUIRE(ab.real() == Approx(ba.real()).epsilon(1e-3f));
        REQUIRE(ab.imag() == Approx(-ba.imag()).epsilon(1e-3f));
    }

    SECTION("norm2sq(A) == real(cdot(A,A))") {
        pgComplex<float> aa = cdot(A, A);
        float normVal = norm(A);
        REQUIRE(normVal * normVal == Approx(aa.real()).epsilon(1e-3f));
        // Imaginary part of cdot(A,A) should be ~0
        REQUIRE(std::abs(aa.imag()) < 1e-2f);
    }

    SECTION("norm vs CPU reference") {
        float expected = 0.0f;
        for (arma::uword i = 0; i < N; i++) {
            expected += A.at(i).real() * A.at(i).real()
                      + A.at(i).imag() * A.at(i).imag();
        }
        expected = std::sqrt(expected);
        float result = norm(A);
        REQUIRE(result == Approx(expected).epsilon(1e-3f));
    }
}

TEST_CASE("pgCol: CPU fallback for small vectors", "[pgCol_metal]") {
    // Below kMinMetalSize threshold (4096) — should still work via CPU
    const arma::uword N = 64;

    pgCol<float> A(N), B(N);
    fillRandom(A, 700);
    fillRandom(B, 800);

    pgCol<float> C = A + B;
    for (arma::uword i = 0; i < N; i++) {
        REQUIRE(C.at(i) == Approx(A.at(i) + B.at(i)).epsilon(1e-5f));
    }

    float s = sum(A);
    float expected = 0.0f;
    for (arma::uword i = 0; i < N; i++) expected += A.at(i);
    REQUIRE(s == Approx(expected).epsilon(1e-5f));
}

TEST_CASE("pgCol<double>: CPU fallback (no Metal dispatch)", "[pgCol_metal]") {
    // double types should always fall through to CPU
    const arma::uword N = 16384;

    pgCol<double> A(N), B(N);
    std::srand(900);
    for (arma::uword i = 0; i < N; i++) {
        A.at(i) = (double)std::rand() / RAND_MAX * 2.0 - 1.0;
        B.at(i) = (double)std::rand() / RAND_MAX * 2.0 - 1.0;
    }

    pgCol<double> C = A + B;
    for (arma::uword i = 0; i < N; i++) {
        REQUIRE(C.at(i) == Approx(A.at(i) + B.at(i)));
    }

    double s = sum(A);
    double expected = 0.0;
    for (arma::uword i = 0; i < N; i++) expected += A.at(i);
    REQUIRE(s == Approx(expected));
}
