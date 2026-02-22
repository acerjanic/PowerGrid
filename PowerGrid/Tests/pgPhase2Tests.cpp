/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file pgPhase2Tests.cpp
/// @brief Phase 2 foundation tests: view mode, column view, subvec, has_nan,
///        conj(pgMat), reshape constructor, real×complex operator%.

#include "catch.hpp"

#include "../Core/PGIncludes.h"
#include "../Core/pgCol.hpp"
#include "../Core/pgMat.hpp"
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>

// Inject NaN via bit manipulation — survives -ffast-math optimizations.
static float make_nan_float() {
    uint32_t bits = 0x7FC00000u; // quiet NaN
    float val;
    memcpy(&val, &bits, sizeof(val));
    return val;
}
static double make_nan_double() {
    uint64_t bits = 0x7FF8000000000000ull; // quiet NaN
    double val;
    memcpy(&val, &bits, sizeof(val));
    return val;
}

// =========================================================================
// pgCol view mode
// =========================================================================

TEST_CASE("pgCol view: wraps external memory", "[pgCol_view]") {
    const arma::uword N = 128;
    pgCol<float> owner(N);
    for (arma::uword i = 0; i < N; i++) owner.at(i) = (float)i;

    pgCol<float> v = pgCol<float>::view(owner.memptr(), N);

    REQUIRE(v.is_view());
    REQUIRE(v.n_elem == N);
    REQUIRE(v.memptr() == owner.memptr());

    // Reading through view matches original
    for (arma::uword i = 0; i < N; i++) {
        REQUIRE(v.at(i) == (float)i);
    }
}

TEST_CASE("pgCol view: write-through modifies external memory", "[pgCol_view]") {
    const arma::uword N = 64;
    pgCol<float> owner(N);
    owner.zeros();

    pgCol<float> v = pgCol<float>::view(owner.memptr(), N);
    for (arma::uword i = 0; i < N; i++) v.at(i) = (float)(i * 2);

    // Owner's memory was modified
    for (arma::uword i = 0; i < N; i++) {
        REQUIRE(owner.at(i) == (float)(i * 2));
    }
}

TEST_CASE("pgCol view: copy from view produces owning pgCol", "[pgCol_view]") {
    const arma::uword N = 32;
    pgCol<float> owner(N);
    for (arma::uword i = 0; i < N; i++) owner.at(i) = (float)i;

    pgCol<float> v = pgCol<float>::view(owner.memptr(), N);
    pgCol<float> copy(v);  // copy constructor

    REQUIRE_FALSE(copy.is_view());
    REQUIRE(copy.n_elem == N);
    REQUIRE(copy.memptr() != owner.memptr());  // separate allocation

    // Values match
    for (arma::uword i = 0; i < N; i++) {
        REQUIRE(copy.at(i) == (float)i);
    }

    // Modifying copy doesn't affect original
    copy.at(0) = 999.0f;
    REQUIRE(owner.at(0) == 0.0f);
}

TEST_CASE("pgCol view: operator%= writes through view", "[pgCol_view]") {
    const arma::uword N = 64;
    pgCol<float> owner(N);
    owner.ones();

    pgCol<float> scale(N);
    for (arma::uword i = 0; i < N; i++) scale.at(i) = 3.0f;

    // Get a view and apply in-place multiply
    pgCol<float> v = pgCol<float>::view(owner.memptr(), N);
    v %= scale;

    // Owner should be all 3s now
    for (arma::uword i = 0; i < N; i++) {
        REQUIRE(owner.at(i) == Approx(3.0f));
    }
}

TEST_CASE("pgCol view: set_size breaks view and allocates fresh", "[pgCol_view]") {
    const arma::uword N = 32;
    pgCol<float> owner(N);
    owner.ones();

    pgCol<float> v = pgCol<float>::view(owner.memptr(), N);
    REQUIRE(v.is_view());

    v.set_size(64);  // Should break view, allocate new memory
    REQUIRE_FALSE(v.is_view());
    REQUIRE(v.n_elem == 64);
    REQUIRE(v.memptr() != owner.memptr());

    // Owner should be untouched (set_size didn't free external memory)
    for (arma::uword i = 0; i < N; i++) {
        REQUIRE(owner.at(i) == Approx(1.0f));
    }
}

TEST_CASE("pgCol view: complex view wraps external memory", "[pgCol_view]") {
    const arma::uword N = 64;
    pgCol<pgComplex<float>> owner(N);
    for (arma::uword i = 0; i < N; i++) {
        owner.at(i) = pgComplex<float>((float)i, (float)(i + 1));
    }

    pgCol<pgComplex<float>> v = pgCol<pgComplex<float>>::view(owner.memptr(), N);
    REQUIRE(v.is_view());
    REQUIRE(v.n_elem == N);
    for (arma::uword i = 0; i < N; i++) {
        REQUIRE(v.at(i).real() == Approx((float)i));
        REQUIRE(v.at(i).imag() == Approx((float)(i + 1)));
    }
}

// =========================================================================
// pgMat column view
// =========================================================================

TEST_CASE("pgMat col(): returns non-owning view", "[pgMat_col_view]") {
    const arma::uword nR = 8, nC = 4;
    pgMat<float> M(nR, nC);
    // Fill with column-major data: M(r,c) = c * 100 + r
    for (arma::uword c = 0; c < nC; c++)
        for (arma::uword r = 0; r < nR; r++)
            M.at(r, c) = (float)(c * 100 + r);

    pgCol<float> col1 = M.col(1);
    REQUIRE(col1.is_view());
    REQUIRE(col1.n_elem == nR);

    for (arma::uword r = 0; r < nR; r++) {
        REQUIRE(col1.at(r) == Approx((float)(100 + r)));
    }
}

TEST_CASE("pgMat col(): write through view modifies matrix", "[pgMat_col_view]") {
    const arma::uword nR = 8, nC = 4;
    pgMat<float> M(nR, nC);
    M.zeros();

    // Write through column view
    {
        pgCol<float> col2 = M.col(2);
        for (arma::uword r = 0; r < nR; r++) {
            col2.at(r) = (float)(r + 1);
        }
    }
    // View is destroyed, but writes persisted
    for (arma::uword r = 0; r < nR; r++) {
        REQUIRE(M.at(r, 2) == Approx((float)(r + 1)));
    }
    // Other columns unchanged
    for (arma::uword r = 0; r < nR; r++) {
        REQUIRE(M.at(r, 0) == Approx(0.0f));
        REQUIRE(M.at(r, 1) == Approx(0.0f));
        REQUIRE(M.at(r, 3) == Approx(0.0f));
    }
}

TEST_CASE("pgMat col(): operator%= modifies matrix in-place", "[pgMat_col_view]") {
    const arma::uword nR = 16, nC = 3;
    pgMat<pgComplex<float>> M(nR, nC);
    // Fill with known data
    for (arma::uword c = 0; c < nC; c++)
        for (arma::uword r = 0; r < nR; r++)
            M.at(r, c) = pgComplex<float>(1.0f, 0.0f);

    // Scale column 1 by (2 + 3i)
    pgCol<pgComplex<float>> scale(nR);
    for (arma::uword r = 0; r < nR; r++)
        scale.at(r) = pgComplex<float>(2.0f, 3.0f);

    M.col(1) %= scale;  // should modify matrix in-place

    for (arma::uword r = 0; r < nR; r++) {
        // (1+0i) * (2+3i) = (2+3i)
        REQUIRE(M.at(r, 1).real() == Approx(2.0f));
        REQUIRE(M.at(r, 1).imag() == Approx(3.0f));
    }
    // Column 0 unchanged
    for (arma::uword r = 0; r < nR; r++) {
        REQUIRE(M.at(r, 0).real() == Approx(1.0f));
        REQUIRE(M.at(r, 0).imag() == Approx(0.0f));
    }
}

TEST_CASE("pgMat set_col() writes column data", "[pgMat_col_view]") {
    const arma::uword nR = 10, nC = 3;
    pgMat<float> M(nR, nC);
    M.zeros();

    pgCol<float> src(nR);
    for (arma::uword r = 0; r < nR; r++) src.at(r) = (float)(r * 10);

    M.set_col(2, src);

    for (arma::uword r = 0; r < nR; r++) {
        REQUIRE(M.at(r, 2) == Approx((float)(r * 10)));
    }
    // Other columns still zero
    for (arma::uword r = 0; r < nR; r++) {
        REQUIRE(M.at(r, 0) == Approx(0.0f));
    }
}

TEST_CASE("pgMat col_copy() returns owning deep copy", "[pgMat_col_view]") {
    const arma::uword nR = 8, nC = 4;
    pgMat<float> M(nR, nC);
    for (arma::uword i = 0; i < M.n_elem; i++) M.at(i) = (float)i;

    pgCol<float> deep = M.col_copy(2);
    REQUIRE_FALSE(deep.is_view());
    REQUIRE(deep.n_elem == nR);
    REQUIRE(deep.memptr() != &M.memptr()[nR * 2]);

    for (arma::uword r = 0; r < nR; r++) {
        REQUIRE(deep.at(r) == Approx(M.at(r, 2)));
    }

    // Modifying deep copy doesn't affect matrix
    deep.at(0) = -999.0f;
    REQUIRE(M.at(0, 2) != Approx(-999.0f));
}

// =========================================================================
// pgCol subvec
// =========================================================================

TEST_CASE("pgCol subvec: returns view of slice", "[pgCol_subvec]") {
    const arma::uword N = 100;
    pgCol<float> owner(N);
    for (arma::uword i = 0; i < N; i++) owner.at(i) = (float)i;

    pgCol<float> sv = owner.subvec(10, 19);
    REQUIRE(sv.is_view());
    REQUIRE(sv.n_elem == 10);

    for (arma::uword i = 0; i < 10; i++) {
        REQUIRE(sv.at(i) == Approx((float)(i + 10)));
    }
}

TEST_CASE("pgCol subvec: write-through to parent", "[pgCol_subvec]") {
    const arma::uword N = 50;
    pgCol<float> owner(N);
    owner.zeros();

    pgCol<float> sv = owner.subvec(20, 29);
    for (arma::uword i = 0; i < 10; i++) sv.at(i) = 42.0f;

    // Check owner's elements
    for (arma::uword i = 0; i < N; i++) {
        if (i >= 20 && i <= 29)
            REQUIRE(owner.at(i) == Approx(42.0f));
        else
            REQUIRE(owner.at(i) == Approx(0.0f));
    }
}

TEST_CASE("pgCol subvec: complex slice", "[pgCol_subvec]") {
    const arma::uword N = 64;
    pgCol<pgComplex<float>> owner(N);
    for (arma::uword i = 0; i < N; i++)
        owner.at(i) = pgComplex<float>((float)i, -(float)i);

    pgCol<pgComplex<float>> sv = owner.subvec(32, 47);
    REQUIRE(sv.n_elem == 16);

    for (arma::uword i = 0; i < 16; i++) {
        REQUIRE(sv.at(i).real() == Approx((float)(i + 32)));
        REQUIRE(sv.at(i).imag() == Approx(-(float)(i + 32)));
    }
}

// =========================================================================
// pgCol has_nan
// =========================================================================

TEST_CASE("pgCol has_nan: no NaN in clean data", "[pgCol_has_nan]") {
    pgCol<float> a(100);
    a.ones();
    REQUIRE_FALSE(a.has_nan());
}

TEST_CASE("pgCol has_nan: detects NaN in float", "[pgCol_has_nan]") {
    pgCol<float> a(100);
    a.ones();
    // Use bit-level NaN injection to survive -ffast-math
    float nanVal = make_nan_float();
    memcpy(&a.memptr()[50], &nanVal, sizeof(float));
    REQUIRE(a.has_nan());
}

TEST_CASE("pgCol has_nan: detects NaN in double", "[pgCol_has_nan]") {
    pgCol<double> a(100);
    a.ones();
    double nanVal = make_nan_double();
    memcpy(&a.memptr()[0], &nanVal, sizeof(double));
    REQUIRE(a.has_nan());
}

TEST_CASE("pgCol has_nan: detects NaN in pgComplex<float> real part", "[pgCol_has_nan]") {
    pgCol<pgComplex<float>> a(100);
    for (arma::uword i = 0; i < 100; i++)
        a.at(i) = pgComplex<float>(1.0f, 2.0f);
    REQUIRE_FALSE(a.has_nan());

    // Inject NaN into the real part of element 42 (first float of interleaved pair)
    float nanVal = make_nan_float();
    float* rawPtr = reinterpret_cast<float*>(a.memptr());
    memcpy(&rawPtr[42 * 2], &nanVal, sizeof(float));  // real part
    REQUIRE(a.has_nan());
}

TEST_CASE("pgCol has_nan: detects NaN in pgComplex<float> imag part", "[pgCol_has_nan]") {
    pgCol<pgComplex<float>> a(100);
    for (arma::uword i = 0; i < 100; i++)
        a.at(i) = pgComplex<float>(1.0f, 2.0f);

    // Inject NaN into the imag part of element 99 (second float of interleaved pair)
    float nanVal = make_nan_float();
    float* rawPtr = reinterpret_cast<float*>(a.memptr());
    memcpy(&rawPtr[99 * 2 + 1], &nanVal, sizeof(float));  // imag part
    REQUIRE(a.has_nan());
}

// =========================================================================
// pgMat conj()
// =========================================================================

TEST_CASE("pgMat conj: element-wise conjugation", "[pgMat_conj]") {
    const arma::uword nR = 4, nC = 3;
    pgMat<pgComplex<float>> M(nR, nC);
    for (arma::uword i = 0; i < M.n_elem; i++)
        M.at(i) = pgComplex<float>((float)i, (float)(i + 1));

    pgMat<pgComplex<float>> Mc = conj(M);
    REQUIRE(Mc.n_rows == nR);
    REQUIRE(Mc.n_cols == nC);

    for (arma::uword i = 0; i < M.n_elem; i++) {
        REQUIRE(Mc.at(i).real() == Approx(M.at(i).real()));
        REQUIRE(Mc.at(i).imag() == Approx(-M.at(i).imag()));
    }
}

// =========================================================================
// pgMat reshape constructor
// =========================================================================

TEST_CASE("pgMat reshape constructor from pgCol", "[pgMat_reshape]") {
    const arma::uword N = 24;
    pgCol<float> v(N);
    for (arma::uword i = 0; i < N; i++) v.at(i) = (float)i;

    pgMat<float> M(v, 6, 4);  // 6 rows, 4 cols
    REQUIRE(M.n_rows == 6);
    REQUIRE(M.n_cols == 4);
    REQUIRE(M.n_elem == 24);

    // Column-major: M(r,c) = v[c*nR + r]
    for (arma::uword c = 0; c < 4; c++)
        for (arma::uword r = 0; r < 6; r++)
            REQUIRE(M.at(r, c) == Approx((float)(c * 6 + r)));
}

// =========================================================================
// pgMat arma::Mat<complex> constructor
// =========================================================================

TEST_CASE("pgMat from arma::Mat<complex<float>>", "[pgMat_arma_ctor]") {
    arma::Mat<std::complex<float>> armM(4, 3);
    for (arma::uword c = 0; c < 3; c++)
        for (arma::uword r = 0; r < 4; r++)
            armM(r, c) = std::complex<float>((float)(r + c), (float)(r - (int)c));

    pgMat<pgComplex<float>> M(armM);
    REQUIRE(M.n_rows == 4);
    REQUIRE(M.n_cols == 3);

    for (arma::uword c = 0; c < 3; c++)
        for (arma::uword r = 0; r < 4; r++) {
            REQUIRE(M.at(r, c).real() == Approx(armM(r, c).real()));
            REQUIRE(M.at(r, c).imag() == Approx(armM(r, c).imag()));
        }
}

// =========================================================================
// Real × complex operator%
// =========================================================================

TEST_CASE("operator%(pgCol<float>, pgCol<pgComplex<float>>): real * complex", "[pgCol_rvec_cmul]") {
    const arma::uword N = 256;
    pgCol<float> W(N);
    pgCol<pgComplex<float>> X(N);

    for (arma::uword i = 0; i < N; i++) {
        W.at(i) = (float)(i + 1) * 0.01f;
        X.at(i) = pgComplex<float>((float)i * 0.1f, -(float)i * 0.05f);
    }

    pgCol<pgComplex<float>> result = W % X;
    REQUIRE(result.n_elem == N);
    REQUIRE_FALSE(result.is_view());

    for (arma::uword i = 0; i < N; i++) {
        float w = W.at(i);
        float expectedRe = w * X.at(i).real();
        float expectedIm = w * X.at(i).imag();
        REQUIRE(result.at(i).real() == Approx(expectedRe).margin(1e-5f));
        REQUIRE(result.at(i).imag() == Approx(expectedIm).margin(1e-5f));
    }
}

TEST_CASE("operator%(pgCol<float>, pgCol<pgComplex<float>>): large N for Metal dispatch", "[pgCol_rvec_cmul]") {
    // N > 4096 should trigger Metal dispatch if METAL_COMPUTE is on
    const arma::uword N = 8192;
    pgCol<float> W(N);
    pgCol<pgComplex<float>> X(N);

    for (arma::uword i = 0; i < N; i++) {
        W.at(i) = 2.0f;
        X.at(i) = pgComplex<float>(3.0f, 4.0f);
    }

    pgCol<pgComplex<float>> result = W % X;

    for (arma::uword i = 0; i < N; i++) {
        // 2 * (3 + 4i) = (6 + 8i)
        REQUIRE(result.at(i).real() == Approx(6.0f));
        REQUIRE(result.at(i).imag() == Approx(8.0f));
    }
}

// =========================================================================
// Aligned allocation (Step 4 verification)
// =========================================================================

#ifdef METAL_COMPUTE
TEST_CASE("pgCol: memory is page-aligned (16384)", "[pgCol_aligned]") {
    pgCol<float> a(1024);
    uintptr_t addr = reinterpret_cast<uintptr_t>(a.memptr());
    REQUIRE((addr % 16384) == 0);
}

TEST_CASE("pgMat: memory is page-aligned (16384)", "[pgMat_aligned]") {
    pgMat<float> M(32, 32);
    uintptr_t addr = reinterpret_cast<uintptr_t>(M.memptr());
    REQUIRE((addr % 16384) == 0);
}

TEST_CASE("pgCol<pgComplex<float>>: memory is page-aligned", "[pgCol_aligned]") {
    pgCol<pgComplex<float>> a(2048);
    uintptr_t addr = reinterpret_cast<uintptr_t>(a.memptr());
    REQUIRE((addr % 16384) == 0);
}
#endif

// =========================================================================
// Move semantics with views
// =========================================================================

TEST_CASE("pgCol view: move constructor transfers view status", "[pgCol_view]") {
    const arma::uword N = 32;
    pgCol<float> owner(N);
    owner.ones();

    pgCol<float> v = pgCol<float>::view(owner.memptr(), N);
    float* origPtr = v.memptr();

    pgCol<float> moved(std::move(v));
    REQUIRE(moved.is_view());
    REQUIRE(moved.memptr() == origPtr);
    REQUIRE(moved.n_elem == N);
}

TEST_CASE("pgCol view: move assignment transfers view status", "[pgCol_view]") {
    const arma::uword N = 32;
    pgCol<float> owner(N);
    owner.ones();

    pgCol<float> v = pgCol<float>::view(owner.memptr(), N);
    float* origPtr = v.memptr();

    pgCol<float> target(10);  // owning
    target = std::move(v);
    REQUIRE(target.is_view());
    REQUIRE(target.memptr() == origPtr);
    REQUIRE(target.n_elem == N);
}
