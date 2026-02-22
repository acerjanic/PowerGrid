#include "catch.hpp"

#ifdef _OPENACC
#include "openacc.h"
#include "accel.h"
#endif

#include "../Core/PGIncludes.h"
#include "../Core/pgCol.hpp"
#include "../Core/pgComplex.hpp"
#include <cmath>
#include <iostream>

TEST_CASE("pgCol<pgComplex<float>>: operators", "[pgCol<pgComplex>>]") {

    // Setup Prerequsites for test
    arma::uword lengthA = 100000;
    arma::uword lengthB = 256 * 256;
    
    pgCol<pgComplex<float>> pgCA(lengthA);
    pgCol<pgComplex<float>> pgCB(lengthB);
    pgCol<pgComplex<float>> pgCC;
    pgCol<pgComplex<float>> pgCD;
    pgCA.zeros();
    pgCB.zeros();

    pgCC.set_size(lengthA);
    pgCD.set_size(lengthB);

    pgCol<pgComplex<float>> pgCZeros(lengthA);
    pgCol<pgComplex<float>> pgDZeros(lengthB);
    pgCZeros.zeros();
    pgDZeros.zeros();

    SECTION( "Test constructor and n_elems " ) {
        REQUIRE(pgCA.n_elem == lengthA);
        REQUIRE(pgCB.n_elem == lengthB);
    }



    SECTION( "Test set_size()" ) {
        REQUIRE(pgCC.n_elem == lengthA);
        REQUIRE(pgCD.n_elem == lengthB);
    }

    SECTION( "Test zeros and sum" ) {
        pgCC.zeros();
        pgCD.zeros();
        REQUIRE(abs(sum(pgCC)) == 0);
        REQUIRE(abs(sum(pgCD)) == 0);
    }

    SECTION( "Test  ones and sum()" ) {
        pgCC.ones();
        pgCD.ones();
        REQUIRE(abs(sum(pgCC)) == (float)lengthA);
        REQUIRE(abs(sum(pgCD)) == (float)lengthB);
    }

    SECTION( "Test Operator+" ) {
        pgCC.ones();
        pgCD.ones();
        REQUIRE(abs(sum( pgCC + pgCC )) == 2 * lengthA);
        REQUIRE(abs(sum( pgCD + pgCD )) == 2 * lengthB);
        
        pgCC.zeros();
        pgCD.zeros();
        REQUIRE(abs(sum( pgCC + pgCC )) == 0);
        REQUIRE(abs(sum( pgCD + pgCD )) == 0);
    }

    SECTION( "Test Operator-" ) {
        pgCC.ones();
        pgCD.ones();
        REQUIRE(abs(sum( pgCC - pgCC )) == 0);
        REQUIRE(abs(sum( pgCD - pgCD )) == 0);
    }

    SECTION( "Test Operator*" ) {
        pgCC.ones();
        pgCD.ones();
        REQUIRE(abs(sum( pgCC % pgCC )) == lengthA);
        REQUIRE(abs(sum( pgCD % pgCD )) == lengthB);
        REQUIRE(abs(sum( pgCC % pgCZeros)) == 0);
        REQUIRE(abs(sum( pgCD % pgDZeros)) == 0);

    }

}

TEST_CASE("pgCol<pgComplex<float>>: extended correctness", "[pgCol<pgComplex> extended]") {

    arma::uword N = 6;
    pgCol<pgComplex<float>> v(N);

    float re_ = 1.5f;
    float im_ = -2.0f;

    SECTION("Element access at() write and read round-trip") {
        v.zeros();
        for (arma::uword i = 0; i < N; i++) {
            v.at(i) = pgComplex<float>(static_cast<float>(i), -static_cast<float>(i));
        }
        for (arma::uword i = 0; i < N; i++) {
            REQUIRE(v.at(i).real() == Approx(static_cast<float>(i)));
            REQUIRE(v.at(i).imag() == Approx(-static_cast<float>(i)));
        }
    }

    SECTION("Conjugate of element: real preserved, imag negated") {
        pgComplex<float> z(re_, im_);
        v.zeros();
        v.at(0) = z;
        pgComplex<float> zConj = conj(v.at(0));
        REQUIRE(zConj.real() == Approx(re_));
        REQUIRE(zConj.imag() == Approx(-im_));
    }

    SECTION("z * conj(z) has zero imaginary part (equals norm)") {
        pgComplex<float> z(re_, im_);
        v.zeros();
        v.at(0) = z;
        pgComplex<float> product = v.at(0) * conj(v.at(0));
        REQUIRE(product.imag() == Approx(0.0f).margin(1e-5f));
        REQUIRE(product.real() == Approx(norm(z)));
    }

    SECTION("sum() real and imag parts match manual accumulation") {
        for (arma::uword i = 0; i < N; i++) {
            v.at(i) = pgComplex<float>(static_cast<float>(i + 1), static_cast<float>(i) * 0.5f);
        }
        // Manual accumulation
        float expectedReal = 0.0f, expectedImag = 0.0f;
        for (arma::uword i = 0; i < N; i++) {
            expectedReal += v.at(i).real();
            expectedImag += v.at(i).imag();
        }
        pgComplex<float> s = sum(v);
        REQUIRE(s.real() == Approx(expectedReal));
        REQUIRE(s.imag() == Approx(expectedImag));
    }

}

