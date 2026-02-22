#include "catch.hpp"

#ifdef _OPENACC
#include "accel.h"
#endif

#include "../Core/PGIncludes.h"
#include "../Core/pgCol.hpp"
#include <cmath>
#include <iostream>

TEST_CASE("pgCol<float>: operators", "[pgCol<float>]") {

    // Setup Prerequsites for test
    arma::uword lengthA = 100000;
    arma::uword lengthB = 256 * 256;
    
    pgCol<float> pgA(lengthA);
    pgCol<float> pgB(lengthB);
    pgCol<float> pgC;
    pgCol<float> pgD;
    
    pgC.set_size(lengthA);
    pgD.set_size(lengthB);

    pgCol<float> pgCZeros(lengthA);
    pgCol<float> pgDZeros(lengthB);
    
    SECTION( "Test constructor and n_elems " ) {

        REQUIRE(pgA.n_elem == lengthA);
        REQUIRE(pgB.n_elem == lengthB);
    } 



    SECTION( "Test set_size()" ) {

        REQUIRE(pgC.n_elem == lengthA);
        REQUIRE(pgD.n_elem == lengthB);
    }

    SECTION( "Test zeros and sum" ) {
        pgC.zeros();
        pgD.zeros();
        REQUIRE(sum(pgC) == 0);
        REQUIRE(sum(pgD) == 0);
    }

    SECTION( "Test ones and sum()" ) {
        pgC.ones();
        pgD.ones();
        REQUIRE(sum(pgC) == (float)lengthA);
        REQUIRE(sum(pgD) == (float)lengthB);
    }

    SECTION( "Test Operator+" ) {
        pgC.ones();
        pgD.ones();
        REQUIRE(sum( pgC + pgC ) == 2 * lengthA);
        REQUIRE(sum( pgD + pgD ) == 2 * lengthB);
        
        pgC.zeros();
        pgD.zeros();
        REQUIRE(sum( pgC + pgC ) == 0);
        REQUIRE(sum( pgD + pgD ) == 0);
    }

    SECTION( "Test Operator-" ) {
        pgC.ones();
        pgD.ones();
        REQUIRE(sum( pgC - pgC ) == 0);
        REQUIRE(sum( pgD - pgD ) == 0);
    }

    SECTION( "Test Operator%" ) {
        pgC.ones();
        pgD.ones();

        pgCZeros.zeros();
        pgDZeros.zeros();
        REQUIRE(sum( pgC % pgC ) == lengthA);
        REQUIRE(sum( pgD % pgD ) == lengthB);
        REQUIRE(sum( pgC % pgCZeros) == 0);
        REQUIRE(sum( pgD % pgDZeros) == 0);

    }

}

TEST_CASE("pgCol<float>: extended correctness", "[pgCol<float> extended]") {

    arma::uword N = 8;
    pgCol<float> v(N);

    SECTION("Element access at() write and read round-trip") {
        v.zeros();
        for (arma::uword i = 0; i < N; i++) {
            v.at(i) = static_cast<float>(i + 1);
        }
        for (arma::uword i = 0; i < N; i++) {
            REQUIRE(v.at(i) == static_cast<float>(i + 1));
        }
    }

    SECTION("operator() matches at() for element access") {
        v.zeros();
        v.at(3) = 42.0f;
        REQUIRE(v(3) == 42.0f);
    }

    SECTION("Copy constructor produces independent deep copy") {
        v.ones();
        pgCol<float> vCopy(v);
        REQUIRE(vCopy.n_elem == N);
        // Modifying copy does not affect original
        vCopy.at(0) = 99.0f;
        REQUIRE(v.at(0) == 1.0f);
    }

    SECTION("Move constructor transfers ownership") {
        v.ones();
        v.at(2) = 5.0f;
        pgCol<float> vMoved(std::move(v));
        REQUIRE(vMoved.n_elem == N);
        REQUIRE(vMoved.at(2) == 5.0f);
        // Source memory pointer is now null
        REQUIRE(v.memptr() == nullptr);
    }

    SECTION("Scalar += increases all elements by the scalar") {
        v.ones();
        v += 2.0f;
        REQUIRE(sum(v) == Approx(static_cast<float>(N) * 3.0f));
        for (arma::uword i = 0; i < N; i++) {
            REQUIRE(v.at(i) == Approx(3.0f));
        }
    }

    SECTION("Scalar -= decreases all elements by the scalar") {
        v.ones();
        v -= 0.5f;
        REQUIRE(sum(v) == Approx(static_cast<float>(N) * 0.5f));
    }

    SECTION("sum() of ones equals n_elem") {
        v.ones();
        REQUIRE(sum(v) == Approx(static_cast<float>(N)));
    }

    SECTION("sum() accuracy on known values") {
        // Load v with 1, 2, 3, ..., N
        for (arma::uword i = 0; i < N; i++) {
            v.at(i) = static_cast<float>(i + 1);
        }
        float expected = static_cast<float>(N * (N + 1)) / 2.0f;
        REQUIRE(sum(v) == Approx(expected));
    }

    SECTION("at() write/read preserves floating-point values accurately") {
        for (arma::uword i = 0; i < N; i++) {
            v.at(i) = static_cast<float>(i) * 0.1f;
        }
        for (arma::uword i = 0; i < N; i++) {
            REQUIRE(v.at(i) == Approx(static_cast<float>(i) * 0.1f));
        }
    }

}
