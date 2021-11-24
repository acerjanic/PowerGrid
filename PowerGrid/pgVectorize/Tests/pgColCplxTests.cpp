#include "../../../Support/Catch/catch.hpp"

#ifdef _OPENACC
#include "openacc.h"
#include "accel.h"
#endif

#include "../../PGIncludes.h"
#include "../pgCol.hpp"
//#include "../pgVectorize/std::complex.hpp"
#include <cmath>
#include <iostream>

TEST_CASE("pgCol<std::complex<float>>: operators", "[pgCol<std::complex>>]") {

    // Setup Prerequsites for test
    arma::uword lengthA = 100;
    arma::uword lengthB = 256;
    
    pgCol<std::complex<float>> pgCA(lengthA);
    pgCol<std::complex<float>> pgCB(lengthB);
    pgCol<std::complex<float>> pgCC;
    pgCol<std::complex<float>> pgCD;
    pgCA.zeros();
    pgCB.zeros();

    pgCC.set_size(lengthA);
    pgCD.set_size(lengthB);

    pgCol<std::complex<float>> pgCZeros(lengthA);
    pgCol<std::complex<float>> pgDZeros(lengthB);
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
        REQUIRE(abs(sum(pgCC)) == (float)0);
        REQUIRE(abs(sum(pgCD)) == (float)0);
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

