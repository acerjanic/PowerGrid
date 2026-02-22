#include "catch.hpp"

#ifdef _OPENACC
#include "accel.h"
#endif

#include "../Core/PGIncludes.h"
#include "../Core/pgMat.hpp"
#include <cmath>
#include <iostream>


TEST_CASE("pgMat<float>: operators", "[pgMat<float>") {

    // Setup Prerequsites for test
    arma::uword lengthA = 256;
    arma::uword lengthB = 1000;
    
    pgMat<float> pgA(lengthA,lengthA);
    pgMat<float> pgB(lengthB,lengthB);
    pgMat<float> pgC;
    pgMat<float> pgD;
    
    pgC.set_size(lengthA,lengthA);
    pgD.set_size(lengthB,lengthB);

    pgMat<float> pgCZeros(lengthA,lengthA);
    pgMat<float> pgDZeros(lengthB,lengthB);
    
    SECTION( "Test constructor and n_elems " ) {

        REQUIRE(pgA.n_elem == lengthA * lengthA);
        REQUIRE(pgB.n_elem == lengthB * lengthB);
    } 

    SECTION( "Test constructor and n_rows " ) {

        REQUIRE(pgA.n_rows == lengthA);
        REQUIRE(pgB.n_rows == lengthB);
    } 

    SECTION( "Test constructor and n_cols " ) {

        REQUIRE(pgA.n_cols == lengthA);
        REQUIRE(pgB.n_cols == lengthB);
    } 

    SECTION( "Test set_size() and n_elem" ) {

        REQUIRE(pgC.n_elem == lengthA * lengthA);
        REQUIRE(pgD.n_elem == lengthB * lengthB);
    }

    SECTION( "Test set_size() and n_rows " ) {

        REQUIRE(pgC.n_rows == lengthA);
        REQUIRE(pgD.n_rows == lengthB);
    } 

    SECTION( "Test set_size() and n_cols " ) {

        REQUIRE(pgC.n_cols == lengthA);
        REQUIRE(pgD.n_cols == lengthB);
    } 


    SECTION( "Test zeros and sum" ) {
        pgC.zeros();
        pgD.zeros();
        REQUIRE(sum(sum(pgC)) == 0);
        REQUIRE(sum(sum(pgD)) == 0);
    }

    SECTION( "Test ones and sum()" ) {
        pgC.ones();
        pgD.ones();
        REQUIRE(sum(sum(pgC)) == (float)(lengthA * lengthA));
        REQUIRE(sum(sum(pgD)) == (float)(lengthB * lengthB));
    }

    SECTION( "Test Operator+" ) {
        pgC.ones();
        pgD.ones();
        REQUIRE(sum(sum( pgC + pgC )) == 2 * lengthA * lengthA);
        REQUIRE(sum(sum( pgD + pgD )) == 2 * lengthB * lengthB);
        
        pgC.zeros();
        pgD.zeros();
        REQUIRE(sum(sum( pgC + pgC )) == 0);
        REQUIRE(sum(sum( pgD + pgD )) == 0);
    }

    SECTION( "Test Operator-" ) {
        pgC.ones();
        pgD.ones();
        REQUIRE(sum(sum( pgC - pgC )) == 0);
        REQUIRE(sum(sum( pgD - pgD )) == 0);
    }

    SECTION( "Test Operator%" ) {
        pgC.ones();
        pgD.ones();

        pgCZeros.zeros();
        pgDZeros.zeros();
        REQUIRE(sum(sum( pgC % pgC )) == lengthA * lengthA);
        REQUIRE(sum(sum( pgD % pgD )) == lengthB * lengthB);
        REQUIRE(sum(sum( pgC % pgCZeros)) == 0);
        REQUIRE(sum(sum( pgD % pgDZeros)) == 0);

    }

}

TEST_CASE("pgMat<float>: extended correctness", "[pgMat<float> extended]") {

    // Use a non-square matrix: 4 rows x 3 cols
    // pgMat(A, B) → n_rows=A, n_cols=B (same convention as Armadillo)
    arma::uword nRows = 4;
    arma::uword nCols = 3;
    pgMat<float> m(nRows, nCols);

    SECTION("Dimensions of non-square matrix") {
        REQUIRE(m.n_rows == nRows);
        REQUIRE(m.n_cols == nCols);
        REQUIRE(m.n_elem == nRows * nCols);
    }

    SECTION("2D element access at(row, col) matches linear at(linear)") {
        m.zeros();
        // Fill with distinct values using 2D indexing
        for (arma::uword c = 0; c < nCols; c++) {
            for (arma::uword r = 0; r < nRows; r++) {
                m.at(r, c) = static_cast<float>(r + c * nRows + 1);
            }
        }
        // Verify consistency with linear (column-major) indexing
        for (arma::uword c = 0; c < nCols; c++) {
            for (arma::uword r = 0; r < nRows; r++) {
                REQUIRE(m.at(r, c) == m.at(r + c * nRows));
            }
        }
    }

    SECTION("col() extracts correct column values") {
        m.zeros();
        // Set column 1 to [10, 20, 30, 40]
        for (arma::uword r = 0; r < nRows; r++) {
            m.at(r, 1) = static_cast<float>((r + 1) * 10);
        }
        pgCol<float> c1 = m.col(1);
        REQUIRE(c1.n_elem == nRows);
        for (arma::uword r = 0; r < nRows; r++) {
            REQUIRE(c1.at(r) == Approx(static_cast<float>((r + 1) * 10)));
        }
    }

    SECTION("vectorise() flattens in column-major order") {
        // Use a 3-rows x 2-cols matrix: pgMat(nCols=2, nRows=3)
        pgMat<float> m32(2, 3);
        // Fill with 1..6 using linear index
        for (arma::uword i = 0; i < 6; i++) {
            m32.at(i) = static_cast<float>(i + 1);
        }
        pgCol<float> v = vectorise(m32);
        REQUIRE(v.n_elem == 6);
        // Column-major: elements in same order as linear index
        for (arma::uword i = 0; i < 6; i++) {
            REQUIRE(v.at(i) == Approx(static_cast<float>(i + 1)));
        }
    }

    SECTION("sum(m, 0) gives column sums: length == n_cols, each entry == sum of column") {
        m.ones();
        pgCol<float> colSums = sum(m, 0);
        REQUIRE(colSums.n_elem == nCols);
        for (arma::uword c = 0; c < nCols; c++) {
            REQUIRE(colSums.at(c) == Approx(static_cast<float>(nRows)));
        }
    }

    SECTION("sum(m, 1) gives row sums: length == n_rows, each entry == sum of row") {
        m.ones();
        pgCol<float> rowSums = sum(m, 1);
        REQUIRE(rowSums.n_elem == nRows);
        for (arma::uword r = 0; r < nRows; r++) {
            REQUIRE(rowSums.at(r) == Approx(static_cast<float>(nCols)));
        }
    }

    SECTION("Copy constructor produces independent deep copy") {
        m.ones();
        pgMat<float> mCopy(m);
        REQUIRE(mCopy.n_rows == nRows);
        REQUIRE(mCopy.n_cols == nCols);
        // Modifying copy does not affect original
        mCopy.at(0) = 99.0f;
        REQUIRE(m.at(0) == 1.0f);
    }

}
