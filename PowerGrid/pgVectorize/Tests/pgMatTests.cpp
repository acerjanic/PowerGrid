#include "../../../Support/Catch/catch.hpp"

#ifdef _OPENACC
#include "accel.h"
#endif

#include "../../PGIncludes.h"
#include "../pgMat.hpp"
#include "../pgSubviewCol.hpp"
#include <cmath>
#include <iostream>


TEST_CASE("pgMat<float>: operators", "pgMat<float>") {

    // Setup Prerequsites for test
    arma::uword lengthA = 64; 
    arma::uword lengthB = 128;
    
    pgMat<float> pgA(lengthA,lengthA);
    pgMat<float> pgB(lengthB,lengthB);
    pgMat<float> pgC;
    pgMat<float> pgD;
    pgMat<float> pgE;
    pgMat<float> pgF;
    
    pgC.set_size(lengthA,lengthA);
    pgD.set_size(lengthB,lengthB);

    pgMat<float> pgCZeros(lengthA,lengthA);
    pgMat<float> pgDZeros(lengthB,lengthB);

    pgE.set_size(lengthA,lengthB);
    pgF.set_size(lengthB,lengthA);

    pgMat<float> pgEZeros(lengthA,lengthB);
    pgMat<float> pgFZeros(lengthB,lengthA);
    
    SECTION( "Square Matrices: Test constructor and n_elems " ) {

        REQUIRE(pgA.n_elem == lengthA * lengthA);
        REQUIRE(pgB.n_elem == lengthB * lengthB);
    } 

    SECTION( "Square Matrices: Test constructor and n_rows " ) {

        REQUIRE(pgA.n_rows == lengthA);
        REQUIRE(pgB.n_rows == lengthB);
    } 

    SECTION( "Square Matrices: Test constructor and n_cols " ) {

        REQUIRE(pgA.n_cols == lengthA);
        REQUIRE(pgB.n_cols == lengthB);
    } 

    SECTION( "Square Matrices: Test set_size() and n_elem" ) {

        REQUIRE(pgC.n_elem == lengthA * lengthA);
        REQUIRE(pgD.n_elem == lengthB * lengthB);
    }

    SECTION( "Square Matrices: Test set_size() and n_rows " ) {

        REQUIRE(pgC.n_rows == lengthA);
        REQUIRE(pgD.n_rows == lengthB);
    } 

    SECTION( "Square Matrices: Test set_size() and n_cols " ) {

        REQUIRE(pgC.n_cols == lengthA);
        REQUIRE(pgD.n_cols == lengthB);
    } 


    SECTION( "Square Matrices: Test zeros and sum" ) {
        pgC.zeros();
        pgD.zeros();
        REQUIRE(sum(sum(pgC)) == 0);
        REQUIRE(sum(sum(pgD)) == 0);
    }

    SECTION( "Square Matrices: Test ones and sum()" ) {
        pgC.ones();
        pgD.ones();
        REQUIRE(sum(sum(pgC)) == (float)(lengthA * lengthA));
        REQUIRE(sum(sum(pgD)) == (float)(lengthB * lengthB));
    }

    SECTION( "Square Matrices: Test Operator+" ) {
        pgC.ones();
        pgD.ones();
        REQUIRE(accu( pgC + pgC ) == 2 * lengthA * lengthA);
        REQUIRE(accu( pgD + pgD ) == 2 * lengthB * lengthB);
        
        pgC.zeros();
        pgD.zeros();
        REQUIRE(accu(pgC + pgC ) == 0);
        REQUIRE(accu( pgD + pgD ) == 0);
    }

    SECTION( "Square Matrices: Test Operator-" ) {
        pgC.ones();
        pgD.ones();
        REQUIRE(accu( pgC - pgC ) == 0);
        REQUIRE(accu( pgD - pgD ) == 0);
    }

    SECTION( "Square Matrices: Test Operator%" ) {
        pgC.ones();
        pgD.ones();

        pgCZeros.zeros();
        pgDZeros.zeros();
        REQUIRE(accu( pgC % pgC ) == lengthA * lengthA);
        REQUIRE(accu( pgD % pgD ) == lengthB * lengthB);
        REQUIRE(accu( pgC % pgCZeros) == 0);
        REQUIRE(accu( pgD % pgDZeros) == 0);
        
    }

    SECTION( "Square Matrices: Test () " ) {
        pgC.zeros();
        pgD.zeros();

        pgCZeros.zeros();
        pgDZeros.zeros();

        #pragma acc parallel loop present(pgC, pgC.mem[0:pgC.n_elem])
        for(arma::uword ii = 0; ii < lengthA; ii++) {
            for(arma::uword jj = 0; jj < lengthA; jj++) {
                
                // Set the values of the matrix to be equal to the element's column index (one indexed)
                pgC(ii,jj) = (jj+1);

            }
        }

        #pragma acc parallel loop present(pgD, pgD.mem[0:pgC.n_elem])
        for(arma::uword ii = 0; ii < lengthB; ii++) {
            for(arma::uword jj = 0; jj < lengthB; jj++) {
                
                // Set the values of the matrix to be equal to the element's row index
                pgD(ii,jj) = (jj+1);

            }
        }

        REQUIRE(accu( pgC ) == lengthA * ((lengthA * (lengthA + 1)) / 2));
        REQUIRE(accu( pgD ) == (((float)lengthB * ((float)lengthB + 1)) / 2.0f)  * (float)lengthB );
        REQUIRE(accu( pgC % pgCZeros) == 0);
        REQUIRE(accu( pgD % pgDZeros) == 0);
        
    }

}
