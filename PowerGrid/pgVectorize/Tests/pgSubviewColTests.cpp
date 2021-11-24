#include "../../../Support/Catch/catch.hpp"

#ifdef _OPENACC
#include "accel.h"
#endif

#include "../../PGIncludes.h"
#include "../pgMat.hpp"
#include "../pgSubviewCol.hpp"
#include <cmath>
#include <iostream>


TEST_CASE("pgSubviewCol<float>: operators", "pgSubviewCol<float>") {

    // Setup Prerequsites for test
    arma::uword lengthA = 32;
    arma::uword lengthB = 256;
    
    pgMat<float> pgA(lengthA,lengthA);
    pgMat<float> pgB(lengthB,lengthB);

    pgMat<float> pgC(lengthA,lengthA);
    pgMat<float> pgD(lengthB,lengthB);

    pgMat<float> pgAZeros(lengthA,lengthA);
    pgMat<float> pgBZeros(lengthB,lengthB);

    pgMat<float> pgCZeros(lengthA,lengthA);
    pgMat<float> pgDZeros(lengthB,lengthB);

    pgA.ones();
    pgB.ones();

    pgC.zeros();
    pgD.zeros();

    pgCZeros.zeros();
    pgDZeros.zeros();

    for(int ii = 0; ii < lengthA; ii++) {
        for(int jj = 0; jj < lengthA; jj++) {
            
            // Set the values of the matrix to be equal to the element's column index (one indexed)
            pgC(ii,jj) = (jj+1);

        }
    }

    for(int ii = 0; ii < lengthB; ii++) {
        for(int jj = 0; jj < lengthB; jj++) {
            
            // Set the values of the matrix to be equal to the element's row index
            pgD(ii,jj) = (ii+1);

        }
    }

    // Generate subviews we will work with.
    
    SECTION( "Test constructor and n_rows " ) {
        pgSubviewCol<float> pgACol = pgA.col(1);    
        pgSubviewCol<float> pgBCol = pgB.col(1);   
        REQUIRE(pgACol.n_rows == lengthA);
        REQUIRE(pgBCol.n_rows == lengthB);

        pgSubviewCol<float> pgACol2 = pgA.col(lengthA - 1);    
        pgSubviewCol<float> pgBCol2 = pgB.col(lengthA - 1);   
        REQUIRE(pgACol2.n_rows == lengthA);
        REQUIRE(pgBCol2.n_rows == lengthB);
    } 
        
    SECTION( "Test constructor and n_elems " ) {
        pgSubviewCol<float> pgACol = pgA.col(1);    
        pgSubviewCol<float> pgBCol = pgB.col(1);   
        REQUIRE(pgACol.n_elem() == lengthA);
        REQUIRE(pgBCol.n_elem() == lengthB);

        pgSubviewCol<float> pgACol2 = pgA.col(lengthA - 1);    
        pgSubviewCol<float> pgBCol2 = pgB.col(lengthA - 1);   
        REQUIRE(pgACol2.n_elem() == lengthA);
        REQUIRE(pgBCol2.n_elem() == lengthB);
    } 

    SECTION( "Test constructor and n_cols " ) {
        pgSubviewCol<float> pgACol = pgA.col(1);    
        pgSubviewCol<float> pgBCol = pgB.col(1);   
        REQUIRE(pgACol.n_cols == 1);
        REQUIRE(pgBCol.n_cols == 1);
    } 

    SECTION( "Test zeros and sum" ) {
        pgSubviewCol<float> pgACol = pgA.col(1);    
        pgSubviewCol<float> pgBCol = pgB.col(1);   
        REQUIRE(sum(pgACol) == lengthA);
        REQUIRE(sum(pgBCol) == lengthB);
    }

    SECTION( "Test ones and sum()" ) {
        pgSubviewCol<float> pgACol = pgA.col(1);    
        pgSubviewCol<float> pgBCol = pgB.col(1);   
        REQUIRE(sum(pgACol) == lengthA);
        REQUIRE(sum(pgBCol) == lengthB);
    }

    SECTION( "Test Operator+" ) {
        pgA.ones();
        pgB.ones();

        pgSubviewCol<float> pgACol = pgA.col(1);    
        pgSubviewCol<float> pgBCol = pgB.col(1);  
        REQUIRE(sum( pgACol + pgACol ) == 2 * lengthA);
        REQUIRE(sum( pgBCol + pgBCol ) == 2 * lengthB);
        
        pgA.zeros();
        pgB.zeros();

        pgSubviewCol<float> pgACol2 = pgA.col(1);    
        pgSubviewCol<float> pgBCol2 = pgB.col(1);  
        REQUIRE(sum( pgACol2 + pgACol2 ) == 0);
        REQUIRE(sum( pgBCol2 + pgBCol2 ) == 0);
    }

    SECTION( "Test Operator-" ) {
        pgC.ones();
        pgD.ones();
        REQUIRE(accu( pgC - pgC ) == 0);
        REQUIRE(accu( pgD - pgD ) == 0);
    }

    SECTION( "Test Operator%" ) {
        pgC.ones();
        pgD.ones();

        pgCZeros.zeros();
        pgDZeros.zeros();
        REQUIRE(accu( pgC % pgC ) == lengthA * lengthA);
        REQUIRE(accu( pgD % pgD ) == lengthB * lengthB);
        REQUIRE(accu( pgC % pgCZeros) == 0);
        REQUIRE(accu( pgD % pgDZeros) == 0);
        
    }

}
