// pgMat.hpp

#ifndef POWER_GRID_pgMat_hpp
#define POWER_GRID_pgMat_hpp

#include "PGIncludes.h"
#include "pgComplex.hpp"
#include "pgCol.hpp"

#ifdef _OPENACC
#include "openacc.h"
#endif

template<typename T>
class pgMat {

private:

T *mem; //Pointer to raw data
bool isInitialized;
bool isOnGPU;

// Number of elements in array
//arma::uword n_elem;

public:

const arma::uword n_elem;
const arma::uword n_rows;
const arma::uword n_cols;
// Constructors

pgMat<T>() :
    isOnGPU(false),
    isInitialized(false),
    mem(NULL),
    n_elem(0),
    n_cols(0),
    n_rows(0) {
        #ifdef _OPENACC
        #pragma acc enter data create(this)
        #endif
    }


pgMat<T>(arma::uword nRows, arma::uword nCols ) :
    isOnGPU(false),
    isInitialized(false),
    mem(NULL),
    n_elem(0),
    n_cols(0),
    n_rows(0) {

    #ifdef _OPENACC
    #pragma acc enter data create(this)
    #endif
    set_size(nCols, nRows);

}

pgMat<T>(arma::Mat<std::complex<T>> &cSCplx) :
    isOnGPU(false),
    isInitialized(false),
    mem(NULL),
    n_elem(0),
    n_cols(0),
    n_rows(0) {

    #ifdef _OPENACC
    #pragma acc enter data create(this)
    #endif
    set_size(cSCplx.n_cols, cSCplx.n_rows);

    for(arma::uword ii = 0; ii < n_elem; ii++) {
        mem[ii] = cSCplx(ii);
    }
    #ifdef _OPENACC
    #pragma acc update device(mem[0:n_elem])
    #endif

}




// Copy Constructor
pgMat<T>(const pgMat<T>& pgA) :
    isOnGPU(false),
    isInitialized(false),
    mem(NULL),
    n_elem(0),
    n_cols(0),
    n_rows(0)  {
    #ifdef _OPENACC
    #pragma acc enter data create(this)
    #endif

    set_size(pgA.n_cols, pgA.n_rows);
    size_t bytes = sizeof(T) * pgA.n_elem;
    #ifdef _OPENACC
        isOnGPU = true;
    #endif
    #ifdef _OPENACC
    #pragma acc update device(this)
    #endif
    memcpy(mem, pgA.memptr(), bytes);
    #ifdef _OPENACC
        acc_memcpy(acc_deviceptr(mem), acc_deviceptr(pgA.memptr()), bytes);
    #endif

}

// Move Constructor
pgMat<T>(pgMat<T>&& pgA) :
    isOnGPU(false),
    isInitialized(false),
    mem(NULL),
    n_elem(0),
    n_cols(0),
    n_rows(0) {

    #ifdef _OPENACC
    #pragma acc enter data create(this)
    #endif
    isInitialized = true;
    isOnGPU = true;
    access::rw(n_elem) = pgA.n_elem;
    access::rw(n_rows) = pgA.n_rows;
    access::rw(n_cols) = pgA.n_cols;
    mem = pgA.memptr();

    #ifdef _OPENACC
        isOnGPU = true;
    #endif
    #ifdef _OPENACC
    #pragma acc update device(this)
    #endif
    #ifdef _OPENACC
        acc_attach((void **) &mem);
    #endif

    pgA.reset_mem();

}



// Default Destructor
~pgMat<T>() {

    #ifdef _OPENACC
        if( acc_deviceptr(mem) != NULL) {
            acc_delete((void *)mem, sizeof(T) * n_elem);
        }
    #endif
    if (mem != NULL) {
        delete[] mem;
    }

    #ifdef _OPENACC
    #pragma acc exit data delete(this)
    #endif

}

T* memptr() const {
    return mem;
}

void reset_mem() {
    #ifdef _OPENACC
        acc_detach((void **) &mem);
    #endif
    mem = NULL;
    isInitialized = false;
    access::rw(n_elem) = 0;
    access::rw(n_rows) = 0;
    access::rw(n_cols) = 0;
    isOnGPU = false;
    #ifdef _OPENACC
    #pragma acc update device(this)
    #endif

}

void set_size(arma::uword nCols, arma::uword nRows) {
    if (isInitialized) {
        if (isOnGPU) {
            #ifdef _OPENACC
            #pragma acc exit data finalize detach(mem) delete(mem[0:n_elem])
            #endif
            isOnGPU  = false;
        }
        delete[] mem;
        mem = NULL;
    }

    isInitialized = true;
    arma::access::rw(n_cols) = nCols;
    arma::access::rw(n_rows) = nRows;
    arma::access::rw(n_elem) = nCols * nRows;

    mem = new T[n_elem];
    #ifdef _OPENACC
        isOnGPU = true;
    #endif
    #ifdef _OPENACC
    #pragma acc update device(this)
    #pragma acc enter data create(mem[0:n_elem])
    #endif

}

void zeros() {
    #ifdef _OPENACC
    #pragma acc parallel loop present(mem[0:n_elem])
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        mem[ii] = T();
    }
}

void ones() {
    #ifdef _OPENACC
    #pragma acc parallel loop present(mem[0:n_elem])
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        mem[ii] = T(1.0);
    }
}

// Conversion from pgMat to arma::Mat
arma::Mat<T> getArma() {
    #ifdef _OPENACC
    #pragma acc update host(mem[0:n_elem])
    #endif
    arma::Mat<T> armaT(mem, n_rows, n_cols, true, false);

    return armaT;
}

// Return a column for use
pgCol<T> col(const arma::uword colIndx) {

    pgCol<T> pgC(n_rows);

    #ifdef _OPENACC
    #pragma acc parallel loop present(pgC, mem[0:n_elem])
    #endif
    for(arma::uword ii = 0; ii < n_rows; ii++ ) {
        pgC.at(ii) = mem[n_rows*colIndx + ii];
    }
    return std::move(pgC);
}

// Operators for element manipulation
// We'll assume .at() is for fast, GPU manipulation
inline
const T at(const arma::uword d) const {
    return mem[d];
}

inline
const T at(const arma::uword rowIdx, const arma::uword colIdx) const {
    return mem[n_rows * colIdx + rowIdx];
}

inline
T& at(const arma::uword d) {
    return mem[d];
}

inline
T& at(const arma::uword rowIdx, const arma::uword colIdx) {
    return mem[n_rows * colIdx + rowIdx];
}

inline
T& operator()(const arma::uword d) {
    //#pragma acc update_host(mem)
    return mem[d];
}

inline
const T operator()(const arma::uword d) const {
    //#pragma acc update_host(mem)
    return mem[d];
}

inline
T& operator()(const arma::uword rowIdx, const arma::uword colIdx) {
    //#pragma acc update_host(mem)
    return mem[n_rows * colIdx + rowIdx];
}

inline
const T operator()(const arma::uword rowIdx, const arma::uword colIdx) const {
    //#pragma acc update_host(mem)
    return mem[n_rows * colIdx + rowIdx];
}

pgMat<T>& operator=(const pgMat<T>& d) {
    #ifdef _OPENACC
    #pragma acc parallel loop present(mem[0:n_elem],d)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] = d.at(ii);
    }
    return *this;
}

pgMat<T>& operator=(pgMat<T>&& d) {

    if (isInitialized) {
        #ifdef _OPENACC
        #pragma acc exit data delete(mem[0:n_elem])
        #endif
        delete[] mem;
        access::rw(n_elem) = 0;
        access::rw(n_cols) = 0;
        access::rw(n_rows) = 0;
        isInitialized = false;
        isOnGPU = false;
    }
    access::rw(n_elem) = d.n_elem;
    access::rw(n_cols) = d.n_cols;
    access::rw(n_rows) = d.n_rows;
    isInitialized = true;
    isOnGPU = true;
    mem = d.memptr();
    #ifdef _OPENACC
    #pragma acc update device(this)
    #endif
    #ifdef _OPENACC
        acc_attach((void **) &mem);
    #endif

    d.reset_mem();

    return *this;
}

pgMat<T>& operator+=(const T& A) {

    #ifdef _OPENACC
    #pragma acc parallel loop present(mem[0:n_elem])
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] += A;
    }
    return *this;
}

pgMat<T>& operator-=(const T& A) {
    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem])
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] -= A;
    }
    return *this;
}

pgMat<T>& operator%=(const T& A) {
    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem])
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] *= A;
    }
    return *this;
}

pgMat<T>& operator/=(const T& A) {
    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem])
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        mem[ii] /= A;
    }
    return *this;
}

template<typename X>
pgMat<T>& operator+=(const pgMat<X> &pgA) {

    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem], pgA)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] += pgA.at(ii);
    }
    return *this;
}

template<typename X>
pgMat<T>& operator-=(const pgMat<X> &pgA) {
    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem], pgA)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] -= pgA.at(ii);
    }
    return *this;
}

template<typename X>
pgMat<T>& operator%=(const pgMat<X> &pgA) {
    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem], pgA)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] *= pgA.at(ii);
    }
    return *this;
}

template<typename X>
pgMat<T>& operator/=(const pgMat<X> &pgA) {
    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem], pgA)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        mem[ii] /= pgA.at(ii);
    }
    return *this;
}

// Operators(pgMat, scalar)
pgMat<T> operator+(const T& B) const {
    pgMat<T> pgC(n_rows, n_cols);

    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem], pgC)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] + B;
    }
    return std::move(pgC);
}

pgMat<T> operator-( const T& B) const {
    pgMat<T> pgC(n_rows, n_cols);

    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem], pgC)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] - B;
    }
    return std::move(pgC);
}

pgMat<T> operator%(const T& B) const {
    pgMat<T> pgC(n_rows, n_cols);

    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem], pgC)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] * B;
    }
    return std::move(pgC);
}

pgMat<T> operator/(const T& B) const {
    pgMat<T> pgC(n_rows, n_cols);

    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem], pgC)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] / B;
    }
    return std::move(pgC);
}

// Operators(pgMat, pgMat)
template<typename X>
pgMat<T> operator+(const pgMat<X>& pgB) const {
    pgMat<T> pgC(n_rows, n_cols);

    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem], pgB, pgC)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] + pgB.at(ii);
    }
    return std::move(pgC);
}

template<typename X>
pgMat<T> operator-(const pgMat<X>& pgB) const {
    pgMat<T> pgC(n_rows, n_cols);

    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem], pgB, pgC)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] - pgB.at(ii);
    }
    return std::move(pgC);
}
template<typename X>
pgMat<T> operator%(const pgMat<X>& pgB) const {
    pgMat<T> pgC(n_rows, n_cols);

    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem], pgB, pgC)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] * pgB.at(ii);
    }

    return std::move(pgC);
}

template<typename X>
pgMat<T> operator/(const pgMat<X>& pgB) const {
    pgMat<T> pgC(n_rows, n_cols);

    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem], pgC)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] / pgB.at(ii);
    }
    return pgC;
}

};

template<typename T>
const pgCol<pgComplex<T>> sum(const pgMat<pgComplex<T>> &pgA, const arma::uword dim = 0) {
    pgCol<T> sumReal;
    pgCol<T> sumImag;

    if (dim == 0) { // Column-wise sums (default)
        sumReal.set_size(pgA.n_cols);
        sumImag.set_size(pgA.n_cols);
        sumReal.zeros();
        sumImag.zeros();

        #ifdef _OPENACC
        #pragma acc parallel loop present(pgA, sumReal)
        #endif
        for(arma::uword jj = 0; jj < pgA.n_cols; jj++) {
            #ifdef _OPENACC
            #pragma acc loop seq
            #endif
            for(arma::uword ii = 0; ii < pgA.n_rows; ii++) {
                sumReal.at(jj) += real(pgA.at(jj * pgA.n_rows + ii ));
            }
        }
        #ifdef _OPENACC
        #pragma acc parallel loop present(pgA, sumImag)
        #endif
        for(arma::uword jj = 0; jj < pgA.n_cols; jj++) {
            #ifdef _OPENACC
            #pragma acc loop seq
            #endif
            for(arma::uword ii = 0; ii < pgA.n_rows; ii++) {
                sumImag.at(jj) += imag(pgA.at(jj * pgA.n_rows + ii ));
            }
        }
    } else if (dim == 1) {
        sumReal.set_size(pgA.n_rows);
        sumImag.set_size(pgA.n_rows);
        sumReal.zeros();
        sumImag.zeros();

        #ifdef _OPENACC
        #pragma acc parallel loop present(pgA, sumReal)
        #endif
        for(arma::uword jj = 0; jj < pgA.n_rows; jj++) {
            #ifdef _OPENACC
            #pragma acc loop seq
            #endif
            for(arma::uword ii = 0; ii < pgA.n_cols; ii++) {
                sumReal.at(jj) += real(pgA.at(jj + pgA.n_rows * ii ));
            }
        }
        #ifdef _OPENACC
        #pragma acc parallel loop present(pgA, sumImag)
        #endif
        for(arma::uword jj = 0; jj < pgA.n_rows; jj++) {
            #ifdef _OPENACC
            #pragma acc loop seq
            #endif
            for(arma::uword ii = 0; ii < pgA.n_cols; ii++) {
                sumImag.at(jj) += imag(pgA.at(jj + pgA.n_rows * ii ));
            }
        }
    } else {
        std::cout << "pgMat::sum Error! Unrecognized dimension: dim = " << dim << std::endl;

    }
    pgComplex<T> J(0,1.0);

    pgCol<pgComplex<T>> out(sumReal.n_elem);
    #ifdef _OPENACC
    #pragma acc parallel loop present(out, sumReal, sumImag)
    #endif
    for(arma::uword jj = 0; jj < sumReal.n_elem; jj++) {
        out.at(jj) = pgComplex<T>(sumReal.at(jj),sumImag.at(jj));
    }

    return std::move(out);

}

template<typename T>
const pgCol<T> sum(const pgMat<T> &pgA, const arma::uword dim = 0) {
    pgCol<T> sumA;

    if (dim == 0) { // Column-wise sums (default)
        sumA.set_size(pgA.n_cols);
        sumA.zeros();
        #ifdef _OPENACC
        #pragma acc parallel loop present(pgA, sumA)
        #endif
        for(arma::uword jj = 0; jj < pgA.n_cols; jj++) {
            #ifdef _OPENACC
            #pragma acc loop seq
            #endif
            for(arma::uword ii = 0; ii < pgA.n_rows; ii++) {
                sumA.at(jj) += pgA.at(jj * pgA.n_rows + ii );
            }
        }

    } else if (dim == 1) { // Row-wise sums
        sumA.set_size(pgA.n_rows);
        sumA.zeros();
        #ifdef _OPENACC
        #pragma acc parallel loop present(pgA, sumA)
        #endif
        for(arma::uword jj = 0; jj < pgA.n_rows; jj++) {
            #ifdef _OPENACC
            #pragma acc loop seq
            #endif
            for(arma::uword ii = 0; ii < pgA.n_cols; ii++) {
                sumA.at(jj) += pgA.at(jj + pgA.n_rows * ii );
            }
        }
    } else {
        std::cout << "pgMat::sum Error! Unrecognized dimension: dim = " << dim << std::endl;
    }
    return std::move(sumA);
}


template<typename T>
const pgCol<T> vectorise(const pgMat<T> &pgA) {
    pgCol<T> vectA(pgA.n_elem);

    #ifdef _OPENACC
    #pragma acc parallel loop present(pgA, vectA)
    #endif
    for(arma::uword ii = 0; ii < pgA.n_elem; ii++) {
        vectA.at(ii) = pgA.at(ii);
    }
    return std::move(vectA);
}

#endif // POWER_GRID_pgMat_hpp
