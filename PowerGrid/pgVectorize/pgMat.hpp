// pgMat.hpp

#ifndef POWER_GRID_pgMat_hpp
#define POWER_GRID_pgMat_hpp

#include "../PGIncludes.h"
//#include "pgComplex.hpp"
#include "pgCol.hpp"
#include "pgSubviewCol.hpp"

#ifdef _OPENACC
#include "openacc.h"
#endif

template<typename T>
class pgMat {

private:


bool isInitialized; 
bool isOnGPU; 

// Number of elements in array
//arma::uword n_elem;

public:
T *mem; //Pointer to raw data
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
        #pragma acc enter data create(this) 
    }


pgMat<T>(arma::uword nRows, arma::uword nCols ) :
    isOnGPU(false),
    isInitialized(false),
    mem(NULL),
    n_elem(0),
    n_cols(0),
    n_rows(0) {
        
    #pragma acc enter data create(this) 
    set_size(nRows, nCols);

}

pgMat<T>(arma::Mat<std::complex<T>> &cSCplx) :
    isOnGPU(false),
    isInitialized(false),
    mem(NULL),
    n_elem(0),
    n_cols(0),
    n_rows(0) {
        
    #pragma acc enter data create(this) 
    set_size(cSCplx.n_elem);

    for(arma::uword ii = 0; ii < n_elem; ii++) {
        mem[ii] = cSCplx(ii);
    }
    #pragma acc update device(mem[0:n_elem])

}




// Copy Constructor
pgMat<T>(const pgMat<T>& pgA) :
    isOnGPU(false),
    isInitialized(false),
    mem(NULL),
    n_elem(0),
    n_cols(0),
    n_rows(0)  {
    #pragma acc enter data create(this) 

    set_size(pgA.n_elem);
    size_t bytes = sizeof(T) * pgA.n_elem;
    #ifdef _OPENACC
        isOnGPU = true;
    #endif
    #pragma acc update device(this) 
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

    #pragma acc enter data create(this) 
    isInitialized = true;
    isOnGPU = true;
    access::rw(n_elem) = pgA.n_elem;
    access::rw(n_rows) = pgA.n_rows;
    access::rw(n_cols) = pgA.n_cols;
    mem = pgA.memptr();

    #ifdef _OPENACC
        isOnGPU = true;
    #endif
    #pragma acc update device(this) 
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

    #pragma acc exit data delete(this)

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
    #pragma acc update device(this)

}

void set_size(arma::uword nRows, arma::uword nCols) {
    if (isInitialized) {
        if (isOnGPU) {
            #pragma acc exit data finalize detach(mem) delete(mem[0:n_elem])
            isOnGPU  = false; 
        }
        delete[] mem;
        mem = NULL;
    }
    
    isInitialized = true;
    arma::access::rw(n_cols) = nCols;
    arma::access::rw(n_rows) = nRows;
    arma::access::rw(n_elem) = nRows * nCols;

    mem = new T[n_elem];
    #ifdef _OPENACC
        isOnGPU = true;
    #endif
    #pragma acc update device(this[0:1])
    #pragma acc enter data create(mem[0:n_elem])

} 

void zeros() {
    #pragma acc parallel loop copyin(this [0:1]) present(mem[0:n_elem])
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        mem[ii] = T();
    }
}

void ones() {
    #pragma acc parallel loop copyin(this [0:1]) present(mem[0:n_elem])
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        mem[ii] = T(1.0);
    }
}

// Conversion from pgMat to arma::Mat
arma::Mat<T> getArma() {
    #pragma acc update host(mem[0:n_elem])
    arma::Mat<T> armaT(mem, n_rows, n_cols, true, false);

    return armaT;
}

template<typename T_> 
arma::Mat<T_> 
    getArma(typename std::enable_if<std::is_same<T_,std::complex<float>>::value, std::nullptr_t>::type = nullptr)
{
#pragma acc update host(mem [0:2*inT.n_elem])
    arma::Col<std::complex<float>> armaT(mem, n_elem, true, false);

    return armaT;
}

template<typename T_> 
arma::Mat<T_> 
    getArma(typename std::enable_if<std::is_same<T_,std::complex<double>>::value, std::nullptr_t>::type = nullptr)
{
#pragma acc update host(mem [0:2*inT.n_elem])
    arma::Col<std::complex<double>> armaT(mem, n_elem, true, false);

    return armaT;
}

// Return a column for use
pgSubviewCol<T> col(const arma::uword colIndx) {

    // Calculate some of the prerequisites
    arma::uword n_elemCol = this->n_rows;
    arma::uword uiColHeader = this->n_rows * colIndx;
    T* pMem = this->mem;
    pgSubviewCol<T> pgC(pMem, this->n_elem, n_elemCol, uiColHeader);
    
    return std::move(pgC);
}

// Operators for element manipulation
// We'll assume .at() is for fast, GPU manipulation
#pragma acc routine seq
inline
const T at(const arma::uword d) const {
    return mem[d];
}

#pragma acc routine seq
inline
const T at(const arma::uword rowIdx, const arma::uword colIdx) const {
    return mem[n_cols * colIdx + rowIdx];
}

#pragma acc routine seq
inline
T& at(const arma::uword d) {
    return mem[d];
}

#pragma acc routine seq
inline
T& at(const arma::uword rowIdx, const arma::uword colIdx) {
    return mem[n_cols * colIdx + rowIdx];
}

#pragma acc routine seq
inline
T& operator()(const arma::uword d) {
    //#pragma acc update_host(mem)
    return mem[d];
}

#pragma acc routine seq
inline
const T operator()(const arma::uword d) const {
    //#pragma acc update_host(mem)
    return mem[d];
}

#pragma acc routine seq
inline
T& operator()(const arma::uword rowIdx, const arma::uword colIdx) {
    //#pragma acc update_host(mem)
    return mem[n_cols * colIdx + rowIdx];
}

#pragma acc routine seq
inline
const T operator()(const arma::uword rowIdx, const arma::uword colIdx) const {
    //#pragma acc update_host(mem)
    return mem[n_cols * colIdx + rowIdx];
}

pgMat<T>& operator=(const pgMat<T>& d) {
    #pragma acc parallel loop present(mem[0:n_elem],d)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] = d.at(ii);
    }
    return *this;
}

pgMat<T>& operator=(pgMat<T>&& d) {

    if (isInitialized) {
        #pragma acc exit data delete(mem[0:n_elem])
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
    #pragma acc update device(this)   
    #ifdef _OPENACC
        acc_attach((void **) &mem);
    #endif

    d.reset_mem();

    return *this;
}

pgMat<T>& operator+=(const T& A) {

    #pragma acc parallel loop copyin(this [0:1]) present(mem[0:n_elem])
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] += A;
    }
    return *this;
}

pgMat<T>& operator-=(const T& A) {
    #pragma acc parallel loop present(this [0:1], mem[0:n_elem])
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] -= A;
    }
    return *this;
}

template<typename X>
pgMat<T>& operator%=(const T& A) {
    #pragma acc parallel loop present(this, mem[0:n_elem])
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] *= A;
    }
    return *this;
}

template<typename X>
pgMat<T>& operator/=(const T& A) {
    #pragma acc parallel loop present(this, mem[0:n_elem])
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        mem[ii] /= A;
    }
    return *this;
}

template<typename X>
pgMat<T>& operator+=(const pgMat<X> &pgA) {

    #pragma acc parallel loop present(this, mem[0:n_elem], pgB, pgA)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] += pgA.at(ii);
    }
    return *this;
}

template<typename X>
pgMat<T>& operator-=(const pgMat<X> &pgA) {
    #pragma acc parallel loop present(this, mem[0:n_elem], pgB, pgA)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] -= pgA.at(ii);
    }
    return *this;
}

template<typename X>
pgMat<T>& operator%=(const pgMat<X> &pgA) {
    #pragma acc parallel loop present(this, mem[0:n_elem], pgB, pgA)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] *= pgA.at(ii);
    }
    return *this;
}

template<typename X>
pgMat<T>& operator/=(const pgMat<X> &pgA) {
    #pragma acc parallel loop present(this, mem[0:n_elem], pgB, pgA)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        mem[ii] /= pgA.at(ii);
    }
    return *this;
}

// Operators(pgMat, scalar)
//template<typename T>
pgMat<T> operator+(const T& B) {
    pgMat<T> pgC(n_elem);

    #pragma acc parallel loop present(this, mem[0:n_elem], pgB, pgC)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] + B;
    }
    return std::move(pgC);
}

//template<typename T>
pgMat<T> operator-( const T& B) {
    pgMat<T> pgC(n_elem);

    #pragma acc parallel loop present(this, mem[0:n_elem], pgB, pgC)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] - B;
    }
    return std::move(pgC);
}

//template<typename T>
pgMat<T> operator%(const T& B) {
    pgMat<T> pgC(n_elem);

    #pragma acc parallel loop present(this, mem[0:n_elem], pgB, SpgC)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] * B;
    }
    return std::move(pgC);
}

//template<typename T>
pgMat<T> operator/(const T& B) {
    pgMat<T> pgC(n_elem);

    #pragma acc parallel loop present(this, mem[0:n_elem], pgB, pgC)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] / B;
    }
    return std::move(pgC);
}

// Operators(pgMat, pgMat)
template<typename X>
pgMat<T> operator+(const pgMat<X>& pgB) {
    pgMat<T> pgC(n_rows, n_cols);

    #pragma acc parallel loop present(this, mem[0:n_elem], pgB, pgC)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] + pgB.at(ii);
    }
    return std::move(pgC);
}

template<typename X>
pgMat<T> operator-(const pgMat<X>& pgB) {
    pgMat<T> pgC(n_rows, n_cols);

    #pragma acc parallel loop present(this, mem[0:n_elem], pgB, pgC)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] - pgB.at(ii);
    }
    return std::move(pgC);
}
template<typename X>
pgMat<T> operator%(const pgMat<X>& pgB) {
    pgMat<T> pgC(n_rows, n_cols);

    #pragma acc parallel loop present(this, mem[0:n_elem], pgB, pgC)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] * pgB.at(ii);
    }

    return std::move(pgC);
}

template<typename X>
pgMat<T> operator/(const pgMat<X>& pgB) {
    pgMat<T> pgC(n_rows, n_cols);

    #pragma acc parallel loop present(this, mem[0:n_elem], pgC)
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] / pgB.at(ii);
    }
    return pgC;
}

};

template<typename T>
const pgCol<std::complex<T>> sum(const pgMat<std::complex<T>> &pgA, const arma::uword dim = 0) {
    pgCol<T> sumReal = {};
    pgCol<T> sumImag = {};

    if (dim == 0) { // Colum-wise sums (default)
        sumReal.set_size(pgA.n_cols);
        sumImag.set_size(pgA.n_cols);
        sumReal.zeros();
        sumReal.zeros();

        #pragma acc parallel loop present(pgA, sumImag)
        for(arma::uword jj = 0; jj < pgA.n__cols; jj++) {
            #pragma acc loop seq
            for(arma::uword ii = 0; ii < pgA.n_rows; ii++) {
                sumReal.at(jj) += real(pgA.at(jj * pgA.n_rows + ii ));
            }
        }
        #pragma acc parallel loop present(pgA, sumImag)
        for(arma::uword jj = 0; jj < pgA.n_cols; jj++) {
            #pragma acc loop seq
            for(arma::uword ii = 0; ii < pgA.n_rows; ii++) {
                sumReal.at(jj) += imag(pgA.at(jj * pgA.n_rows + ii ));
            }
        }
    } else if (dim == 1) {
        sumReal.set_size(pgA.n_rows);
        sumImag.set_size(pgA.n_rows);
        sumReal.zeros();
        sumReal.zeros();

        #pragma acc parallel loop present(pgA, sumImag)
        for(arma::uword jj = 0; jj < pgA.n__rows; jj++) {
            #pragma acc loop seq
            for(arma::uword ii = 0; ii < pgA.n_cols; ii++) {
                sumReal.at(jj) += real(pgA.at(jj + pgA.n_rows * ii ));
            }
        }
        #pragma acc parallel loop present(pgA, sumImag)
        for(arma::uword jj = 0; jj < pgA.n_rows; jj++) {
            #pragma acc loop seq
            for(arma::uword ii = 0; ii < pgA.n_cols; ii++) {
                sumReal.at(jj) += imag(pgA.at(jj + pgA.n_rows * ii ));
            }
        }
    } else {
        std::cout << "pgMat::sum Error! Unrecognized dimension: dim = " << dim << std::endl;       

    }
    std::complex<T> J(0,1.0);

    pgCol<std::complex<T>> out(sumReal.n_elem);
    #pragma acc parallel loop present(sum, sumReal, sumImag)
    for(arma::uword jj = 0; jj < sumReal.n_elem; jj++) {
        out.at(jj) = std::complex<T>(sumReal.at(jj),sumImag.at(jj));
    }
    
    return std::move(out);

}

template<typename T>
const pgCol<T> sum(const pgMat<T> &pgA, const arma::uword dim = 0) {
    pgCol<T> sumA;

    if (dim == 0) { // Column-wise sums (default)
        sumA.set_size(pgA.n_cols);
        sumA.zeros();
        #pragma acc parallel loop present(pgA, sumA)
        for(arma::uword jj = 0; jj < pgA.n_cols; jj++) {
            #pragma acc loop seq
            for(arma::uword ii = 0; ii < pgA.n_rows; ii++) {
                sumA.at(jj) += pgA.at(jj * pgA.n_rows + ii );
            }
        }

    } else if (dim == 1) { // Row-wise sums 
        sumA.set_size(pgA.n_rows);
        sumA.zeros();
        #pragma acc parallel loop present(pgA, sumA)
        for(arma::uword jj = 0; jj < pgA.n_rows; jj++) {
            #pragma acc loop seq
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
const T accu(const pgMat<T> &pgA) {
    T sumA = 0;
    #pragma acc parallel loop present(pgA, pgA.mem[0:pgA.n_elem]) reduction(+ : sumA)
    for(arma::uword ii = 0; ii < pgA.n_elem; ii++) {
        sumA += pgA.at(ii);
    }

    return sumA;
}

template<typename T>
const pgCol<T> vectorise(const pgMat<T> &pgA) {
    pgCol<T> vectA(pgA.n_elem);

    #pragma acc parallel loop present(pgA, vectA, pgA.mem[0:pgA.n_elem], vectA.mem[0:vectA.n_elem])
    for(arma::uword ii = 0; ii < pgA.n_elem; ii++) {
        vectA.at(ii) = pgA.at(ii);
    }
    return std::move(vectA);
}

#endif // POWER_GRID_pgMat_hpp
