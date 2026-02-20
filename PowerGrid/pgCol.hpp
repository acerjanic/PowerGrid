// pgCol.hpp

#ifndef POWER_GRID_pgCol_hpp
#define POWER_GRID_pgCol_hpp

#include "PGIncludes.h"
#include "pgComplex.hpp"
#include <type_traits>

#ifdef _OPENACC
#include "openacc.h"
#include "accel.h"
#endif

#ifdef METAL_COMPUTE
#include "Metal/MetalVectorOps_dispatch.hpp"
#endif


template<typename T>
class pgCol {

private:

T *mem; //Pointer to raw data
bool isInitialized;
bool isOnGPU;
bool isCopy;
// Number of elements in array
//arma::uword n_elem;

public:

const arma::uword n_elem;

// Constructors

pgCol<T>() :
    isOnGPU(false),
    isInitialized(false),
    mem(NULL),
    isCopy(false),
    n_elem(0) {
        #ifdef _OPENACC
        #pragma acc enter data copyin(this)
        #endif
    }


pgCol<T>(arma::uword length) :
    isOnGPU(false),
    isInitialized(false),
    mem(NULL),
    isCopy(false),
    n_elem(0) {
    #ifdef _OPENACC
    #pragma acc enter data create(this)
    #endif
    set_size(length);
}


// Constructor from arma::Col<T> for non-complex types (float, double, etc.)
template <typename U = T,
          typename std::enable_if<!std::is_same<U, pgComplex<float>>::value &&
                                  !std::is_same<U, pgComplex<double>>::value,
                                  int>::type = 0>
pgCol(arma::Col<T> &cSCplx) :
    isOnGPU(false),
    isInitialized(false),
    mem(NULL),
    isCopy(false),
    n_elem(0) {
    std::cout << "Entering pgCol<T> copy from Armadillo constructor" << std::endl;

    #ifdef _OPENACC
    #pragma acc enter data create(this)
    #endif
    set_size(cSCplx.n_elem);
    memcpy(this->mem, cSCplx.memptr(), sizeof(T) * cSCplx.n_elem);

    #ifdef _OPENACC
    #pragma acc update device(mem[0:n_elem])
    #endif

}

// Constructor from arma::Col<complex<float>> for pgComplex<float>
template <typename U = T,
          typename std::enable_if<std::is_same<U, pgComplex<float>>::value,
                                  int>::type = 0>
pgCol(arma::Col<std::complex<float>> &cSCplx) :
    isOnGPU(false),
    isInitialized(false),
    mem(NULL),
    isCopy(false),
    n_elem(0) {
    std::cout << "Entering pgCol<T> copy from Armadillo cplx constructor" << std::endl;

    #ifdef _OPENACC
    #pragma acc enter data create(this)
    #endif
    set_size(cSCplx.n_elem);
    memcpy(this->mem, reinterpret_cast<float *>(cSCplx.memptr()), sizeof(T) * cSCplx.n_elem);

    #ifdef _OPENACC
    #pragma acc update device(mem[0:n_elem])
    #endif

}

// Constructor from arma::Col<complex<double>> for pgComplex<double>
template <typename U = T,
          typename std::enable_if<std::is_same<U, pgComplex<double>>::value,
                                  int>::type = 0>
pgCol(arma::Col<std::complex<double>> &cSCplx) :
    isOnGPU(false),
    isInitialized(false),
    mem(NULL),
    isCopy(false),
    n_elem(0) {
    std::cout << "Entering pgCol<T> copy from Armadillo cplx constructor" << std::endl;

    #ifdef _OPENACC
    #pragma acc enter data create(this)
    #endif
    set_size(cSCplx.n_elem);
    memcpy(this->mem, reinterpret_cast<double *>(cSCplx.memptr()), sizeof(T) * cSCplx.n_elem);

    #ifdef _OPENACC
    #pragma acc update device(mem[0:n_elem])
    #endif

}



// Copy Constructor
pgCol<T>(const pgCol<T>& pgA) :
    isOnGPU(false),
    isInitialized(false),
    mem(NULL),
    isCopy(true),
    n_elem(0) {
    #ifdef _OPENACC
    #pragma acc enter data create(this)
    #endif
    set_size(pgA.n_elem);
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
pgCol<T>(pgCol<T>&& pgA) :
    isOnGPU(false),
    isInitialized(false),
    mem(NULL),
    isCopy(false),
    n_elem(0) {
    #ifdef _OPENACC
    #pragma acc enter data create(this)
    #endif
    access::rw(n_elem) = pgA.n_elem;
    mem = pgA.memptr();
    isInitialized = true;
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
~pgCol<T>() {

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
    isOnGPU = false;
    #ifdef _OPENACC
    #pragma acc update device(this)
    #endif

}

void set_size(arma::uword length) {
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
    arma::access::rw(n_elem) = length;
    //mem = (T*)malloc(sizeof(T) * n_elem);
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
        mem[ii] = T(0.0);
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

// Conversion from pgCol to arma::Col — non-complex types
template <typename U = T,
          typename std::enable_if<!std::is_same<U, pgComplex<float>>::value &&
                                  !std::is_same<U, pgComplex<double>>::value,
                                  int>::type = 0>
arma::Col<T> getArma() {
    #ifdef _OPENACC
    #pragma acc update host(mem[0:n_elem])
    #endif
    arma::Col<T> armaT(mem, n_elem, true, false);

    return armaT;
}

// Conversion from pgCol<pgComplex<double>> to arma::Col<complex<double>>
template <typename U = T,
          typename std::enable_if<std::is_same<U, pgComplex<double>>::value,
                                  int>::type = 0>
arma::Col<std::complex<double>> getArma() {
    #ifdef _OPENACC
    #pragma acc update host(mem[0:n_elem])
    #endif
    arma::Col<std::complex<double>> armaT(reinterpret_cast<std::complex<double> *>(mem), n_elem, true, false);

    return armaT;
}

// Conversion from pgCol<pgComplex<float>> to arma::Col<complex<float>>
template <typename U = T,
          typename std::enable_if<std::is_same<U, pgComplex<float>>::value,
                                  int>::type = 0>
arma::Col<std::complex<float>> getArma() {
    #ifdef _OPENACC
    #pragma acc update host(mem[0:n_elem])
    #endif
    arma::Col<std::complex<float>> armaT(reinterpret_cast<std::complex<float> *>(mem), n_elem, true, false);

    return armaT;
}

// Operators for element manipulation
// We'll assume .at() is for fast, GPU manipulation
inline
const T at(const arma::uword d) const {
    return mem[d];
}

inline
T& at(const arma::uword d) {
    return mem[d];
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

pgCol<T>& operator=(const pgCol<T>& d) {
    #ifdef _OPENACC
    #pragma acc parallel loop present(mem[0:n_elem],d)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] = d.at(ii);
    }
    return *this;
}

pgCol<T>& operator=(pgCol<T>&& d) {
    //size_t bytes = sizeof(T) * d.n_elem;
    if (isInitialized) {
        #ifdef _OPENACC
        #pragma acc exit data delete(mem[0:n_elem])
        #endif
        delete[] mem;
        access::rw(n_elem) = 0;
        isInitialized = false;
        isOnGPU = false;
    }
    access::rw(n_elem) = d.n_elem;
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

pgCol<T>& operator+=(const T& A) {
    #ifdef METAL_COMPUTE
    if constexpr (std::is_same<T, float>::value) {
        if (pg_metal::try_metal_add_scalar(mem, A, mem, n_elem))
            return *this;
    }
    #endif
    #ifdef _OPENACC
    #pragma acc parallel loop present(mem[0:n_elem])
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] += A;
    }
    return *this;
}

pgCol<T>& operator-=(const T& A) {
    #ifdef METAL_COMPUTE
    if constexpr (std::is_same<T, float>::value) {
        if (pg_metal::try_metal_sub_scalar(mem, A, mem, n_elem))
            return *this;
    }
    #endif
    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem])
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] -= A;
    }
    return *this;
}

pgCol<T>& operator%=(const T& A) {
    #ifdef METAL_COMPUTE
    if constexpr (std::is_same<T, float>::value) {
        if (pg_metal::try_metal_mul_scalar(mem, A, mem, n_elem))
            return *this;
    } else if constexpr (std::is_same<T, pgComplex<float>>::value) {
        if (pg_metal::try_metal_mul_scalar(mem, A, mem, n_elem))
            return *this;
    }
    #endif
    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem])
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] *= A;
    }
    return *this;
}

pgCol<T>& operator/=(const T& A) {
    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem])
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        mem[ii] /= A;
    }
    return *this;
}

template<typename X>
pgCol<T>& operator+=(const pgCol<X> &pgA) {
    #ifdef METAL_COMPUTE
    if constexpr (pg_metal::is_metal_type<T>::value && std::is_same<T, X>::value) {
        if (pg_metal::try_metal_add(mem, pgA.memptr(), mem, n_elem))
            return *this;
    }
    #endif
    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem], pgA)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] += pgA.at(ii);
    }
    return *this;
}

template<typename X>
pgCol<T>& operator-=(const pgCol<X> &pgA) {
    #ifdef METAL_COMPUTE
    if constexpr (pg_metal::is_metal_type<T>::value && std::is_same<T, X>::value) {
        if (pg_metal::try_metal_sub(mem, pgA.memptr(), mem, n_elem))
            return *this;
    }
    #endif
    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem], pgA)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] -= pgA.at(ii);
    }
    return *this;
}

template<typename X>
pgCol<T>& operator%=(const pgCol<X> &pgA) {
    #ifdef METAL_COMPUTE
    if constexpr (pg_metal::is_metal_type<T>::value && std::is_same<T, X>::value) {
        if (pg_metal::try_metal_mul(mem, pgA.memptr(), mem, n_elem))
            return *this;
    }
    #endif
    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem], pgA)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        this->mem[ii] *= pgA.at(ii);
    }
    return *this;
}

template<typename X>
pgCol<T>& operator/=(const pgCol<X> &pgA) {
    #ifdef METAL_COMPUTE
    if constexpr (pg_metal::is_metal_type<T>::value && std::is_same<T, X>::value) {
        if (pg_metal::try_metal_div(mem, pgA.memptr(), mem, n_elem))
            return *this;
    }
    #endif
    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem], pgA)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        mem[ii] /= pgA.at(ii);
    }
    return *this;
}

// Operators(pgCol, scalar)
const pgCol<T> operator+(const T& B) const {
    pgCol<T> pgC(n_elem);
    #ifdef METAL_COMPUTE
    if constexpr (std::is_same<T, float>::value) {
        if (pg_metal::try_metal_add_scalar(mem, B, pgC.memptr(), n_elem))
            return std::move(pgC);
    }
    #endif
    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem], pgC)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] + B;
    }
    return std::move(pgC);
}

const pgCol<T> operator-( const T& B) const {
    pgCol<T> pgC(n_elem);
    #ifdef METAL_COMPUTE
    if constexpr (std::is_same<T, float>::value) {
        if (pg_metal::try_metal_sub_scalar(mem, B, pgC.memptr(), n_elem))
            return std::move(pgC);
    }
    #endif
    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem], pgC)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] - B;
    }
    return std::move(pgC);
}

const pgCol<T> operator%(const T& B) const {
    pgCol<T> pgC(n_elem);
    #ifdef METAL_COMPUTE
    if constexpr (std::is_same<T, float>::value) {
        if (pg_metal::try_metal_mul_scalar(mem, B, pgC.memptr(), n_elem))
            return std::move(pgC);
    } else if constexpr (std::is_same<T, pgComplex<float>>::value) {
        if (pg_metal::try_metal_mul_scalar(mem, B, pgC.memptr(), n_elem))
            return std::move(pgC);
    }
    #endif
    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem], pgC)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] * B;
    }
    return std::move(pgC);
}

const pgCol<T> operator/(const T& B) const {
    pgCol<T> pgC(n_elem);
    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem], pgC)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] / B;
    }
    return std::move(pgC);
}

// Operators(pgCol, pgCol)
template<typename X>
const pgCol<T> operator+(const pgCol<X>& pgB) const {
    pgCol<T> pgC(n_elem);
    #ifdef METAL_COMPUTE
    if constexpr (pg_metal::is_metal_type<T>::value && std::is_same<T, X>::value) {
        if (pg_metal::try_metal_add(mem, pgB.memptr(), pgC.memptr(), n_elem))
            return std::move(pgC);
    }
    #endif
    #ifdef _OPENACC
    #pragma acc parallel loop present( mem[0:n_elem], pgB, pgC)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] + pgB.at(ii);
    }
    return std::move(pgC);
}

template<typename X>
const pgCol<T> operator-(const pgCol<X>& pgB) const {
    pgCol<T> pgC(n_elem);
    #ifdef METAL_COMPUTE
    if constexpr (pg_metal::is_metal_type<T>::value && std::is_same<T, X>::value) {
        if (pg_metal::try_metal_sub(mem, pgB.memptr(), pgC.memptr(), n_elem))
            return std::move(pgC);
    }
    #endif
    #ifdef _OPENACC
    #pragma acc parallel loop present( mem[0:n_elem], pgB, pgC)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] - pgB.at(ii);
    }
    return std::move(pgC);
}
template<typename X>
const pgCol<T> operator%(const pgCol<X>& pgB) const {
    pgCol<T> pgC(n_elem);
    #ifdef METAL_COMPUTE
    if constexpr (pg_metal::is_metal_type<T>::value && std::is_same<T, X>::value) {
        if (pg_metal::try_metal_mul(mem, pgB.memptr(), pgC.memptr(), n_elem))
            return std::move(pgC);
    }
    #endif
    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem], pgB, pgC)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] * pgB.at(ii);
    }

    return std::move(pgC);
}

template<typename X>
const pgCol<T> operator/(const pgCol<X>& pgB) const {
    pgCol<T> pgC(n_elem);
    #ifdef METAL_COMPUTE
    if constexpr (pg_metal::is_metal_type<T>::value && std::is_same<T, X>::value) {
        if (pg_metal::try_metal_div(mem, pgB.memptr(), pgC.memptr(), n_elem))
            return std::move(pgC);
    }
    #endif
    #ifdef _OPENACC
    #pragma acc parallel loop present(this, mem[0:n_elem], pgC)
    #endif
    for(arma::uword ii = 0; ii < n_elem; ii++) {
        pgC.at(ii) = mem[ii] / pgB.at(ii);
    }
    return std::move(pgC);
}

};

template<typename T>
const pgComplex<T> sum(const pgCol<pgComplex<T>> &pgA) {
    T sumReal = {};
    T sumImag = {};

    #ifdef _OPENACC
    #pragma acc parallel loop present(pgA) reduction(+:sumReal)
    #endif
    for(arma::uword ii = 0; ii < pgA.n_elem; ii++) {
        sumReal += real(pgA.at(ii));
    }

    #ifdef _OPENACC
    #pragma acc parallel loop present(pgA) reduction(+:sumImag)
    #endif
    for(arma::uword ii = 0; ii < pgA.n_elem; ii++) {
        sumImag += imag(pgA.at(ii));
    }

    return pgComplex<T>(sumReal, sumImag);

}

template<typename T>
const T sum(const pgCol<T> &pgA) {
    #ifdef METAL_COMPUTE
    if constexpr (std::is_same<T, float>::value) {
        T result;
        if (pg_metal::try_metal_sum(pgA.memptr(), &result, pgA.n_elem))
            return result;
    }
    #endif
    T sumA = {};

    #ifdef _OPENACC
    #pragma acc parallel loop present(pgA) reduction(+:sumA)
    #endif
    for(arma::uword ii = 0; ii < pgA.n_elem; ii++) {
        sumA += pgA.at(ii);
    }
    return sumA;
}

// =========================================================================
// Free functions needed by the PCG solver and operator pipelines
// =========================================================================

/// Complex dot product: sum(conj(A[i]) * B[i]).
template<typename T>
pgComplex<T> cdot(const pgCol<pgComplex<T>>& A, const pgCol<pgComplex<T>>& B) {
    #ifdef METAL_COMPUTE
    if constexpr (std::is_same<T, float>::value) {
        pgComplex<T> result;
        if (pg_metal::try_metal_cdot(A.memptr(), B.memptr(), &result, A.n_elem))
            return result;
    }
    #endif
    T re = T(0), im = T(0);
    for (arma::uword ii = 0; ii < A.n_elem; ii++) {
        pgComplex<T> ca = conj(A.at(ii));
        pgComplex<T> prod = ca * B.at(ii);
        re += prod.real();
        im += prod.imag();
    }
    return pgComplex<T>(re, im);
}

/// L2 norm of a complex vector: sqrt(sum(|A[i]|^2)).
template<typename T>
T norm(const pgCol<pgComplex<T>>& A) {
    #ifdef METAL_COMPUTE
    if constexpr (std::is_same<T, float>::value) {
        T n2;
        if (pg_metal::try_metal_norm2sq(A.memptr(), &n2, A.n_elem))
            return std::sqrt(n2);
    }
    #endif
    T acc = T(0);
    for (arma::uword ii = 0; ii < A.n_elem; ii++) {
        acc += A.at(ii).real() * A.at(ii).real()
             + A.at(ii).imag() * A.at(ii).imag();
    }
    return std::sqrt(acc);
}

/// L2 norm of a real vector.
template<typename T>
T norm(const pgCol<T>& A) {
    #ifdef METAL_COMPUTE
    if constexpr (std::is_same<T, float>::value) {
        // norm^2 = sum(A[i]^2).  Use metal_vec_sum on A.*A via scratch.
        // For now, fall through to CPU — real norm is less critical.
    }
    #endif
    T acc = T(0);
    for (arma::uword ii = 0; ii < A.n_elem; ii++) {
        acc += A.at(ii) * A.at(ii);
    }
    return std::sqrt(acc);
}

/// Element-wise complex conjugate.
template<typename T>
pgCol<pgComplex<T>> conj(const pgCol<pgComplex<T>>& A) {
    pgCol<pgComplex<T>> out(A.n_elem);
    for (arma::uword ii = 0; ii < A.n_elem; ii++) {
        out.at(ii) = conj(A.at(ii));
    }
    return out;
}

/// Element-wise absolute value (magnitude) of complex vector.
template<typename T>
pgCol<T> abs(const pgCol<pgComplex<T>>& A) {
    pgCol<T> out(A.n_elem);
    for (arma::uword ii = 0; ii < A.n_elem; ii++) {
        out.at(ii) = abs(A.at(ii));
    }
    return out;
}

/// Extract real parts of complex vector.
template<typename T>
pgCol<T> real(const pgCol<pgComplex<T>>& A) {
    pgCol<T> out(A.n_elem);
    for (arma::uword ii = 0; ii < A.n_elem; ii++) {
        out.at(ii) = A.at(ii).real();
    }
    return out;
}

/// Extract imaginary parts of complex vector.
template<typename T>
pgCol<T> imag(const pgCol<pgComplex<T>>& A) {
    pgCol<T> out(A.n_elem);
    for (arma::uword ii = 0; ii < A.n_elem; ii++) {
        out.at(ii) = A.at(ii).imag();
    }
    return out;
}

#endif //POWER_GRID_pgCol_hpp
