// pgCol.hpp

#ifndef POWER_GRID_pgCol_hpp
#define POWER_GRID_pgCol_hpp

//#include "../PGIncludes.h"
//#include "pgComplex.hpp"
#include <type_traits>
#include <armadillo>

#ifdef _OPENACC
#include "accel.h"
#include "openacc.h"
#include <cuda_runtime_api.h>
#include <cuda.h>
#endif

// Handy template aliases to fill in gaps in C++11
template< bool Condition, typename T = void >
using enable_if_t = typename std::enable_if<Condition, T>::type;


template <typename T>
class pgCol {

private:

    bool isInitialized;
    bool isOnGPU;
    bool isCopy;
    // Number of elements in array
    //arma::uword n_elem;

public:
    const arma::uword n_elem;
    const arma::uword n_rows;
    T* mem; //Pointer to raw data
    // Constructors

    pgCol<T>()
        : isOnGPU(false)
        , isInitialized(false)
        , mem(NULL)
        , isCopy(false)
        , n_elem(0)
        , n_rows(0) 
        {
#pragma acc enter data copyin(this)
        }

        pgCol<T>(arma::uword length)
        : isOnGPU(false)
        , isInitialized(false)
        , mem(NULL)
        , isCopy(false)
        , n_elem(0)
        , n_rows(0) 
    {
#pragma acc enter data create(this)
        set_size(length);
    }

    pgCol(const arma::Col<T>& cSCplx)
        : isOnGPU(false)
        , isInitialized(false)
        , mem(NULL)
        , isCopy(false)
        , n_elem(0)
        , n_rows(0) 
    {
        #pragma acc enter data copyin(this[0:1])

        this->set_size(cSCplx.n_elem);
        memcpy(this->mem, cSCplx.memptr(), sizeof(T) * cSCplx.n_elem);
        #pragma acc update device(mem[0:cSCplx.n_elem])

    }

    pgCol(const T* pMemory, arma::uword n_elements)
        : isOnGPU(false)
        , isInitialized(false)
        , mem(NULL)
        , isCopy(false)
        , n_elem(0)
        , n_rows(0) 
    {
        #pragma acc enter data copyin(this[0:1])
        this->set_size(n_elements);
        if(acc_is_present((void *)pMemory, sizeof(T) * n_elements)) {
            #pragma acc host_data use_device(pMemory, mem) 
            { 
                cudaMemcpy(this->mem, pMemory, sizeof(T) * n_elements, cudaMemcpyDeviceToDevice); 
            } 
        } else {
            memcpy(this->mem, pMemory, sizeof(T) * n_elements);
            #pragma acc update device(mem[0:n_elements])
        }


    }
/*
    template<typename T_ = T>
    pgCol(arma::Col<T_> &cSCplx,
          typename std::enable_if<std::is_same< T_, std::complex<float>>::value, std::nullptr_t>::type = nullptr)
        : isOnGPU(false)
        , isInitialized(false)
        , mem(NULL)
        , isCopy(false)
        , n_elem(0)
    {

#pragma acc enter data create(this[0:1])
        set_size(cSCplx.n_elem);
        memcpy(this->mem, reinterpret_cast<float*>(cSCplx.memptr()), 2 * sizeof(float) * cSCplx.n_elem);

#pragma acc update device(this[0:1], mem [0:n_elem])
    }

    template<typename T_ = T>
    pgCol(arma::Col<T_>& cSCplx,
          typename std::enable_if<std::is_same<T_,std::complex<double>>::value, std::nullptr_t>::type = nullptr)
        : isOnGPU(false)
        , isInitialized(false)
        , mem(NULL)
        , isCopy(false)
        , n_elem(0)
    {
        std::cout << "Entering pgCol<T> copy from Armadillo cplx constructor" << std::endl;

#pragma acc enter data create(this)
        set_size(cSCplx.n_elem);
        memcpy(this->mem, reinterpret_cast<double*>(cSCplx.memptr()), 2 * sizeof(double) * cSCplx.n_elem);

#pragma acc update device(this[0:1], mem [0:n_elem])
    }
    */


    // Copy Constructor
    pgCol<T>(const pgCol<T>& pgA)
        : isOnGPU(false)
        , isInitialized(false)
        , mem(NULL)
        , isCopy(true)
        , n_elem(0)
        , n_rows(0) 
    {
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
    pgCol<T>(pgCol<T>&& pgA)
        : isOnGPU(false)
        , isInitialized(false)
        , mem(NULL)
        , isCopy(false)
        , n_elem(0)
        , n_rows(0) 
    {
#pragma acc enter data create(this)
        access::rw(n_elem) = pgA.n_elem;
        mem = pgA.memptr();
        isInitialized = true;
#ifdef _OPENACC
        isOnGPU = true;
#endif

#pragma acc update device(this[0:1])
#ifdef _OPENACC
        acc_attach((void**)&mem);
#endif

        pgA.reset_mem();
    }

    // Default Destructor
    ~pgCol<T>()
    {

        #ifdef _OPENACC
            if (acc_deviceptr(mem) != NULL) {
                acc_delete((void*)mem, sizeof(T) * n_elem);
            }
        #endif
        if (mem != NULL) {
            delete[] mem;
        }

        #pragma acc exit data delete (this[0:1])
    }

    T* memptr() const
    {
        return mem;
    }

    T* memptrGPU() const
    {
        #pragma acc host_data use_device(mem)
        {
            //T *memVal = mem;
            return mem;
        }
        
    }

    void reset_mem()
    {
#ifdef _OPENACC
        acc_detach((void**)&mem);
#endif
        mem = NULL;
        isInitialized = false;
        access::rw(n_elem) = 0;
        isOnGPU = false;
#pragma acc update device(this) 
    }

    void set_size(arma::uword length)
    {
        #ifdef _OPENACC
            if (acc_deviceptr(mem) != NULL) {
                acc_delete((void*)mem, sizeof(T) * n_elem);
            }
        #endif
        if (mem != NULL) {
            delete[] mem;
        }
        
        isOnGPU = false;
    
        mem = NULL;
        
        
        arma::access::rw(n_elem) = length;
        arma::access::rw(n_rows) = length;
        mem = new T[n_elem];
        #ifdef _OPENACC
            isOnGPU = true;
        #endif
        #pragma acc update device(this [0:1])
        #pragma acc enter data create(mem [0:n_elem])
        isInitialized = true;
    }

    void zeros()
    {
        #pragma acc parallel loop copyin(this [0:1]) present(mem [0:n_elem])
        for (arma::uword ii = 0; ii < n_elem; ii++) {
            mem[ii] = T(0.0);
        }
    }

    void ones()
    {
        #pragma acc parallel loop copyin(this [0:1]) present(mem [0:n_elem])
        for (arma::uword ii = 0; ii < n_elem; ii++) {
            mem[ii] = T(1.0);
        }
    }

    void ones(arma::uword n_rows)
    {
        set_size(n_rows, 1);

        #pragma acc parallel loop copyin(this [0:1]) present(mem [0:n_elem])
        for (arma::uword ii = 0; ii < n_elem; ii++) {
            mem[ii] = T(1.0);
        }
    }
    
    arma::Col<T> getArma()
    {
        //std::cout << "Memory Address of mem = " << std::hex << mem << std::dec << std::endl;
        //std::cout << "GPU Memory Address of mem = " << std::hex << acc_deviceptr(mem) << std::dec << std::endl;
        #pragma acc wait
        #pragma acc update self(this->mem [0:this->n_elem])
        arma::Col<T> armaT(this->mem, this->n_elem);

        return armaT;
    }
/*
template<typename T_> 
arma::Col<T_> 
    getArma(typename std::enable_if<std::is_same<T_,std::complex<float>>::value, std::nullptr_t>::type = nullptr)
{
    #pragma acc update host(mem [0:2*inT.n_elem])
    arma::Col<std::complex<float>> armaT(mem, n_elem, true, false);

    return armaT;
}

template<typename T_> 
arma::Col<T_> 
    getArma(typename std::enable_if<std::is_same<T_,std::complex<double>>::value, std::nullptr_t>::type = nullptr)
{
    #pragma acc update host(mem [0:2*inT.n_elem])
    arma::Col<std::complex<double>> armaT(mem, n_elem, true, false);

    return armaT;
}
*/
    // Operators for element manipulation
    // We'll assume .at() is for fast, GPU manipulation
    #pragma acc routine seq
    inline const T at(const arma::uword d) const
    {
        return mem[d];
    }

    #pragma acc routine seq
    inline T& at(const arma::uword d)
    {
        return mem[d];
    }

    #pragma acc routine seq
    inline T& operator()(const arma::uword d)
    {
        //#pragma acc update_host(mem)
        return mem[d];
    }

    #pragma acc routine seq
    inline const T operator()(const arma::uword d) const
    {
        //#pragma acc update_host(mem)
        return mem[d];
    }

    pgCol<T>& operator=(const pgCol<T>& d)
    {
        #pragma acc parallel loop present(mem [0:n_elem], d)
        for (arma::uword ii = 0; ii < n_elem; ii++) {
            this->mem[ii] = d.at(ii);
        }
        return *this;
    }

    pgCol<T>& operator=(pgCol<T>&& d)
    {
        //size_t bytes = sizeof(T) * d.n_elem;
        if (isInitialized) {
#pragma acc exit data delete (mem [0:n_elem])
            delete[] mem;
            access::rw(n_elem) = 0;
            isInitialized = false;
            isOnGPU = false;
        }
        access::rw(n_elem) = d.n_elem;
        isInitialized = true;
        isOnGPU = true;
        mem = d.memptr();
#pragma acc update device(this)
#ifdef _OPENACC
        acc_attach((void**)&mem);
#endif

        d.reset_mem();

        return *this;
    }

    pgCol<T>& operator+=(const T& A)
    {

#pragma acc parallel loop present(mem [0:n_elem])
        for (arma::uword ii = 0; ii < n_elem; ii++) {
            this->mem[ii] += A;
        }
        return *this;
    }

    pgCol<T>& operator-=(const T& A)
    {
#pragma acc parallel loop present(this, mem [0:n_elem])
        for (arma::uword ii = 0; ii < n_elem; ii++) {
            this->mem[ii] -= A;
        }
        return *this;
    }

    template <typename X>
    pgCol<T>& operator%=(const T& A)
    {
#pragma acc parallel loop present(this, mem [0:n_elem])
        for (arma::uword ii = 0; ii < n_elem; ii++) {
            this->mem[ii] *= A;
        }
        return *this;
    }

    template <typename X>
    pgCol<T>& operator/=(const T& A)
    {
#pragma acc parallel loop present(this, mem [0:n_elem])
        for (arma::uword ii = 0; ii < n_elem; ii++) {
            mem[ii] /= A;
        }
        return *this;
    }

    template <typename X>
    pgCol<T>& operator+=(const pgCol<X>& pgA)
    {

#pragma acc parallel loop present(this, mem [0:n_elem], pgA)
        for (arma::uword ii = 0; ii < n_elem; ii++) {
            this->mem[ii] += pgA.at(ii);
        }
        return *this;
    }

    template <typename X>
    pgCol<T>& operator-=(const pgCol<X>& pgA)
    {
#pragma acc parallel loop present(this, mem [0:n_elem], pgA)
        for (arma::uword ii = 0; ii < n_elem; ii++) {
            this->mem[ii] -= pgA.at(ii);
        }
        return *this;
    }

    template <typename X>
    pgCol<T>& operator%=(const pgCol<X>& pgA)
    {
#pragma acc parallel loop present(this, mem [0:n_elem], pgA)
        for (arma::uword ii = 0; ii < n_elem; ii++) {
            this->mem[ii] *= pgA.at(ii);
        }
        return *this;
    }

    template <typename X>
    pgCol<T>& operator/=(const pgCol<X>& pgA)
    {
#pragma acc parallel loop present(this, mem [0:n_elem], pgA)
        for (arma::uword ii = 0; ii < n_elem; ii++) {
            mem[ii] /= pgA.at(ii);
        }
        return *this;
    }

    // Operators(pgCol, scalar)
    //template <typename T>
    const pgCol<T> operator+(const T& B) const
    {
        pgCol<T> pgC(n_elem);

#pragma acc parallel loop present(this, mem [0:n_elem], pgC)
        for (arma::uword ii = 0; ii < n_elem; ii++) {
            pgC.at(ii) = mem[ii] + B;
        }
        return std::move(pgC);
    }

    const T min() const {
        T minVal = 0;
        #pragma acc parallel loop present(this[0:1], mem[0:n_elem]) reduction(min: minVal)
        for(int ii = 0; ii < n_elem; ii++) {
            minVal = std::min(minVal, this->mem[ii]);
        }

        return minVal;
    }

    const T max() const {
        T maxVal = 0;
        #pragma acc parallel loop present(this[0:1], mem[0:n_elem]) reduction(max: maxVal)
        for(int ii = 0; ii < n_elem; ii++) {
            maxVal = std::max(maxVal, this->mem[ii]);
        }
        
        return maxVal;
    }

    //template <typename T>
    const pgCol<T> operator-(const T& B) const
    {
        pgCol<T> pgC(n_elem);

#pragma acc parallel loop present(this, mem [0:n_elem], pgC)
        for (arma::uword ii = 0; ii < n_elem; ii++) {
            pgC.at(ii) = mem[ii] - B;
        }
        return std::move(pgC);
    }

    //template <typename T>
    const pgCol<T> operator%(const T& B) const
    {
        pgCol<T> pgC(n_elem);

#pragma acc parallel loop present(this, mem [0:n_elem], SpgC)
        for (arma::uword ii = 0; ii < n_elem; ii++) {
            pgC.at(ii) = mem[ii] * B;
        }
        return std::move(pgC);
    }

    //template <typename T>
    const pgCol<T> operator/(const T& B) const
    {
        pgCol<T> pgC(n_elem);

#pragma acc parallel loop present(this, mem [0:n_elem], pgC)
        for (arma::uword ii = 0; ii < n_elem; ii++) {
            pgC.at(ii) = mem[ii] / B;
        }
        return std::move(pgC);
    }

    // Operators(pgCol, pgCol)
    template <typename X>
    const pgCol<T> operator+(const pgCol<X>& pgB) const
    {
        pgCol<T> pgC(n_elem);

#pragma acc parallel loop copyin(this [0:1]) present(mem [0:n_elem], pgB, pgC)
        for (arma::uword ii = 0; ii < n_elem; ii++) {
            pgC.at(ii) = mem[ii] + pgB.at(ii);
        }
        return std::move(pgC);
    }

    template <typename X>
    const pgCol<T> operator-(const pgCol<X>& pgB) const
    {
        pgCol<T> pgC(n_elem);

#pragma acc parallel loop copyin(this [0:1]) present(mem [0:n_elem], pgB, pgC)
        for (arma::uword ii = 0; ii < n_elem; ii++) {
            pgC.at(ii) = mem[ii] - pgB.at(ii);
        }
        return std::move(pgC);
    }
    template <typename X>
    const pgCol<T> operator%(const pgCol<X>& pgB) const
    {
        pgCol<T> pgC(n_elem);

#pragma acc parallel loop present(this, mem [0:n_elem], pgB, pgC)
        for (arma::uword ii = 0; ii < n_elem; ii++) {
            pgC.at(ii) = mem[ii] * pgB.at(ii);
        }

        return std::move(pgC);
    }

    tcd emplate <typename X>
    const pgCol<T> operator/(const pgCol<X>& pgB) const
    {
        pgCol<T> pgC(n_elem);

#pragma acc parallel loop present(this, mem [0:n_elem], pgC)
        for (arma::uword ii = 0; ii < n_elem; ii++) {
            pgC.at(ii) = mem[ii] / pgB.at(ii);
        }
        return std::move(pgC);
    }
};

template <typename T>
const std::complex<T> sum(const pgCol<std::complex<T>>& pgA)
{
    T sumReal = 0;
    T sumImag = 0;

#pragma acc parallel loop present(pgA, pgA.mem[0:pgA.n_elem]) reduction(+ \
                                                 : sumReal)
    for (arma::uword ii = 0; ii < pgA.n_elem; ii++) {
        sumReal += real(pgA.at(ii));
    }

#pragma acc parallel loop present(pgA, pgA.mem[0:pgA.n_elem]) reduction(+ \
                                                 : sumImag)
    for (arma::uword ii = 0; ii < pgA.n_elem; ii++) {
        sumImag += imag(pgA.at(ii));
    }

    return std::complex<T>(sumReal, sumImag);
}

template <typename T>
const T sum(const pgCol<T>& pgA)
{

    T sumA = 0;

#pragma acc parallel loop present(pgA) reduction(+ \
                                                 : sumA) copy(sumA) 
    for (arma::uword ii = 0; ii < pgA.n_elem; ii++) {
        sumA += pgA.at(ii);
    }

    return sumA;
}

template <typename T1, typename T2>
pgCol<typename T1> operator*(T2 lhval, pgCol const & rhval) {
    
}


#endif //POWER_GRID_pgCol_hpp
