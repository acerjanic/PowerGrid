// pgSubviewCol.hpp
#ifndef POWER_GRID_pgSubviewCol_hpp
#define POWER_GRID_pgSubviewCol_hpp

#include "../PGIncludes.h"
#include "pgCol.hpp"

template<typename T>
class pgSubviewCol {

private:


arma::uword uiColHeader; // Number of the first entry of the column in question
bool isInitialized; 
bool isOnGPU; 
//arma::uword uiColHead; // Index of first element of the subview column in the original array from the original Mat object.

public:
T *mem; //Pointer to raw data from source object
// Number of elements in array
arma::uword n_elemOrig;
arma::uword n_rows; // Number of elements in the column;
const arma::uword n_cols = 1;

// This should return the length of the column we are working with, not the number of entries in the original matrix
arma::uword n_elem() const{
    return n_rows;
}

pgSubviewCol<T>(T* memptr, arma::uword n_elem, arma::uword n_rows, arma::uword uiColHead) {
    this->mem = memptr;
    this->n_elemOrig = n_elem;
    this->n_rows = n_rows;
    this->uiColHeader = uiColHead;
    #pragma acc enter data copyin(this[0:1])
}

~pgSubviewCol<T>(){};

void zeros() {
    #pragma acc parallel loop present(mem[0:n_elemOrig])
    for(arma::uword ii = 0; ii < n_rows; ii++) {
        // Note that we modify this to start use the original pointer to avoid OpenACC issues with the present table
        // And that we only work over the span of one column
        mem[uiColHeader + ii] = T();
    }
}

void ones() {
    #pragma acc parallel loop present(mem[0:n_elemOrig])
    for(arma::uword ii = 0; ii < n_rows; ii++) {
        mem[uiColHeader + ii] = T(1.0);
    }
}

    #pragma acc routine seq
    inline const T at(const arma::uword d) const
    {
        return mem[uiColHeader + d];
    }

    #pragma acc routine seq
    inline T& at(const arma::uword d)
    {
        return mem[uiColHeader + d];
    }

    #pragma acc routine seq
    inline T& operator()(const arma::uword d)
    {
        //#pragma acc update_host(mem)
        return mem[uiColHeaer + d];
    }

    #pragma acc routine seq
    inline const T operator()(const arma::uword d) const
    {
        //#pragma acc update_host(mem)
        return mem[uiColHeader + d];
    }

    pgSubviewCol<T>& operator=(const pgCol<T>& d)
    {
        T* pDmem = d.memptr();
        //TODO: Fix this OpenACC data directive, can't rely on referencing d as it is an object
#pragma acc parallel loop present(mem [0:n_elemOrig], pDmem[0:n_rows])
        for (arma::uword ii = 0; ii < n_rows; ii++) {
            this->mem[uiColHeader + ii] = pDmem[ii];
        }
        return *this;
    }

/*
// I don't think this is necessary. This is changing the size of the elements, and that's not a legal operation if we're working on a subview
    pgCol<T>& operator=(pgCol<T>&& d)
    {
        //size_t bytes = sizeof(T) * d.n_elemOrig;
        if (isInitialized) {
#pragma acc exit data delete (mem [0:n_elemOrig])
            delete[] mem;
            access::rw(n_elemOrig) = 0;
            isInitialized = false;
            isOnGPU = false;
        }
        access::rw(n_elemOrig) = d.n_elemOrig;
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
*/

    pgCol<T>& operator+=(const T& A)
    {

#pragma acc parallel loop present(mem [0:n_elemOrig], A)
        for (arma::uword ii = 0; ii < n_rows; ii++) {
            this->mem[uiColHeader + ii] += A;
        }
        return *this;
    }

    pgCol<T>& operator-=(const T& A)
    {
#pragma acc parallel loop present(this, mem [0:n_elemOrig], A)
        for (arma::uword ii = 0; ii < n_rows; ii++) {
            this->mem[uiColHeader + ii] -= A;
        }
        return *this;
    }

    template <typename X>
    pgCol<T>& operator%=(const T& A)
    {
#pragma acc parallel loop present(this, mem [0:n_elemOrig], A)
        for (arma::uword ii = 0; ii < n_rows; ii++) {
            this->mem[uiColHeader + ii] *= A;
        }
        return *this;
    }

    template <typename X>
    pgCol<T>& operator/=(const T& A)
    {
#pragma acc parallel loop present(this, mem [0:n_elemOrig], A)
        for (arma::uword ii = 0; ii < n_rows; ii++) {
            mem[uiColHeader + ii] /= A;
        }
        return *this;
    }

    template <typename X>
    pgSubviewCol<T>& operator+=(const pgCol<X>& pgA)
    {
        T* pAmem = pgA.memptr();
// TODO: Fix this OpenACC data directive, can't refer to pgA like this.
#pragma acc parallel loop present(this, mem [0:n_elemOrig], pAmem[0:n_rows])
        for (arma::uword ii = 0; ii < n_rows; ii++) {
            this->mem[uiColHeader + ii] += pAmem[ii];
        }
        return *this;
    }

    template <typename X>
    pgCol<T>& operator-=(const pgCol<X>& pgA)
    {

        T* pAmem = pgA.memptr();

#pragma acc parallel loop present(this, mem [0:n_elemOrig], pAmem[0:n_rows])
        for (arma::uword ii = 0; ii < n_rows; ii++) {
            this->mem[uiColHeader + ii] -= pAmem[ii];
        }
        return *this;
    }

    template <typename X>
    pgCol<T>& operator%=(const pgCol<X>& pgA)
    {
        T* pAmem = pgA.memptr();

#pragma acc parallel loop present(this, mem [0:n_elemOrig], pAmem[0:n_rows])
        for (arma::uword ii = 0; ii < n_rows; ii++) {
            this->mem[uiColHeader + ii] *= pAmem[ii];
        }
        return *this;
    }

    template <typename X>
    pgCol<T>& operator/=(const pgCol<X>& pgA)
    {
        T* pAmem = pgA.memptr();

#pragma acc parallel loop present(this, mem [0:n_elemOrig], pAmem [0:n_rows])
        for (arma::uword ii = 0; ii < n_rows; ii++) {
            mem[uiColHeader + ii] /= pAmem[ii];
        }
        return *this;
    }

    // Operators(pgCol, scalar)
    //template <typename T>
    const pgCol<T> operator+(const T& B) const
    {
        pgCol<T> pgC(n_rows);
        T* pCmem = pgC.memptr();

#pragma acc parallel loop present(this, mem [0:n_elemOrig], pCmem [0:n_rows], B)
        for (arma::uword ii = 0; ii < n_rows; ii++) {
            pCmem[ii] = mem[uiColHeader + ii] + B;
        }
        return std::move(pgC);
    }

    //template <typename T>
    const pgCol<T> operator-(const T& B) const
    {
        pgCol<T> pgC(n_rows);
        T* pCmem = pgC.memptr();

#pragma acc parallel loop present(this, mem [0:n_elemOrig], pCmem[0:n_rows], B)
        for (arma::uword ii = 0; ii < n_rows; ii++) {
            pCmem[ii] = mem[uiColHeader + ii] - B;
        }
        return std::move(pgC);
    }

    //template <typename T>
    const pgCol<T> operator%(const T& B) const
    {
        pgCol<T> pgC(n_rows);
        T* pCmem = pgC.memptr();

#pragma acc parallel loop present(this, mem [0:n_elemOrig], pCmem[0:n_rows], B)
        for (arma::uword ii = 0; ii < n_rows; ii++) {
            pCmem[ii] = mem[uiColHeader + ii] * B;
        }
        return std::move(pgC);
    }

    //template <typename T>
    const pgCol<T> operator/(const T& B) const
    {
        pgCol<T> pgC(n_rows);
        T* pCmem = pgC.memptr();
        T* pBmem = pgB.memptr();

#pragma acc parallel loop present(this, mem [0:n_elemOrig], pCmem [0:n_rows], pBmem [0:n_rows], B)
        for (arma::uword ii = 0; ii < n_rows; ii++) {
            pgC.at(ii) = mem[uiColHeader + ii] / B;
        }
        return std::move(pgC);
    }

    // Operators(pgSubviewCol, pgCol)
    template <typename X>
    const pgCol<T> operator+(const pgCol<X>& pgB) const
    {
        pgCol<T> pgC(n_rows);
        T* pCmem = pgC.memptr();
        T* pBmem = pgB.memptr();
#pragma acc parallel loop present(mem [0:n_elemOrig], pCmem [0:n_rows], pBmem [0:n_rows])
        for (arma::uword ii = 0; ii < n_rows; ii++) {
            pCmem[ii] = mem[uiColHeader + ii] + pBmem[ii];
        }
        return std::move(pgC);
    }

    template <typename X>
    const pgCol<T> operator-(const pgCol<X>& pgB) const
    {
        pgCol<T> pgC(n_rows);
        T* pCmem = pgC.memptr();
        T* pBmem = pgB.memptr();

#pragma acc parallel loop present(mem [0:n_elemOrig], pCmem [0:n_rows], pBmem [0:n_rows])
        for (arma::uword ii = 0; ii < n_rows; ii++) {
            pCmem[ii] = mem[uiColHeader + ii] - pBmem[ii];
        }
        return std::move(pgC);
    }

    template <typename X>
    const pgCol<T> operator%(const pgCol<X>& pgB) const
    {
        pgCol<T> pgC(n_rows);
        //Obtain memory pointers to manage the openACC data directives
        T* pCmem = pgC.memptr();
        T* pBmem = pgB.memptr();

#pragma acc parallel loop copyin(this [0:1]) present( mem [0:n_elemOrig], pCmem [0:n_rows], pBmem [0:n_rows])
        for (arma::uword ii = 0; ii < n_rows; ii++) {
            pCmem[ii] = mem[uiColHeader + ii] * pBmem[ii];
        }

        return std::move(pgC);
    }

    template <typename X>
    const pgCol<T> operator/(const pgCol<X>& pgB) const
    {
        pgCol<T> pgC(n_rows);
        //Obtain memory pointers to manage the openACC data directives
        T* pCmem = pgC.memptr();
        T* pBmem = pgB.memptr();

#pragma acc parallel loop copyin(this [0:1]) present(mem [0:n_elemOrig], pCmem [0:n_rows], pBmem [0:n_rows])
        for (arma::uword ii = 0; ii < n_rows; ii++) {
            pCmem[ii] = mem[uiColHeader + ii] / pBmem[ii];
        }
        return std::move(pgC);
    }

    // Operators(pgSubviewCol, pgSubviewCol)
    template <typename X>
    const pgCol<T> operator+(const pgSubviewCol<X>& pgB) const
    {
        pgCol<T> pgC(n_rows);
        T* pCmem = pgC.memptr();
        T* pBmem = pgB.mem;
        //Need the size of the memory from the second operand just for OpenACC purposes
        arma::uword n_elem2 = pgB.n_elemOrig;
#pragma acc parallel loop copyin(this [0:1]) present(mem [0:n_elemOrig], pCmem [0:n_rows], pBmem [0:n_elem2])
        for (arma::uword ii = 0; ii < n_rows; ii++) {
            pCmem[ii] = mem[uiColHeader + ii] + pBmem[ii];
        }
        return std::move(pgC);
    }

    template <typename X>
    const pgCol<T> operator-(const pgSubviewCol<X>& pgB) const
    {
        pgCol<T> pgC(n_rows);
        T* pCmem = pgC.memptr();
        T* pBmem = pgB.mem;
        //Need the size of the memory from the second operand just for OpenACC purposes
        arma::uword n_elem2 = pgB.n_elemOrig;
#pragma acc parallel loop copyin(this [0:1]) present(mem [0:n_elemOrig], pCmem [0:n_rows], pBmem [0:n_elem2])
        for (arma::uword ii = 0; ii < n_rows; ii++) {
            pCmem[ii] = mem[uiColHeader + ii] - pBmem[ii];
        }
        return std::move(pgC);
    }

    template <typename X>
    const pgCol<T> operator%(const pgSubviewCol<X>& pgB) const
    {
        pgCol<T> pgC(n_rows);
        //Obtain memory pointers to manage the openACC data directives
        T* pCmem = pgC.memptr();
        T* pBmem = pgB.mem;
        //Need the size of the memory from the second operand just for OpenACC purposes
        arma::uword n_elem2 = pgB.n_elemOrig;
#pragma acc parallel loop copyin(this [0:1]) present(mem [0:n_elemOrig], pCmem [0:n_rows], pBmem [0:n_elem2])
        for (arma::uword ii = 0; ii < n_rows; ii++) {
            pCmem[ii] = mem[uiColHeader + ii] * pBmem[ii];
        }

        return std::move(pgC);
    }

    template <typename X>
    const pgCol<T> operator/(const pgSubviewCol<X>& pgB) const
    {
        pgCol<T> pgC(n_rows);
        //Obtain memory pointers to manage the openACC data directives
        T* pCmem = pgC.memptr();
        T* pBmem = pgB.mem;
        //Need the size of the memory from the second operand just for OpenACC purposes
        arma::uword n_elem2 = pgB.n_elemOrig;
#pragma acc parallel loop copyin(this [0:1]) present(mem [0:n_elemOrig], pCmem [0:n_rows], pBmem [0:n_elem2])
        for (arma::uword ii = 0; ii < n_rows; ii++) {
            pCmem[ii] = mem[uiColHeader + ii] / pBmem[ii];
        }
        return std::move(pgC);
    }


};

template <typename T>
const std::complex<T> sum(const pgSubviewCol<std::complex<T>>& pgA)
{
    T sumReal = 0;
    T sumImag = 0;

#pragma acc parallel loop present(pgA) present(pgA.mem[0:pgA.n_elemOrig]) reduction(+ \
                                                 : sumReal)
    for (arma::uword ii = 0; ii < pgA.n_rows; ii++) {
        sumReal += real(pgA.at(ii));
    }

#pragma acc parallel loop present(pgA) present(pgA.mem[0:pgA.n_elemOrig]) reduction(+ \
                                                 : sumImag)
    for (arma::uword ii = 0; ii < pgA.n_rows; ii++) {
        sumImag += imag(pgA.at(ii));
    }

    return std::complex<T>(sumReal, sumImag);
}

template <typename T>
const T sum(const pgSubviewCol<T>& pgA)
{
    T sumA = 0;

#pragma acc parallel loop copyin(pgA) present(pgA.mem[0:pgA.n_elemOrig]) reduction(+ : sumA)
    for (arma::uword ii = 0; ii < pgA.n_rows; ii++) {
        sumA += pgA.at(ii);
    }
    return sumA;
}

#endif //POWER_GRID_pgSubviewCol_hpp

