// pgSubview_Col.hpp

#include "../PGIncludes.h"

template<typename T>
class pgSubviewCol<T> {

private:

T *mem; //Pointer to raw data from source object
unsigned int uiColHeader; // Number of the first entry of the column in question
bool isInitialized; 
bool isOnGPU; 

// Number of elements in array
unsigned int n_elem;
unsigned int n_elemCol; // Number of elements in the column;
unsigned int uiColHead; // Index of first element of the subview column in the original array from the original Mat object.
public:

// This should return the length of the column we are working with, not the number of entries in the original matrix
unsigned int n_elem() const{
    return n_elemCol;
}

template<typename T_ = T>
    pgSubviewCol(T* memptr, arma::uword n_elem, arma::uword n_elemCol, arma::uword uiColHead)
    {
        this.mem = memptr;
        this.n_elem = n_elem;
        this.n_elemCol = n_elemCol;
        this.uiColHead = uiColHead;
    }

    ~pgCol<T>()
    {
        // Do not free the memory as this is just a subview of the original object and it's aliased.
    }

void zeros() {
    #pragma acc parallel loop present(mem[0:n_elem])
    for(unsigned int ii = 0; ii < n_elemCol; ii++) {
        // Note that we modify this to start use the original pointer to avoid OpenACC issues with the present table
        // And that we only work over the span of one column
        mem[uiColHeader + ii] = T();
    }
}

void ones() {
    #pragma acc parallel loop present(mem[0:n_elem])
    for(unsigned int ii = 0; ii < n_elemCol; ii++) {
        mem[uiColHeader + ii] = T(1.0);
    }
}


    inline const T at(const arma::uword d) const
    {
        return mem[uiColHeader + d];
    }

    inline T& at(const arma::uword d)
    {
        return mem[uiColHeader + d];
    }

    inline T& operator()(const arma::uword d)
    {
        //#pragma acc update_host(mem)
        return mem[uiColHeaer + d];
    }

    inline const T operator()(const arma::uword d) const
    {
        //#pragma acc update_host(mem)
        return mem[uiColHeader + d];
    }

    pgCol<T>& operator=(const pgCol<T>& d)
    {
        T* pDmem = d.memptr();
        //TODO: Fix this OpenACC data directive, can't rely on referencing d as it is an object
#pragma acc parallel loop present(mem [0:n_elem], pDmem[0:n_elemCol])
        for (arma::uword ii = 0; ii < n_elemCol; ii++) {
            this->mem[uiColHeader + ii] = pDmem[ii];
        }
        return *this;
    }

/*
// I don't think this is necessary. This is changing the size of the elements, and that's not a legal operation if we're working on a subview
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
*/
    pgCol<T>& operator+=(const T& A)
    {

#pragma acc parallel loop present(mem [0:n_elem], A)
        for (arma::uword ii = 0; ii < n_elemCol; ii++) {
            this->mem[uiColHeader + ii] += A;
        }
        return *this;
    }

    pgCol<T>& operator-=(const T& A)
    {
#pragma acc parallel loop present(this, mem [0:n_elem], A)
        for (arma::uword ii = 0; ii < n_elemCol; ii++) {
            this->mem[uiColHeader + ii] -= A;
        }
        return *this;
    }

    template <typename X>
    pgCol<T>& operator%=(const T& A)
    {
#pragma acc parallel loop present(this, mem [0:n_elem], A)
        for (arma::uword ii = 0; ii < n_elemCol; ii++) {
            this->mem[uiColHeader + ii] *= A;
        }
        return *this;
    }

    template <typename X>
    pgCol<T>& operator/=(const T& A)
    {
#pragma acc parallel loop present(this, mem [0:n_elem], A)
        for (arma::uword ii = 0; ii < n_elemCol; ii++) {
            mem[uiColHeader + ii] /= A;
        }
        return *this;
    }

    template <typename X>
    pgCol<T>& operator+=(const pgCol<X>& pgA)
    {
        T* pAmem = pgA.memptr();
// TODO: Fix this OpenACC data directive, can't refer to pgA like this.
#pragma acc parallel loop present(this, mem [0:n_elem], pAmem[0:n_elemCol])
        for (arma::uword ii = 0; ii < n_elemCol; ii++) {
            this->mem[uiColHeader + ii] += pAmem[ii];
        }
        return *this;
    }

    template <typename X>
    pgCol<T>& operator-=(const pgCol<X>& pgA)
    {

        T* pAmem = pgA.memptr();

#pragma acc parallel loop present(this, mem [0:n_elem], pAmem[0:n_elemCol])
        for (arma::uword ii = 0; ii < n_elemCol; ii++) {
            this->mem[uiColHeader + ii] -= pAmem[ii];
        }
        return *this;
    }

    template <typename X>
    pgCol<T>& operator%=(const pgCol<X>& pgA)
    {
        T* pAmem = pgA.memptr();

#pragma acc parallel loop present(this, mem [0:n_elem], pAmem[0:n_elemCol])
        for (arma::uword ii = 0; ii < n_elemCol; ii++) {
            this->mem[uiColHeader + ii] *= pAmem[ii];
        }
        return *this;
    }

    template <typename X>
    pgCol<T>& operator/=(const pgCol<X>& pgA)
    {
        T* pAmem = pgA.memptr();

#pragma acc parallel loop present(this, mem [0:n_elem], pAmem [0:n_elemCol])
        for (arma::uword ii = 0; ii < n_elemCol; ii++) {
            mem[uiColHeader + ii] /= pAmem[ii];
        }
        return *this;
    }

    // Operators(pgCol, scalar)
    //template <typename T>
    const pgCol<T> operator+(const T& B) const
    {
        pgCol<T> pgC(n_elemCol);
        T* pCmem = pgC.memptr();

#pragma acc parallel loop present(this, mem [0:n_elem], pCmem [0:n_elemCol], B)
        for (arma::uword ii = 0; ii < n_elemCol; ii++) {
            pCmem[ii] = mem[uiColHeader + ii] + B;
        }
        return std::move(pgC);
    }

    //template <typename T>
    const pgCol<T> operator-(const T& B) const
    {
        pgCol<T> pgC(n_elemCol);
        T* pCmem = pgC.memptr();

#pragma acc parallel loop present(this, mem [0:n_elem], pCmem[0:n_elemCol], B)
        for (arma::uword ii = 0; ii < n_elemCol; ii++) {
            pCmem[ii] = mem[uiColHeader + ii] - B;
        }
        return std::move(pgC);
    }

    //template <typename T>
    const pgCol<T> operator%(const T& B) const
    {
        pgCol<T> pgC(n_elemCol);
        T* pCmem = pgC.memptr();

#pragma acc parallel loop present(this, mem [0:n_elem], pCmem[0:n_elemCol], B)
        for (arma::uword ii = 0; ii < n_elemCol; ii++) {
            pCmem[ii] = mem[uiColHeader + ii] * B;
        }
        return std::move(pgC);
    }

    //template <typename T>
    const pgCol<T> operator/(const T& B) const
    {
        pgCol<T> pgC(n_elemCol);
        T* pCmem = pgC.memptr();
        T* pBmem = pgB.memptr();

#pragma acc parallel loop present(this, mem [0:n_elem], pCmem [0:n_elemCol], pBmem [0:n_elemCol], B)
        for (arma::uword ii = 0; ii < n_elemCol; ii++) {
            pgC.at(ii) = mem[uiColHeader + ii] / B;
        }
        return std::move(pgC);
    }

    // Operators(pgCol, pgCol)
    template <typename X>
    const pgCol<T> operator+(const pgCol<X>& pgB) const
    {
        pgCol<T> pgC(n_elemCol);
        T* pCmem = pgC.memptr();
        T* pBmem = pgB.memptr();
#pragma acc parallel loop present(mem [0:n_elem], pCmem [0:n_elemCol], pBmem [0:n_elemCol])
        for (arma::uword ii = 0; ii < n_elemCol; ii++) {
            pCmem[ii] = mem[uiColHeader + ii] + pBmem[ii];
        }
        return std::move(pgC);
    }

    template <typename X>
    const pgCol<T> operator-(const pgCol<X>& pgB) const
    {
        pgCol<T> pgC(n_elemCol);
        T* pCmem = pgC.memptr();
        T* pBmem = pgB.memptr();

#pragma acc parallel loop present(mem [0:n_elem], pCmem [0:n_elemCol], pBmem [0:n_elemCol])
        for (arma::uword ii = 0; ii < n_elemCol; ii++) {
            pCmem[ii] = mem[uiColHeader + ii] - pBmem[ii];
        }
        return std::move(pgC);
    }

    template <typename X>
    const pgCol<T> operator%(const pgCol<X>& pgB) const
    {
        pgCol<T> pgC(n_elemCol);
        //Obtain memory pointers to manage the openACC data directives
        T* pCmem = pgC.memptr();
        T* pBmem = pgB.memptr();

#pragma acc parallel loop present(this, mem [0:n_elem], pCmem [0:n_elemCol], pBmem [0:n_elemCol])
        for (arma::uword ii = 0; ii < n_elemCol; ii++) {
            pCmem[ii] = mem[uiColHeader + ii] * pBmem[ii];
        }

        return std::move(pgC);
    }

    template <typename X>
    const pgCol<T> operator/(const pgCol<X>& pgB) const
    {
        pgCol<T> pgC(n_elemCol);
        //Obtain memory pointers to manage the openACC data directives
        T* pCmem = pgC.memptr();
        T* pBmem = pgB.memptr();

#pragma acc parallel loop present(this, mem [0:n_elem], pCmem [0:n_elemCol], pBmem [0:n_elemCol])
        for (arma::uword ii = 0; ii < n_elemCol; ii++) {
            pCmem[ii] = mem[uiColHeader + ii] / pBmem[ii];
        }
        return std::move(pgC);
    }

};


