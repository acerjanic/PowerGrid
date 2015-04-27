//This file is distributed under a BSD 3-Clause license.
// See Ceemple_BSD_license.txt for details.

#ifndef CEEMPLEMATIO_H
#define CEEMPLEMATIO_H

#include "CeempleArmadillo.h"
#include "matio.h"

using namespace arma;

template <class T>struct MatioTypes {
    static const matio_types DataType = MAT_T_UNKNOWN;
};
template <> struct MatioTypes<std::complex<double>> {
    static const matio_types DataType = MAT_T_DOUBLE;
    static const matio_flags FlagType = MAT_F_COMPLEX;
};
template <> struct MatioTypes<std::complex<float>> {
    static const matio_types DataType = MAT_T_SINGLE;
    static const matio_flags FlagType = MAT_F_COMPLEX;
};
template <> struct MatioTypes<uint8_t> {
    static const matio_types DataType = MAT_T_UINT8;
};
template <> struct MatioTypes<float> {
    static const matio_types DataType = MAT_T_SINGLE;
};
template <> struct MatioTypes<double> {
    static const matio_types DataType = MAT_T_DOUBLE;
};



// Return true on success.
template <typename T>
bool
savemat(const std::string &FileName, const std::string &VarName,
        const T &Var,
        const typename arma_not_cx<typename T::elem_type>::result* junk = 0) {
    cout << "Calling real savemat" << endl;
    arma_ignore(junk);
    mat_t *matfp = Mat_CreateVer(FileName.c_str(), NULL, MAT_FT_MAT5);
    if (!matfp) {
        fprintf(stderr, "saveMat: could not create the file '%s'.\n",
                FileName.c_str());
        return false;
    }
    size_t dims[2] = {Var.n_rows, Var.n_cols};
    matvar_t *matvar =
            Mat_VarCreate(VarName.c_str(), MAT_C_DOUBLE, MatioTypes<typename T::elem_type>::DataType, 2,
                          dims, (void *)Var.memptr(), 0);
    if (!matvar) {
        fprintf(stderr, "saveMat: error creating variable '%s'.\n.",
                VarName.c_str());
        Mat_Close(matfp);
        return false;
    }
    Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);
    Mat_VarFree(matvar);
    Mat_Close(matfp);
    return true;
}

// Return true on success.
template <typename T1>
inline
bool
savemat(const std::string &FileName, const std::string &VarName,
        const T1 &Var,
        const typename arma_cx_only<typename T1::elem_type>::result* junk = 0) {
    arma_extra_debug_sigprint();
    arma_ignore(junk);
    cout << "Calling complex savemat" << endl;

    mat_t *matfp = Mat_CreateVer(FileName.c_str(), NULL, MAT_FT_MAT5);
    if (!matfp) {
        fprintf(stderr, "saveMat: could not create the file '%s'.\n",
                FileName.c_str());
        return false;
    }
    arma::Mat<cx_double> temp = conv_to<arma::Mat<cx_double>>::from(Var);
    size_t dims[2] = {temp.n_rows, temp.n_cols};
    arma::Mat<double> reTemp = arma::real(temp);
    arma::Mat<double> imTemp = arma::imag(temp);
    struct mat_complex_split_t z = {reTemp.memptr(),imTemp.memptr()};

    std::cout << "MatioTypes<T1>::DataType = " << MatioTypes<typename T1::elem_type>::DataType << std::endl;

    matvar_t *matvar =
            Mat_VarCreate(VarName.c_str(), MAT_C_DOUBLE, MatioTypes<typename T1::elem_type>::DataType, 2,
                          dims, &z, MAT_F_COMPLEX);
    if (!matvar) {
        fprintf(stderr, "saveMat: error creating variable '%s'.\n.",
                VarName.c_str());
        Mat_Close(matfp);
        return false;
    }
    Mat_VarWrite(matfp, matvar, MAT_COMPRESSION_ZLIB);
    Mat_VarFree(matvar);
    Mat_Close(matfp);
    return true;
}


template <typename T>
bool
loadmat(const std::string &FileName, const std::string &VarName,
        T *Var,
        const typename arma_not_cx<typename T::elem_type>::result* junk = 0) {
    arma_ignore(junk);
    cout << "Calling real loadmat" << endl;

    mat_t *matfp = Mat_Open(FileName.c_str(), MAT_ACC_RDONLY);
    if (!matfp) {
        fprintf(stderr, "loadMat: could not open the file '%s'.\n",
                FileName.c_str());
        return false;
    }
    matvar_t *matvar = Mat_VarRead(matfp, VarName.c_str());
    if (!matvar) {
        fprintf(stderr, "loadMat: variable '%s' not found in file '%s'.\n",
                VarName.c_str(), FileName.c_str());
        Mat_Close(matfp);
        return false;
    }
    bool matvar2D = (matvar->rank == 2) && (matvar->class_type == MAT_C_DOUBLE);
    if (matvar2D) {
        unsigned Rows = matvar->dims[0], Cols = matvar->dims[1];

        arma::Mat<double> out = arma::Mat<double>(Rows, Cols);
        (*Var) = conv_to<T>::from(out);
        memcpy(Var->memptr(), matvar->data, Rows * Cols * sizeof(typename T::elem_type));
    } else {
        fprintf(stderr,
                "loadMat: Variable '%s' is not 2-dimensional double matrix.\n"
                        "rank =** %d class_type = %d\n",
                VarName.c_str(), matvar->rank, matvar->class_type);
    }
    Mat_VarFree(matvar);
    Mat_Close(matfp);
    return matvar2D;
}

template <typename T1>
bool
loadmat(const std::string &FileName, const std::string &VarName,
        T1 *Var,
        const typename arma_cx_only<typename T1::elem_type>::result* junk = 0) {
    arma_ignore(junk);
    cout << "Calling complex loadmat" << endl;

    mat_t *matfp = Mat_Open(FileName.c_str(), MAT_ACC_RDONLY);
    if (!matfp) {
        fprintf(stderr, "loadMat: could not open the file '%s'.\n",
                FileName.c_str());
        return false;
    }
    matvar_t *matvar = Mat_VarRead(matfp, VarName.c_str());
    if (!matvar) {
        fprintf(stderr, "loadMat: variable '%s' not found in file '%s'.\n",
                VarName.c_str(), FileName.c_str());
        Mat_Close(matfp);
        return false;
    }
    bool matvar2D = (matvar->rank == 2) && (matvar->class_type == MAT_C_DOUBLE);
    if (matvar2D) {
        unsigned Rows = matvar->dims[0], Cols = matvar->dims[1];
        //(*Var) = arma::Mat<T1>(Rows, Cols);

        //If complex data, we get a pointer to a mat_complex_split_t where we can get real and imaginary pointers
        mat_complex_split_t* z = (mat_complex_split_t*)matvar->data;

        arma::Mat<double> re((double*)z->Re,Rows,Cols,false,true);
        arma::Mat<double> im((double*)z->Im,Rows,Cols,false,true);
        arma::Mat<arma::cx_double> out(re,im);
        (*Var) = conv_to<T1>::from(out);
        //memcpy(Var->memptr(), matvar->data, Rows * Cols * sizeof(T));
    } else {
        fprintf(stderr,
                "loadMat: Variable '%s' is not 2-dimensional double matrix.\n"
                        "rank = %d class_type = %d\n",
                VarName.c_str(), matvar->rank, matvar->class_type);
    }
    Mat_VarFree(matvar);
    Mat_Close(matfp);
    return matvar2D;
}

#endif // CEEMPLEMATIO_H//
