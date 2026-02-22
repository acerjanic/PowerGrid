/*
   (C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
   All rights reserved.

   See LICENSE.txt for the University of Illinois/NCSA Open Source license.

   Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
 */

/*****************************************************************************

    File Name   [TimeSegmentation.cpp]

    Synopsis    [Wrappers to the cuFFT library supporting single and double
        precision for GPU accelerated FFTs]

    Description []

    Revision    [0.1.0; Giang-Chau Ngo, BIOE UIUC]

    Date        [4/19/2016]

*****************************************************************************/
#include "TimeSegmentation.h"
// This using field correction by time segmentation
// The data is corrected to time 0 with reference to the time vector passed
//
//

// Safe scalar min/max that work correctly under -ffast-math.
// Armadillo's vectorized min()/max() can return garbage when compiled with
// -Ofast/-ffast-math due to auto-vectorized SIMD reductions.
namespace {
template<typename T>
T safe_min(const arma::Col<T>& v) {
    T m = v(0);
    for (arma::uword i = 1; i < v.n_elem; i++) {
        if (v(i) < m) m = v(i);
    }
    return m;
}
template<typename T>
T safe_max(const arma::Col<T>& v) {
    T m = v(0);
    for (arma::uword i = 1; i < v.n_elem; i++) {
        if (v(i) > m) m = v(i);
    }
    return m;
}
} // anonymous namespace

// We are using two template types at the moment. One for the type of data to be
// processed (ie Col<cx_double>) and one for the type of G object (ie
// Gfft<Col<cx_double>>
// T1 is the data type for complex, T2 is the data type for real data
template <typename T1, typename Tobj>
TimeSegmentation<T1, Tobj>::TimeSegmentation(Tobj& G, Col<T1> map_in,
    Col<T1> timeVec_in, uword a,
    uword b, uword nTimeSegs, uword interptype,
    uword shots)
{
    //cout << "Entering TimeSegmenatation Class constructor" << endl;
    n1 = a; // Data size
    n2 = b; // Image size
    // L == 0 or 1 is ambiguous. In matlab we would often specific L = 0 to mean no time segmentation.
    if ((nTimeSegs == 0) || (nTimeSegs == 1)) {
        L = 1;
    } else {
        L = nTimeSegs;
    } // number of time segments
    type = interptype; // type of time segmentation performed
    Nshots = shots; // number of shots
    obj = &G;
    fieldMap = map_in;
    //cout << "N1 = " << n1 << endl;
    //cout << "N2 = " << n2 << endl;
    //cout << "L = " << L << endl;

    outData.set_size(n1, L);
    outImg.set_size(n2, L);
    tempD.set_size(n2, L);
    tempAD.set_size(n1, L);

    AA.set_size(n1, L + 1); // time segments weights
    timeVec = timeVec_in;
    T_min = safe_min(timeVec);
    T1 rangt = safe_max(timeVec) - T_min;
    tau = (rangt + datum::eps) / L;
    timeVec = timeVec - T_min;

    uword NOneShot = n1 / Nshots;
    if (L <= 1) {
        L = 1;
        tau = 0;
        AA.ones();
        Wo.ones(n2, 1);
        WoH.ones(n2, 1);
    } else {
        Mat<complex<T1>> tempAA(NOneShot, L);
        if (type == 1) { // Hanning interpolator
            cout << "Using Hanning window temporal interpolator" << endl;
            cout << "Field Map size = " << this->fieldMap.n_rows << endl;
            for (unsigned int ii = 0; ii < L; ii++) {
                for (unsigned int jj = 0; jj < NOneShot; jj++) {
                    if ((abs(timeVec(jj) - ((ii)*tau))) <= tau) {
                        tempAA(jj, ii) = 0.5 + 0.5 * std::cos((datum::pi) * (timeVec(jj) - ((ii)*tau)) / tau);
                    } else {
                        tempAA(jj, ii) = 0.0;
                    }
                }
            }
            AA = repmat(tempAA, Nshots, 1);
            //AA.save("AA.csv", arma::csv_ascii);
        } else if (type == 2) { // Min-max interpolator: Exact LS interpolator

            cout << "Using exact LS minmax temporal interpolator" << endl;

            Mat<complex<T1>> Ltp;
            Ltp.ones(1, L);
            Col<complex<T1>> ggtp;
            ggtp.ones(n2, 1);
            Mat<complex<T1>> gg;
            gg = exp(i * fieldMap * tau) * Ltp;
            Mat<complex<T1>> iGTGGT;
            iGTGGT.set_size(L + 1, n2);
            Mat<complex<T1>> gl;
            gl.zeros(n2, L);

            for (unsigned int ii = 0; ii < L; ii++) {
                for (unsigned int jj = 0; jj < n2; jj++) {
                    gl(jj, ii) = pow(gg(jj, ii), (T1)(ii + 1));
                }
            }

            Mat<complex<T1>> G;
            G.set_size(n2, L);

            for (unsigned int jj = 0; jj < L; jj++) {
                if (jj == 0) {
                    G.col(jj) = ggtp;
                } else {
                    G.col(jj) = gl.col(jj - 1);
                }
            }

            Col<complex<T1>> glsum;
            Mat<complex<T1>> GTG;
            GTG.zeros(L, L);
            GTG.diag(0) += n2;
            glsum = sum(gl.t(), 1);
            Mat<complex<T1>> GTGtp(L, L);
            for (unsigned int ii = 0; ii < (L - 1); ii++) {
                GTGtp.zeros();
                GTGtp.diag(-(T1)(ii + 1)) += glsum(ii);
                GTGtp.diag((T1)(ii + 1)) += std::conj(glsum(ii));
                GTG = GTG + GTGtp;
            }

            T1 rcn = 1 / cond(GTG);
            if (rcn > 10 * 2e-16) { // condition number of GTG
                iGTGGT = inv(GTG) * G.t();

            } else {
                iGTGGT = pinv(GTG) * G.t(); // pseudo inverse
            }

            Mat<complex<T1>> iGTGGTtp;
            Mat<complex<T1>> ftp;
            Col<complex<T1>> res, temp;

            for (unsigned int ii = 0; ii < NOneShot; ii++) {
                ftp = exp(i * fieldMap * timeVec(ii));
                res = iGTGGT * ftp;
                tempAA.row(ii) = res.t();
            }
            AA = repmat(tempAA, Nshots, 1);
        } else if (type == 3) { // Approximate minmax estimator
            // Estimate histogram of field map
            auto start = std::chrono::high_resolution_clock::now();

            std::cout << "Using approximate minmax temporal interpolator." << std::endl;
            int numBins = 256;
            Col<uword> fm_histo = hist(vectorise(fieldMap), numBins).eval();

            T1 minFM = min(vectorise(fieldMap));
            T1 maxFM = max(vectorise(fieldMap));

            T1 rangeFM = maxFM - minFM;
            std::cout << "rangeFM = " << rangeFM << std::endl;

            Col<T1> bin_edges = linspace<Col<T1>>(minFM, maxFM, numBins);

            int KK = floor(2.0 * datum::pi / ((rangeFM / (T1)numBins) * tau));
            std::cout << "KK = " << KK << std::endl;
            T1 dwn = 2 * datum::pi / (KK * tau);

            //cx_vec ftwe_ap = vectorise(fft(fm_histo % (exp(i * dwn * regspace(0,(numBins-1))*tau*L)),KK));
            Col<CxT1> temp = fm_histo % (exp(i * dwn * regspace<Col<CxT1>>(0, (numBins - 1)) * tau * L));
            Col<CxT1> ftwe_ap = vectorise(calcFFT1D(temp, KK));

            ftwe_ap = (exp(-i * (minFM + dwn / 2) * tau * (regspace<Col<CxT1>>(0, KK - 1) - L)) % ftwe_ap);

            Mat<CxT1> GTGap_ap = zeros<Mat<CxT1>>(L + 1, L + 1);

            for (int ii = 1; ii <= (2 * L + 1); ii++) {
                GTGap_ap += diagmat(ftwe_ap(ii - 1) * ones<Mat<CxT1>>(L + 1 - abs((L + 1) - (ii)), 1), -(L + 1) + (ii));
            }
            Mat<CxT1> iGTGap_ap;

            if (rcond(GTGap_ap.st()) > 10 * datum::eps) {
                iGTGap_ap = inv(GTGap_ap.st());
            } else {
                iGTGap_ap = pinv(GTGap_ap.st().eval());
                std::cout << "Used pinv instead" << std::endl;
            }
            Col<CxT1> ftc_ap, GTc_ap;

            //#pragma omp parallel for shared(tempAA) private(ftc_ap, GTc_ap)
            for (uword ii = 0; ii < NOneShot; ii++) {
                temp = fm_histo % exp(i * regspace<Col<CxT1>>(0, (numBins - 1)) * dwn * timeVec(ii));
                ftc_ap = vectorise(calcFFT1D(temp, KK));
                ftc_ap = exp(i * (minFM + dwn / 2) * (timeVec(ii) - (tau * regspace<Col<CxT1>>(0, KK - 1)))) % ftc_ap;
                GTc_ap = ftc_ap(span(0, L));
                tempAA.row(ii) = conv_to<Row<CxT1>>::from((iGTGap_ap * GTc_ap).t());
            }
            AA = repmat(tempAA, Nshots, 1);
            auto finish = std::chrono::high_resolution_clock::now();

            std::chrono::duration<float> elapsed = finish - start;
            std::cout << "Time Segmentation Histo Elapsed time: " << elapsed.count() << " s\n";
        }

        Wo.set_size(n2, L + 1);
        WoH.set_size(n2, L + 1);
        for (unsigned int ii = 0; ii < L + 1; ii++) {
            Wo.col(ii) = exp(-i * (this->fieldMap) * ((ii) * this->tau + this->T_min));
            WoH.col(ii) = exp(i * (this->fieldMap) * ((ii) * this->tau + this->T_min));
        }
    }

#ifdef METAL_COMPUTE
    if constexpr (std::is_same<T1, float>::value) {
        AA_pg = pgMat<pgComplex<T1>>(AA);
        Wo_pg = pgMat<pgComplex<T1>>(Wo);
        WoH_pg = pgMat<pgComplex<T1>>(WoH);
        outData_pg = pgMat<pgComplex<T1>>(n1, L);
        outImg_pg = pgMat<pgComplex<T1>>(n2, L);
        tempD_pg = pgMat<pgComplex<T1>>(n2, L);
        tempAD_pg = pgMat<pgComplex<T1>>(n1, L);
    }
#endif
    //cout << "Exiting class constructor." << endl;
}

// Use FFT
template <typename T1, typename Tobj>
inline Col<complex<T1>> TimeSegmentation<T1, Tobj>::
    calcFFT1D(const Col<CxT1>& d, uword K) const
{
    arma::Col<CxT1> data;
    // Start by dealing with zero padding or trimming
    if (K < d.n_rows) {
        data = d(span(0, K - 1));
    } else {
        data.set_size(K, 1);
        data.zeros();
        data(span(0, d.n_rows - 1)) = d;
    }
    data.eval();
    // Now we can do an FF without any additional
    // zero padding or trimming
    T1* cplxRawData = reinterpret_cast<T1*>(data.memptr());

    fft1dCPU<T1>(cplxRawData, K);

    return data;
}

// Overloaded operators go here

// Forward transformation is *
// d is the vector of data of type T1, note it is const, so we don't modify it
// directly rather return another vector of type T1
template <typename T1, typename Tobj>
Col<complex<T1>> TimeSegmentation<T1, Tobj>::
operator*(const Col<complex<T1>>& d) const
{
    RANGE(__FUNCTION__)

    Tobj* G = this->obj;

#ifdef METAL_COMPUTE
    if constexpr (std::is_same<T1, float>::value) {
        // Metal path: pgCol/pgMat with GPU-dispatched element-wise ops
        pgCol<pgComplex<T1>> d_pg(d);

        // Copy Wo columns into tempD and weight by d
        for (unsigned int ii = 0; ii < this->L; ii++) {
            // Copy Wo column ii into tempD, then multiply by d
            pgCol<pgComplex<T1>> col_copy = Wo_pg.col_copy(ii);
            col_copy %= d_pg;
            tempD_pg.set_col(ii, col_copy);
        }

        // Forward NUFFT per segment (pgCol overload — no arma conversion)
        for (unsigned int ii = 0; ii < this->L; ii++) {
            pgCol<pgComplex<T1>> seg = tempD_pg.col_copy(ii);
            outData_pg.set_col(ii, (*G) * seg);
        }

        // Weight by interpolation coefficients
        for (unsigned int ii = 0; ii < this->L; ii++) {
            outData_pg.col(ii) %= AA_pg.col(ii);
        }

        pgCol<pgComplex<T1>> sumVec = sum(outData_pg, 1);
        return sumVec.getArma();
    }
#endif

    // Armadillo path (double, or non-Metal builds)
    tempD = Wo;

    for (unsigned int ii = 0; ii < this->L; ii++) {
        tempD.col(ii) %= d;
    }

    for (unsigned int ii = 0; ii < this->L; ii++) {
        outData.col(ii) = (*G * tempD.col(ii));
    }

    for (unsigned int ii = 0; ii < this->L; ii++) {
        outData.col(ii) %= AA.col(ii);
    }

    return sum(outData, 1);
}

template <typename T1, typename Tobj>
Col<complex<T1>> TimeSegmentation<T1, Tobj>::
operator/(const Col<complex<T1>>& d) const
{
    RANGE(__FUNCTION__)

    Tobj* G = this->obj;

#ifdef METAL_COMPUTE
    if constexpr (std::is_same<T1, float>::value) {
        // Metal path: pgCol/pgMat with GPU-dispatched element-wise ops
        pgCol<pgComplex<T1>> d_pg(d);

        // Conjugate AA and weight by data
        pgMat<pgComplex<T1>> conjAA_pg = conj(AA_pg);
        for (unsigned int ii = 0; ii < this->L; ii++) {
            pgCol<pgComplex<T1>> col_copy = conjAA_pg.col_copy(ii);
            col_copy %= d_pg;
            tempAD_pg.set_col(ii, col_copy);
        }

        // Adjoint NUFFT per segment (pgCol overload — no arma conversion)
        for (unsigned int ii = 0; ii < this->L; ii++) {
            pgCol<pgComplex<T1>> seg = tempAD_pg.col_copy(ii);
            outImg_pg.set_col(ii, (*G) / seg);
        }

        // Weight by conjugate field map
        for (unsigned int ii = 0; ii < this->L; ii++) {
            outImg_pg.col(ii) %= WoH_pg.col(ii);
        }

        pgCol<pgComplex<T1>> sumVec = sum(outImg_pg, 1);
        return sumVec.getArma();
    }
#endif

    // Armadillo path (double, or non-Metal builds)
    tempAD = conj(AA);

    for (unsigned int ii = 0; ii < this->L; ii++) {
        tempAD.col(ii) %= d;
    }

    for (unsigned int ii = 0; ii < this->L; ii++) {
        outImg.col(ii) = ((*G) / tempAD.col(ii));
    }

    for (unsigned int ii = 0; ii < this->L; ii++) {
        outImg.col(ii) %= WoH.col(ii);
    }

    return sum(outImg, 1);
}

// pgCol forward: image → k-space with field correction
template <typename T1, typename Tobj>
pgCol<pgComplex<T1>> TimeSegmentation<T1, Tobj>::
operator*(const pgCol<pgComplex<T1>>& d) const
{
#ifdef METAL_COMPUTE
    if constexpr (std::is_same<T1, float>::value) {
        Tobj* G = this->obj;

        for (unsigned int ii = 0; ii < this->L; ii++) {
            pgCol<pgComplex<T1>> col_copy = Wo_pg.col_copy(ii);
            col_copy %= d;
            tempD_pg.set_col(ii, col_copy);
        }

        for (unsigned int ii = 0; ii < this->L; ii++) {
            pgCol<pgComplex<T1>> seg = tempD_pg.col_copy(ii);
            outData_pg.set_col(ii, (*G) * seg);
        }

        for (unsigned int ii = 0; ii < this->L; ii++) {
            outData_pg.col(ii) %= AA_pg.col(ii);
        }

        return sum(outData_pg, 1);
    }
#endif
    // Fallback: convert to arma, call arma operator, convert back
    Col<CxT1> armaResult = this->operator*(d.getArma());
    return pgCol<pgComplex<T1>>(armaResult);
}

// pgCol adjoint: k-space → image with field correction
template <typename T1, typename Tobj>
pgCol<pgComplex<T1>> TimeSegmentation<T1, Tobj>::
operator/(const pgCol<pgComplex<T1>>& d) const
{
#ifdef METAL_COMPUTE
    if constexpr (std::is_same<T1, float>::value) {
        Tobj* G = this->obj;

        pgMat<pgComplex<T1>> conjAA_pg = conj(AA_pg);
        for (unsigned int ii = 0; ii < this->L; ii++) {
            pgCol<pgComplex<T1>> col_copy = conjAA_pg.col_copy(ii);
            col_copy %= d;
            tempAD_pg.set_col(ii, col_copy);
        }

        for (unsigned int ii = 0; ii < this->L; ii++) {
            pgCol<pgComplex<T1>> seg = tempAD_pg.col_copy(ii);
            outImg_pg.set_col(ii, (*G) / seg);
        }

        for (unsigned int ii = 0; ii < this->L; ii++) {
            outImg_pg.col(ii) %= WoH_pg.col(ii);
        }

        return sum(outImg_pg, 1);
    }
#endif
    // Fallback: convert to arma, call arma operator, convert back
    Col<CxT1> armaResult = this->operator/(d.getArma());
    return pgCol<pgComplex<T1>>(armaResult);
}

// Explicit Instantiations
template class TimeSegmentation<float, Gnufft<float>>;
//template Col<complex<float>> TimeSegmentation<float, Gnufft<float> >::operator/(arma::Col<std::complex<float> > const &) const;
//template Col<complex<float>> TimeSegmentation<float, Gnufft<float> >::operator*(arma::Col<std::complex<float> > const &) const;
template class TimeSegmentation<float, Gdft<float>>;
template class TimeSegmentation<double, Gnufft<double>>;
//template Col<complex<double>> TimeSegmentation<double, Gnufft<double> >::operator/(arma::Col<std::complex<double>> const &) const;
//template Col<complex<double>> TimeSegmentation<double, Gnufft<double> >::operator*(arma::Col<std::complex<double>> const &) const;
template class TimeSegmentation<double, Gdft<double>>;
