/*
(C) Copyright 2015-2016 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/*****************************************************************************

    File Name   [fftGPU.hpp]

    Synopsis    [Wrappers to the cuFFT library supporting single and double
                                precision for GPU accelerated FFTs]

    Description []

    Revision    [0.2.0; Alex Cerjanic, BIOE UIUC]
                [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [12/13/2016]

 *****************************************************************************/

/// @file fftGPU.h
/// @brief cuFFT-based GPU FFT wrappers for 1-D, 2-D, and 3-D transforms.
///
/// All functions are only available when compiled with `OPENACC_GPU` defined.
/// Each function operates on a GPU-resident interleaved real/imag array via
/// an OpenACC stream and a pre-allocated cuFFT plan.

// Based on
// https://www.olcf.ornl.gov/tutorials/mixing-openacc-with-gpu-libraries/

#ifndef PowerGrid_fftGPU_hpp
#define PowerGrid_fftGPU_hpp

#ifdef OPENACC_GPU
#include "cufft.h"
#include <complex>
#include <type_traits>

// Like Armadillo, we're using SFINAE here to choose between float and double.
// (Maybe FP16 some day in the future)
// We need enable_if to choose which version to run based on the type of the
// template parameter.

/// @brief In-place 1-D forward FFT on GPU (single precision).
///
/// @tparam T1    Must be `float`.
/// @param d_data GPU pointer to interleaved real/imag float array, length 2·nx.
/// @param nx     Number of complex samples.
/// @param stream OpenACC/CUDA stream pointer.
/// @param plan   Pointer to a pre-allocated cuFFT plan.
template <typename T1, typename std::enable_if<std::is_same<T1, float>::value,
                                               int>::type = 0>
void fft1dGPU(T1 *d_data, int nx, void *stream, cufftHandle *plan);

/// @brief In-place 2-D inverse FFT on GPU (single precision).
///
/// @tparam T1    Must be `float`.
/// @param d_data GPU pointer to interleaved real/imag float array, length 2·nx·ny.
/// @param nx     Transform size in x.
/// @param ny     Transform size in y.
/// @param stream OpenACC/CUDA stream pointer.
/// @param plan   Pointer to a pre-allocated cuFFT plan.
template <typename T1, typename std::enable_if<std::is_same<T1, float>::value,
                                               int>::type = 0>
void ifft2dGPU(T1 *d_data, int nx, int ny, void *stream, cufftHandle *plan);

/// @brief In-place 2-D forward FFT on GPU (single precision).
///
/// @tparam T1    Must be `float`.
/// @param d_data GPU pointer to interleaved real/imag float array, length 2·nx·ny.
/// @param nx     Transform size in x.
/// @param ny     Transform size in y.
/// @param stream OpenACC/CUDA stream pointer.
/// @param plan   Pointer to a pre-allocated cuFFT plan.
template <typename T1, typename std::enable_if<std::is_same<T1, float>::value,
                                               int>::type = 0>
void fft2dGPU(T1 *d_data, int nx, int ny, void *stream, cufftHandle *plan);

/// @brief In-place 3-D inverse FFT on GPU (single precision).
///
/// @tparam T1    Must be `float`.
/// @param d_data GPU pointer to interleaved real/imag float array.
/// @param nx     Transform size in x.
/// @param ny     Transform size in y.
/// @param nz     Transform size in z.
/// @param stream OpenACC/CUDA stream pointer.
/// @param plan   Pointer to a pre-allocated cuFFT plan.
template <typename T1, typename std::enable_if<std::is_same<T1, float>::value,
                                               int>::type = 0>
void ifft3dGPU(T1 *d_data, int nx, int ny, int nz, void *stream, cufftHandle *plan);

/// @brief In-place 3-D forward FFT on GPU (single precision).
///
/// @tparam T1    Must be `float`.
/// @param d_data GPU pointer to interleaved real/imag float array.
/// @param nx     Transform size in x.
/// @param ny     Transform size in y.
/// @param nz     Transform size in z.
/// @param stream OpenACC/CUDA stream pointer.
/// @param plan   Pointer to a pre-allocated cuFFT plan.
template <typename T1, typename std::enable_if<std::is_same<T1, float>::value,
                                               int>::type = 0>
void fft3dGPU(T1 *d_data, int nx, int ny, int nz, void *stream, cufftHandle *plan);

/// @brief In-place 1-D forward FFT on GPU (double precision).
///
/// @tparam T1    Must be `double`.
/// @param d_data GPU pointer to interleaved real/imag double array, length 2·nx.
/// @param nx     Number of complex samples.
/// @param stream OpenACC/CUDA stream pointer.
/// @param plan   Pointer to a pre-allocated cuFFT plan.
template <typename T1, typename std::enable_if<std::is_same<T1, double>::value,
                                               int>::type = 0>
void fft1dGPU(T1 *d_data, int nx, void *stream, cufftHandle *plan);

/// @brief In-place 2-D inverse FFT on GPU (double precision).
///
/// @tparam T1    Must be `double`.
/// @param d_data GPU pointer to interleaved real/imag double array, length 2·nx·ny.
/// @param nx     Transform size in x.
/// @param ny     Transform size in y.
/// @param stream OpenACC/CUDA stream pointer.
/// @param plan   Pointer to a pre-allocated cuFFT plan.
template <typename T1, typename std::enable_if<std::is_same<T1, double>::value,
                                               int>::type = 0>
void ifft2dGPU(T1 *d_data, int nx, int ny, void *stream, cufftHandle *plan);

/// @brief In-place 2-D forward FFT on GPU (double precision).
///
/// @tparam T1    Must be `double`.
/// @param d_data GPU pointer to interleaved real/imag double array, length 2·nx·ny.
/// @param nx     Transform size in x.
/// @param ny     Transform size in y.
/// @param stream OpenACC/CUDA stream pointer.
/// @param plan   Pointer to a pre-allocated cuFFT plan.
template <typename T1, typename std::enable_if<std::is_same<T1, double>::value,
                                               int>::type = 0>
void fft2dGPU(T1 *d_data, int nx, int ny, void *stream, cufftHandle *plan);

/// @brief In-place 3-D inverse FFT on GPU (double precision).
///
/// @tparam T1    Must be `double`.
/// @param d_data GPU pointer to interleaved real/imag double array.
/// @param nx     Transform size in x.
/// @param ny     Transform size in y.
/// @param nz     Transform size in z.
/// @param stream OpenACC/CUDA stream pointer.
/// @param plan   Pointer to a pre-allocated cuFFT plan.
template <typename T1, typename std::enable_if<std::is_same<T1, double>::value,
                                               int>::type = 0>
void ifft3dGPU(T1 *d_data, int nx, int ny, int nz, void *stream, cufftHandle *plan);

/// @brief In-place 3-D forward FFT on GPU (double precision).
///
/// @tparam T1    Must be `double`.
/// @param d_data GPU pointer to interleaved real/imag double array.
/// @param nx     Transform size in x.
/// @param ny     Transform size in y.
/// @param nz     Transform size in z.
/// @param stream OpenACC/CUDA stream pointer.
/// @param plan   Pointer to a pre-allocated cuFFT plan.
template <typename T1, typename std::enable_if<std::is_same<T1, double>::value,
                                               int>::type = 0>
void fft3dGPU(T1 *d_data, int nx, int ny, int nz, void *stream, cufftHandle *plan);

// Explicit Instantiations
extern template void fft1dGPU<float>(float *, int, void *, cufftHandle *);
extern template void fft1dGPU<double>(double *, int, void *, cufftHandle *);

extern template void ifft2dGPU<float>(float *, int, int, void *, cufftHandle *);
extern template void ifft2dGPU<double>(double *, int, int, void *, cufftHandle *);

extern template void fft2dGPU<float>(float *, int, int, void *, cufftHandle *);
extern template void fft2dGPU<double>(double *, int, int, void *, cufftHandle *);

extern template void ifft3dGPU<float>(float *, int, int, int, void *, cufftHandle *);
extern template void ifft3dGPU<double>(double *, int, int, int, void *, cufftHandle *);

extern template void fft3dGPU<float>(float *, int, int, int, void *, cufftHandle *);
extern template void fft3dGPU<double>(double *, int, int, int, void *, cufftHandle *);

#endif //OPENACC_GPU
#endif // PowerGrid_fftGPU_hpp
