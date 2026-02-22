/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file AccelerateDispatch.hpp
/// @brief Inline dispatch helpers calling Apple Accelerate (vDSP/BLAS) for
///        pgCol element-wise operations.
///
/// Replaces the Metal GPU dispatch for element-wise ops, which benchmarking
/// showed to be 2-100x slower than Accelerate at all vector sizes due to
/// GPU command buffer overhead.
///
/// Guarded by __APPLE__ so that even METAL_COMPUTE=OFF builds on macOS get
/// SIMD-optimized pgCol operations.

#ifndef POWER_GRID_AccelerateDispatch_hpp
#define POWER_GRID_AccelerateDispatch_hpp

#ifdef __APPLE__

#include <Accelerate/Accelerate.h>
#include "pgComplex.hpp"
#include <type_traits>
#include <cstring>

namespace pg_accel {

/// Minimum element count for Accelerate dispatch.  Accelerate has near-zero
/// overhead, but very small vectors (< 32 elements) may not benefit.
static constexpr size_t kMinAccelSize = 32;

// -----------------------------------------------------------------------
// Type traits — same types as the old Metal dispatch
// -----------------------------------------------------------------------
template<typename T> struct is_accel_type : std::false_type {};
template<> struct is_accel_type<float> : std::true_type {};
template<> struct is_accel_type<pgComplex<float>> : std::true_type {};

// -----------------------------------------------------------------------
// Element-wise vector-vector: float
// -----------------------------------------------------------------------
inline bool try_accel_add(const float* A, const float* B, float* C, size_t n) {
    if (n < kMinAccelSize) return false;
    vDSP_vadd(A, 1, B, 1, C, 1, (vDSP_Length)n);
    return true;
}
inline bool try_accel_sub(const float* A, const float* B, float* C, size_t n) {
    if (n < kMinAccelSize) return false;
    // vDSP_vsub computes C = B - A (reversed!), so pass args swapped
    vDSP_vsub(B, 1, A, 1, C, 1, (vDSP_Length)n);
    return true;
}
inline bool try_accel_mul(const float* A, const float* B, float* C, size_t n) {
    if (n < kMinAccelSize) return false;
    vDSP_vmul(A, 1, B, 1, C, 1, (vDSP_Length)n);
    return true;
}
inline bool try_accel_div(const float* A, const float* B, float* C, size_t n) {
    if (n < kMinAccelSize) return false;
    // vDSP_vdiv computes C = B / A (denominator first!), so pass args swapped
    vDSP_vdiv(B, 1, A, 1, C, 1, (vDSP_Length)n);
    return true;
}

// -----------------------------------------------------------------------
// Element-wise vector-vector: pgComplex<float>
// Interleaved layout: 2n floats for n complex elements.
// For add/sub, treat as 2n real floats.
// For mul/div, use split-complex vDSP_zvmul/manual.
// -----------------------------------------------------------------------
inline bool try_accel_add(const pgComplex<float>* A, const pgComplex<float>* B,
                           pgComplex<float>* C, size_t n) {
    if (n < kMinAccelSize) return false;
    // Interleaved complex: just add 2n floats
    vDSP_vadd(reinterpret_cast<const float*>(A), 1,
              reinterpret_cast<const float*>(B), 1,
              reinterpret_cast<float*>(C), 1, (vDSP_Length)(2 * n));
    return true;
}
inline bool try_accel_sub(const pgComplex<float>* A, const pgComplex<float>* B,
                           pgComplex<float>* C, size_t n) {
    if (n < kMinAccelSize) return false;
    vDSP_vsub(reinterpret_cast<const float*>(B), 1,
              reinterpret_cast<const float*>(A), 1,
              reinterpret_cast<float*>(C), 1, (vDSP_Length)(2 * n));
    return true;
}
inline bool try_accel_mul(const pgComplex<float>* A, const pgComplex<float>* B,
                           pgComplex<float>* C, size_t n) {
    if (n < kMinAccelSize) return false;
    // Complex multiply using split-complex vDSP_zvmul.
    // Interleaved → split, multiply, split → interleaved.
    // vDSP_ctoz / vDSP_ztoc handle the conversion with stride.
    DSPSplitComplex scA, scB, scC;
    // Use stride-2 views into the interleaved arrays.
    // vDSP_zvmul with stride works on split-complex, so we need temp buffers.
    // For performance, use a stack buffer for small N, heap for large.
    const size_t stackThreshold = 16384; // 64KB on stack (4 * 16K floats)
    float stackBuf[stackThreshold * 4] __attribute__((aligned(16)));
    float* heap = nullptr;
    float* buf;
    if (n <= stackThreshold) {
        buf = stackBuf;
    } else {
        heap = (float*)malloc(4 * n * sizeof(float));
        buf = heap;
    }
    scA.realp = buf;
    scA.imagp = buf + n;
    scB.realp = buf + 2 * n;
    scB.imagp = buf + 3 * n;

    // Deinterleave
    vDSP_ctoz((const DSPComplex*)A, 2, &scA, 1, (vDSP_Length)n);
    vDSP_ctoz((const DSPComplex*)B, 2, &scB, 1, (vDSP_Length)n);

    // Can write result into scA (reuse buffer)
    scC.realp = scA.realp;
    scC.imagp = scA.imagp;
    vDSP_zvmul(&scA, 1, &scB, 1, &scC, 1, (vDSP_Length)n, 1); // 1 = no conjugate

    // Re-interleave
    vDSP_ztoc(&scC, 1, (DSPComplex*)C, 2, (vDSP_Length)n);

    if (heap) free(heap);
    return true;
}
inline bool try_accel_div(const pgComplex<float>* A, const pgComplex<float>* B,
                           pgComplex<float>* C, size_t n) {
    if (n < kMinAccelSize) return false;
    // Complex division: C = A / B
    // vDSP_zvdiv computes C = B / A (denominator first), so swap args
    const size_t stackThreshold = 16384;
    float stackBuf[stackThreshold * 4] __attribute__((aligned(16)));
    float* heap = nullptr;
    float* buf;
    if (n <= stackThreshold) {
        buf = stackBuf;
    } else {
        heap = (float*)malloc(4 * n * sizeof(float));
        buf = heap;
    }
    DSPSplitComplex scA = {buf, buf + n};
    DSPSplitComplex scB = {buf + 2*n, buf + 3*n};
    vDSP_ctoz((const DSPComplex*)A, 2, &scA, 1, (vDSP_Length)n);
    vDSP_ctoz((const DSPComplex*)B, 2, &scB, 1, (vDSP_Length)n);

    DSPSplitComplex scC = {scA.realp, scA.imagp};
    // vDSP_zvdiv: C = B / A, so pass denominator (B) as first arg
    vDSP_zvdiv(&scB, 1, &scA, 1, &scC, 1, (vDSP_Length)n);

    vDSP_ztoc(&scC, 1, (DSPComplex*)C, 2, (vDSP_Length)n);
    if (heap) free(heap);
    return true;
}

// -----------------------------------------------------------------------
// Scalar operations: float
// -----------------------------------------------------------------------
inline bool try_accel_add_scalar(const float* A, float s, float* C, size_t n) {
    if (n < kMinAccelSize) return false;
    vDSP_vsadd(A, 1, &s, C, 1, (vDSP_Length)n);
    return true;
}
inline bool try_accel_sub_scalar(const float* A, float s, float* C, size_t n) {
    float neg_s = -s;
    return try_accel_add_scalar(A, neg_s, C, n);
}
inline bool try_accel_mul_scalar(const float* A, float s, float* C, size_t n) {
    if (n < kMinAccelSize) return false;
    vDSP_vsmul(A, 1, &s, C, 1, (vDSP_Length)n);
    return true;
}

// -----------------------------------------------------------------------
// Scalar operations: pgComplex<float>
// Complex scalar multiply: C[i] = alpha * A[i]
// -----------------------------------------------------------------------
inline bool try_accel_mul_scalar(const pgComplex<float>* A,
                                  const pgComplex<float>& s,
                                  pgComplex<float>* C, size_t n) {
    if (n < kMinAccelSize) return false;
    // (sr + i*si) * (ar + i*ai) = (sr*ar - si*ai) + i*(sr*ai + si*ar)
    // Use vDSP on interleaved: multiply real parts by sr, imag parts by sr,
    // then adjust for cross terms.
    // Simplest correct approach: deinterleave, scale, reinterleave.
    const size_t stackThreshold = 16384;
    float stackBuf[stackThreshold * 2] __attribute__((aligned(16)));
    float* heap = nullptr;
    float* buf;
    if (n <= stackThreshold) {
        buf = stackBuf;
    } else {
        heap = (float*)malloc(2 * n * sizeof(float));
        buf = heap;
    }
    DSPSplitComplex scA = {buf, buf + n};
    vDSP_ctoz((const DSPComplex*)A, 2, &scA, 1, (vDSP_Length)n);

    // Scale by complex scalar using split complex
    DSPSplitComplex scS;
    float sr = s.real(), si = s.imag();
    scS.realp = &sr;
    scS.imagp = &si;
    // vDSP_zvzsml: C = A * s (element-wise by scalar)
    DSPSplitComplex scC = {scA.realp, scA.imagp}; // reuse buffer
    vDSP_zvzsml(&scA, 1, &scS, &scC, 1, (vDSP_Length)n);

    vDSP_ztoc(&scC, 1, (DSPComplex*)C, 2, (vDSP_Length)n);
    if (heap) free(heap);
    return true;
}

// -----------------------------------------------------------------------
// AXPY: C = A + alpha * B  (complex)
// -----------------------------------------------------------------------
inline bool try_accel_axpy(const pgComplex<float>* A,
                            const pgComplex<float>* B,
                            pgComplex<float>* C,
                            const pgComplex<float>& alpha, size_t n) {
    if (n < kMinAccelSize) return false;
    // Two-step: tmp = alpha * B, then C = A + tmp
    const size_t stackThreshold = 8192;
    float stackBuf[stackThreshold * 6] __attribute__((aligned(16)));
    float* heap = nullptr;
    float* buf;
    if (n <= stackThreshold) {
        buf = stackBuf;
    } else {
        heap = (float*)malloc(6 * n * sizeof(float));
        buf = heap;
    }
    DSPSplitComplex scA = {buf, buf + n};
    DSPSplitComplex scB = {buf + 2*n, buf + 3*n};
    DSPSplitComplex scC = {buf + 4*n, buf + 5*n};
    vDSP_ctoz((const DSPComplex*)A, 2, &scA, 1, (vDSP_Length)n);
    vDSP_ctoz((const DSPComplex*)B, 2, &scB, 1, (vDSP_Length)n);

    // scC = alpha * scB
    DSPSplitComplex scAlpha;
    float ar = alpha.real(), ai = alpha.imag();
    scAlpha.realp = &ar;
    scAlpha.imagp = &ai;
    vDSP_zvzsml(&scB, 1, &scAlpha, &scC, 1, (vDSP_Length)n);

    // scC = scA + scC
    vDSP_vadd(scA.realp, 1, scC.realp, 1, scC.realp, 1, (vDSP_Length)n);
    vDSP_vadd(scA.imagp, 1, scC.imagp, 1, scC.imagp, 1, (vDSP_Length)n);

    vDSP_ztoc(&scC, 1, (DSPComplex*)C, 2, (vDSP_Length)n);
    if (heap) free(heap);
    return true;
}

// -----------------------------------------------------------------------
// Real weight * complex vector: C = W .* X
// W has n floats, X and C have n pgComplex<float> (2n floats).
// -----------------------------------------------------------------------
inline bool try_accel_rvec_cmul(const float* W,
                                 const pgComplex<float>* X,
                                 pgComplex<float>* C, size_t n) {
    if (n < kMinAccelSize) return false;
    // W[i] * (Xr[i] + i*Xi[i]) = W[i]*Xr[i] + i*W[i]*Xi[i]
    // Multiply real parts and imag parts separately using stride-2 access
    const float* Xf = reinterpret_cast<const float*>(X);
    float* Cf = reinterpret_cast<float*>(C);
    // Real parts: Cf[2*i] = W[i] * Xf[2*i]
    vDSP_vmul(W, 1, Xf, 2, Cf, 2, (vDSP_Length)n);
    // Imag parts: Cf[2*i+1] = W[i] * Xf[2*i+1]
    vDSP_vmul(W, 1, Xf + 1, 2, Cf + 1, 2, (vDSP_Length)n);
    return true;
}

// -----------------------------------------------------------------------
// Reductions: float
// -----------------------------------------------------------------------
inline bool try_accel_sum(const float* A, float* result, size_t n) {
    if (n < kMinAccelSize) return false;
    vDSP_sve(A, 1, result, (vDSP_Length)n);
    return true;
}

// -----------------------------------------------------------------------
// Reductions: pgComplex<float>
// -----------------------------------------------------------------------
inline bool try_accel_cdot(const pgComplex<float>* A,
                            const pgComplex<float>* B,
                            pgComplex<float>* result, size_t n) {
    if (n < kMinAccelSize) return false;
    // cblas_cdotc_sub computes conj(A) . B
    // pgComplex<float> is layout-compatible with std::complex<float>
    std::complex<float> r;
    cblas_cdotc_sub((int)n,
                    reinterpret_cast<const std::complex<float>*>(A), 1,
                    reinterpret_cast<const std::complex<float>*>(B), 1,
                    &r);
    *result = pgComplex<float>(r.real(), r.imag());
    return true;
}

inline bool try_accel_norm2sq(const pgComplex<float>* A, float* result, size_t n) {
    if (n < kMinAccelSize) return false;
    // cblas_scnrm2 returns sqrt(sum(|A[i]|^2)), so square it
    float nrm = cblas_scnrm2((int)n, reinterpret_cast<const std::complex<float>*>(A), 1);
    *result = nrm * nrm;
    return true;
}

} // namespace pg_accel

#endif // __APPLE__
#endif // POWER_GRID_AccelerateDispatch_hpp
