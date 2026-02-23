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
#include "Core/pgComplex.hpp"
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
    // Complex multiply: (a+bi)(c+di) = (ac-bd) + (ad+bc)i
    // Read inputs with stride-2, compute into contiguous temp buffers,
    // then interleave to output once.  Safe when C aliases A or B.
    const float* Af = reinterpret_cast<const float*>(A);
    const float* Bf = reinterpret_cast<const float*>(B);

    // 3 temp buffers: t_re (real result), t_im (imag result), t_tmp (scratch)
    const size_t stackThreshold = 10922;
    float stackBuf[stackThreshold * 3] __attribute__((aligned(16)));
    float* heap = nullptr;
    float* t_re, *t_im, *t_tmp;
    if (n <= stackThreshold) {
        t_re = stackBuf;
        t_im = stackBuf + n;
        t_tmp = stackBuf + 2 * n;
    } else {
        heap = (float*)malloc(3 * n * sizeof(float));
        t_re = heap;
        t_im = heap + n;
        t_tmp = heap + 2 * n;
    }

    // Real: ac - bd
    vDSP_vmul(Af, 2, Bf, 2, t_re, 1, (vDSP_Length)n);           // t_re = ac
    vDSP_vmul(Af + 1, 2, Bf + 1, 2, t_tmp, 1, (vDSP_Length)n);  // t_tmp = bd
    vDSP_vsub(t_tmp, 1, t_re, 1, t_re, 1, (vDSP_Length)n);      // t_re = ac - bd

    // Imag: ad + bc
    vDSP_vmul(Af, 2, Bf + 1, 2, t_im, 1, (vDSP_Length)n);       // t_im = ad
    vDSP_vmul(Af + 1, 2, Bf, 2, t_tmp, 1, (vDSP_Length)n);      // t_tmp = bc
    vDSP_vadd(t_im, 1, t_tmp, 1, t_im, 1, (vDSP_Length)n);      // t_im = ad + bc

    // Single interleave: [t_re, t_im] → interleaved C
    DSPSplitComplex sc = {t_re, t_im};
    vDSP_ztoc(&sc, 1, (DSPComplex*)C, 2, (vDSP_Length)n);

    if (heap) free(heap);
    return true;
}
inline bool try_accel_div(const pgComplex<float>* A, const pgComplex<float>* B,
                           pgComplex<float>* C, size_t n) {
    if (n < kMinAccelSize) return false;
    // Complex division: C = A / B = (a+bi)/(c+di)
    //   Cr = (ac + bd) / (c² + d²)
    //   Ci = (bc - ad) / (c² + d²)
    // Read inputs with stride-2, compute into contiguous temp buffers,
    // then interleave to output once.  Safe when C aliases A or B.
    const float* Af = reinterpret_cast<const float*>(A);
    const float* Bf = reinterpret_cast<const float*>(B);

    // 4 temp buffers: t_re, t_im, t_tmp, denom
    const size_t stackThreshold = 8192;
    float stackBuf[stackThreshold * 4] __attribute__((aligned(16)));
    float* heap = nullptr;
    float* t_re, *t_im, *t_tmp, *denom;
    if (n <= stackThreshold) {
        t_re = stackBuf;
        t_im = stackBuf + n;
        t_tmp = stackBuf + 2 * n;
        denom = stackBuf + 3 * n;
    } else {
        heap = (float*)malloc(4 * n * sizeof(float));
        t_re = heap;
        t_im = heap + n;
        t_tmp = heap + 2 * n;
        denom = heap + 3 * n;
    }

    // denom = c² + d²
    vDSP_vmul(Bf, 2, Bf, 2, t_re, 1, (vDSP_Length)n);           // t_re = c²
    vDSP_vmul(Bf + 1, 2, Bf + 1, 2, t_tmp, 1, (vDSP_Length)n);  // t_tmp = d²
    vDSP_vadd(t_re, 1, t_tmp, 1, denom, 1, (vDSP_Length)n);     // denom = c² + d²

    // Real: (ac + bd) / denom
    vDSP_vmul(Af, 2, Bf, 2, t_re, 1, (vDSP_Length)n);           // t_re = ac
    vDSP_vmul(Af + 1, 2, Bf + 1, 2, t_tmp, 1, (vDSP_Length)n);  // t_tmp = bd
    vDSP_vadd(t_re, 1, t_tmp, 1, t_re, 1, (vDSP_Length)n);      // t_re = ac + bd
    vDSP_vdiv(denom, 1, t_re, 1, t_re, 1, (vDSP_Length)n);      // t_re = (ac+bd)/denom

    // Imag: (bc - ad) / denom
    vDSP_vmul(Af + 1, 2, Bf, 2, t_im, 1, (vDSP_Length)n);       // t_im = bc
    vDSP_vmul(Af, 2, Bf + 1, 2, t_tmp, 1, (vDSP_Length)n);      // t_tmp = ad
    vDSP_vsub(t_tmp, 1, t_im, 1, t_im, 1, (vDSP_Length)n);      // t_im = bc - ad
    vDSP_vdiv(denom, 1, t_im, 1, t_im, 1, (vDSP_Length)n);      // t_im = (bc-ad)/denom

    // Single interleave: [t_re, t_im] → interleaved C
    DSPSplitComplex sc = {t_re, t_im};
    vDSP_ztoc(&sc, 1, (DSPComplex*)C, 2, (vDSP_Length)n);

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
    // (sr + si*i) * (a + b*i) = (sr*a - si*b) + (sr*b + si*a)*i
    // Read inputs with stride-2, compute into contiguous temp buffers,
    // then interleave to output once.  Safe when C aliases A.
    const float* Af = reinterpret_cast<const float*>(A);
    float sr = s.real(), si = s.imag();

    // 3 temp buffers: t_re, t_im, t_tmp
    const size_t stackThreshold = 10922;
    float stackBuf[stackThreshold * 3] __attribute__((aligned(16)));
    float* heap = nullptr;
    float* t_re, *t_im, *t_tmp;
    if (n <= stackThreshold) {
        t_re = stackBuf;
        t_im = stackBuf + n;
        t_tmp = stackBuf + 2 * n;
    } else {
        heap = (float*)malloc(3 * n * sizeof(float));
        t_re = heap;
        t_im = heap + n;
        t_tmp = heap + 2 * n;
    }

    // Real: sr*a - si*b
    vDSP_vsmul(Af, 2, &sr, t_re, 1, (vDSP_Length)n);            // t_re = sr*a
    vDSP_vsmul(Af + 1, 2, &si, t_tmp, 1, (vDSP_Length)n);       // t_tmp = si*b
    vDSP_vsub(t_tmp, 1, t_re, 1, t_re, 1, (vDSP_Length)n);      // t_re = sr*a - si*b

    // Imag: sr*b + si*a
    vDSP_vsmul(Af + 1, 2, &sr, t_im, 1, (vDSP_Length)n);        // t_im = sr*b
    vDSP_vsmul(Af, 2, &si, t_tmp, 1, (vDSP_Length)n);           // t_tmp = si*a
    vDSP_vadd(t_im, 1, t_tmp, 1, t_im, 1, (vDSP_Length)n);      // t_im = sr*b + si*a

    // Single interleave: [t_re, t_im] → interleaved C
    DSPSplitComplex sc = {t_re, t_im};
    vDSP_ztoc(&sc, 1, (DSPComplex*)C, 2, (vDSP_Length)n);

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
    // C = A + alpha * B, where alpha = (ar + ai*i), all complex vectors.
    // Compute alpha*B real/imag into contiguous temps, add A's real/imag
    // (stride-2 gather), then interleave to C.  Safe when C aliases A or B.
    const float* Af = reinterpret_cast<const float*>(A);
    const float* Bf = reinterpret_cast<const float*>(B);
    float ar = alpha.real(), ai = alpha.imag();

    // 3 temp buffers: t_re, t_im, t_tmp
    const size_t stackThreshold = 10922;
    float stackBuf[stackThreshold * 3] __attribute__((aligned(16)));
    float* heap = nullptr;
    float* t_re, *t_im, *t_tmp;
    if (n <= stackThreshold) {
        t_re = stackBuf;
        t_im = stackBuf + n;
        t_tmp = stackBuf + 2 * n;
    } else {
        heap = (float*)malloc(3 * n * sizeof(float));
        t_re = heap;
        t_im = heap + n;
        t_tmp = heap + 2 * n;
    }

    // Real of alpha*B: ar*Br - ai*Bi
    vDSP_vsmul(Bf, 2, &ar, t_re, 1, (vDSP_Length)n);            // t_re = ar*Br
    vDSP_vsmul(Bf + 1, 2, &ai, t_tmp, 1, (vDSP_Length)n);       // t_tmp = ai*Bi
    vDSP_vsub(t_tmp, 1, t_re, 1, t_re, 1, (vDSP_Length)n);      // t_re = ar*Br - ai*Bi

    // Imag of alpha*B: ar*Bi + ai*Br
    vDSP_vsmul(Bf + 1, 2, &ar, t_im, 1, (vDSP_Length)n);        // t_im = ar*Bi
    vDSP_vsmul(Bf, 2, &ai, t_tmp, 1, (vDSP_Length)n);           // t_tmp = ai*Br
    vDSP_vadd(t_im, 1, t_tmp, 1, t_im, 1, (vDSP_Length)n);      // t_im = ar*Bi + ai*Br

    // Add A's real and imag parts (gathered with stride-2)
    // vDSP_vadd with stride-2 on one input, stride-1 on other
    vDSP_vadd(Af, 2, t_re, 1, t_re, 1, (vDSP_Length)n);         // t_re += Ar
    vDSP_vadd(Af + 1, 2, t_im, 1, t_im, 1, (vDSP_Length)n);     // t_im += Ai

    // Single interleave: [t_re, t_im] → interleaved C
    DSPSplitComplex sc = {t_re, t_im};
    vDSP_ztoc(&sc, 1, (DSPComplex*)C, 2, (vDSP_Length)n);

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
