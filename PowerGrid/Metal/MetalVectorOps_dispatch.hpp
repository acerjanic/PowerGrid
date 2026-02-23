/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file MetalVectorOps_dispatch.hpp
/// @brief Inline type-trait dispatch helpers called by pgCol/pgMat operators.
///
/// Provides overloaded try_metal_* functions for float and pgComplex<float>.
/// Each returns true if Metal handled the operation, false for CPU fallback.
/// For unsupported types (double, pgComplex<double>), the is_metal_type trait
/// is false, and `if constexpr` eliminates the Metal code path entirely.
///
/// Prefers zero-copy dispatch when pointers are page-aligned (pgCol uses
/// aligned_alloc(16384)).  Falls back to memcpy dispatch otherwise.
/// The minimum element count is configurable via PG_METAL_MIN_SIZE env var.

#ifndef POWER_GRID_MetalVectorOps_dispatch_hpp
#define POWER_GRID_MetalVectorOps_dispatch_hpp

#ifdef METAL_COMPUTE

#include "Metal/MetalVectorOps.h"
#include "Core/pgComplex.hpp"
#include <type_traits>
#include <cstdlib>
#include <cstdint>

namespace pg_metal {

/// Compile-time default minimum element count for Metal dispatch.
static constexpr size_t kMinMetalSize = 4096;

/// Runtime-configurable threshold.  Override with PG_METAL_MIN_SIZE env var.
inline size_t metalMinSize() {
    static size_t val = [] {
        if (const char* env = std::getenv("PG_METAL_MIN_SIZE"))
            return std::max(size_t(256), size_t(std::atol(env)));
        return kMinMetalSize;
    }();
    return val;
}

// -----------------------------------------------------------------------
// Type traits
// -----------------------------------------------------------------------
template<typename T> struct is_metal_type : std::false_type {};
template<> struct is_metal_type<float> : std::true_type {};
template<> struct is_metal_type<pgComplex<float>> : std::true_type {};

inline MetalVectorContext* ctx() { return metal_vector_get_context(); }
inline bool available() { return ctx() != nullptr; }

/// Check if pointer is page-aligned (16384 = Apple Silicon page size).
inline bool isPageAligned(const void* p) {
    return (reinterpret_cast<uintptr_t>(p) & 16383) == 0;
}

// -----------------------------------------------------------------------
// Element-wise vector-vector: float
// -----------------------------------------------------------------------
inline bool try_metal_add(const float* A, const float* B, float* C, size_t n) {
    if (n < metalMinSize()) return false;
    auto* c = ctx(); if (!c) return false;
    if (isPageAligned(A) && isPageAligned(B) && isPageAligned(C))
        metal_vec_add_zc(c, const_cast<float*>(A), const_cast<float*>(B), C, n);
    else
        metal_vec_add(c, A, B, C, n);
    return true;
}
inline bool try_metal_sub(const float* A, const float* B, float* C, size_t n) {
    if (n < metalMinSize()) return false;
    auto* c = ctx(); if (!c) return false;
    if (isPageAligned(A) && isPageAligned(B) && isPageAligned(C))
        metal_vec_sub_zc(c, const_cast<float*>(A), const_cast<float*>(B), C, n);
    else
        metal_vec_sub(c, A, B, C, n);
    return true;
}
inline bool try_metal_mul(const float* A, const float* B, float* C, size_t n) {
    if (n < metalMinSize()) return false;
    auto* c = ctx(); if (!c) return false;
    if (isPageAligned(A) && isPageAligned(B) && isPageAligned(C))
        metal_vec_mul_zc(c, const_cast<float*>(A), const_cast<float*>(B), C, n);
    else
        metal_vec_mul(c, A, B, C, n);
    return true;
}
inline bool try_metal_div(const float* A, const float* B, float* C, size_t n) {
    if (n < metalMinSize()) return false;
    auto* c = ctx(); if (!c) return false;
    if (isPageAligned(A) && isPageAligned(B) && isPageAligned(C))
        metal_vec_div_zc(c, const_cast<float*>(A), const_cast<float*>(B), C, n);
    else
        metal_vec_div(c, A, B, C, n);
    return true;
}

// -----------------------------------------------------------------------
// Element-wise vector-vector: pgComplex<float>
// -----------------------------------------------------------------------
inline bool try_metal_add(const pgComplex<float>* A, const pgComplex<float>* B,
                           pgComplex<float>* C, size_t n) {
    if (n < metalMinSize()) return false;
    auto* c = ctx(); if (!c) return false;
    if (isPageAligned(A) && isPageAligned(B) && isPageAligned(C))
        metal_cvec_add_zc(c, const_cast<float*>(reinterpret_cast<const float*>(A)),
                             const_cast<float*>(reinterpret_cast<const float*>(B)),
                             reinterpret_cast<float*>(C), n);
    else
        metal_cvec_add(c, reinterpret_cast<const float*>(A),
                          reinterpret_cast<const float*>(B),
                          reinterpret_cast<float*>(C), n);
    return true;
}
inline bool try_metal_sub(const pgComplex<float>* A, const pgComplex<float>* B,
                           pgComplex<float>* C, size_t n) {
    if (n < metalMinSize()) return false;
    auto* c = ctx(); if (!c) return false;
    if (isPageAligned(A) && isPageAligned(B) && isPageAligned(C))
        metal_cvec_sub_zc(c, const_cast<float*>(reinterpret_cast<const float*>(A)),
                             const_cast<float*>(reinterpret_cast<const float*>(B)),
                             reinterpret_cast<float*>(C), n);
    else
        metal_cvec_sub(c, reinterpret_cast<const float*>(A),
                          reinterpret_cast<const float*>(B),
                          reinterpret_cast<float*>(C), n);
    return true;
}
inline bool try_metal_mul(const pgComplex<float>* A, const pgComplex<float>* B,
                           pgComplex<float>* C, size_t n) {
    if (n < metalMinSize()) return false;
    auto* c = ctx(); if (!c) return false;
    if (isPageAligned(A) && isPageAligned(B) && isPageAligned(C))
        metal_cvec_mul_zc(c, const_cast<float*>(reinterpret_cast<const float*>(A)),
                             const_cast<float*>(reinterpret_cast<const float*>(B)),
                             reinterpret_cast<float*>(C), n);
    else
        metal_cvec_mul(c, reinterpret_cast<const float*>(A),
                          reinterpret_cast<const float*>(B),
                          reinterpret_cast<float*>(C), n);
    return true;
}
inline bool try_metal_div(const pgComplex<float>* A, const pgComplex<float>* B,
                           pgComplex<float>* C, size_t n) {
    if (n < metalMinSize()) return false;
    auto* c = ctx(); if (!c) return false;
    if (isPageAligned(A) && isPageAligned(B) && isPageAligned(C))
        metal_cvec_div_zc(c, const_cast<float*>(reinterpret_cast<const float*>(A)),
                             const_cast<float*>(reinterpret_cast<const float*>(B)),
                             reinterpret_cast<float*>(C), n);
    else
        metal_cvec_div(c, reinterpret_cast<const float*>(A),
                          reinterpret_cast<const float*>(B),
                          reinterpret_cast<float*>(C), n);
    return true;
}

// -----------------------------------------------------------------------
// Scalar operations: float
// -----------------------------------------------------------------------
inline bool try_metal_add_scalar(const float* A, float s, float* C, size_t n) {
    if (n < metalMinSize()) return false;
    auto* c = ctx(); if (!c) return false;
    if (isPageAligned(A) && isPageAligned(C))
        metal_vec_add_scalar_zc(c, const_cast<float*>(A), s, C, n);
    else
        metal_vec_add_scalar(c, A, s, C, n);
    return true;
}
inline bool try_metal_sub_scalar(const float* A, float s, float* C, size_t n) {
    return try_metal_add_scalar(A, -s, C, n);
}
inline bool try_metal_mul_scalar(const float* A, float s, float* C, size_t n) {
    if (n < metalMinSize()) return false;
    auto* c = ctx(); if (!c) return false;
    if (isPageAligned(A) && isPageAligned(C))
        metal_vec_mul_scalar_zc(c, const_cast<float*>(A), s, C, n);
    else
        metal_vec_mul_scalar(c, A, s, C, n);
    return true;
}

// -----------------------------------------------------------------------
// Scalar operations: pgComplex<float>
// -----------------------------------------------------------------------
inline bool try_metal_mul_scalar(const pgComplex<float>* A,
                                  const pgComplex<float>& s,
                                  pgComplex<float>* C, size_t n) {
    if (n < metalMinSize()) return false;
    auto* c = ctx(); if (!c) return false;
    if (isPageAligned(A) && isPageAligned(C))
        metal_cvec_mul_scalar_zc(c, const_cast<float*>(reinterpret_cast<const float*>(A)),
                                  s.real(), s.imag(),
                                  reinterpret_cast<float*>(C), n);
    else
        metal_cvec_mul_scalar(c, reinterpret_cast<const float*>(A),
                              s.real(), s.imag(),
                              reinterpret_cast<float*>(C), n);
    return true;
}

// -----------------------------------------------------------------------
// AXPY: C = A + alpha * B  (complex)
// -----------------------------------------------------------------------
inline bool try_metal_axpy(const pgComplex<float>* A,
                            const pgComplex<float>* B,
                            pgComplex<float>* C,
                            const pgComplex<float>& alpha, size_t n) {
    if (n < metalMinSize()) return false;
    auto* c = ctx(); if (!c) return false;
    if (isPageAligned(A) && isPageAligned(B) && isPageAligned(C))
        metal_cvec_axpy_zc(c, const_cast<float*>(reinterpret_cast<const float*>(A)),
                              const_cast<float*>(reinterpret_cast<const float*>(B)),
                              reinterpret_cast<float*>(C),
                              alpha.real(), alpha.imag(), n);
    else
        metal_cvec_axpy(c, reinterpret_cast<const float*>(A),
                           reinterpret_cast<const float*>(B),
                           reinterpret_cast<float*>(C),
                           alpha.real(), alpha.imag(), n);
    return true;
}

// -----------------------------------------------------------------------
// Real weight * complex vector: C = W .* X
// -----------------------------------------------------------------------
inline bool try_metal_rvec_cmul(const float* W,
                                 const pgComplex<float>* X,
                                 pgComplex<float>* C, size_t n) {
    if (n < metalMinSize()) return false;
    auto* c = ctx(); if (!c) return false;
    if (isPageAligned(W) && isPageAligned(X) && isPageAligned(C))
        metal_rvec_cmul_zc(c, const_cast<float*>(W),
                              const_cast<float*>(reinterpret_cast<const float*>(X)),
                              reinterpret_cast<float*>(C), n);
    else
        metal_rvec_cmul(c, W,
                        reinterpret_cast<const float*>(X),
                        reinterpret_cast<float*>(C), n);
    return true;
}

// -----------------------------------------------------------------------
// Reductions: float
// -----------------------------------------------------------------------
inline bool try_metal_sum(const float* A, float* result, size_t n) {
    if (n < metalMinSize()) return false;
    auto* c = ctx(); if (!c) return false;
    if (isPageAligned(A))
        *result = metal_vec_sum_zc(c, const_cast<float*>(A), n);
    else
        *result = metal_vec_sum(c, A, n);
    return true;
}

// -----------------------------------------------------------------------
// Reductions: pgComplex<float>
// -----------------------------------------------------------------------
inline bool try_metal_cdot(const pgComplex<float>* A,
                            const pgComplex<float>* B,
                            pgComplex<float>* result, size_t n) {
    if (n < metalMinSize()) return false;
    auto* c = ctx(); if (!c) return false;
    float re, im;
    if (isPageAligned(A) && isPageAligned(B))
        metal_cvec_cdot_zc(c, const_cast<float*>(reinterpret_cast<const float*>(A)),
                              const_cast<float*>(reinterpret_cast<const float*>(B)),
                              &re, &im, n);
    else
        metal_cvec_cdot(c, reinterpret_cast<const float*>(A),
                           reinterpret_cast<const float*>(B),
                           &re, &im, n);
    *result = pgComplex<float>(re, im);
    return true;
}

inline bool try_metal_norm2sq(const pgComplex<float>* A, float* result, size_t n) {
    if (n < metalMinSize()) return false;
    auto* c = ctx(); if (!c) return false;
    if (isPageAligned(A))
        *result = metal_cvec_norm2sq_zc(c, const_cast<float*>(reinterpret_cast<const float*>(A)), n);
    else
        *result = metal_cvec_norm2sq(c, reinterpret_cast<const float*>(A), n);
    return true;
}

} // namespace pg_metal

#endif // METAL_COMPUTE
#endif // POWER_GRID_MetalVectorOps_dispatch_hpp
