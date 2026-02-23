/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                    MRFIL Research Groups
               University of Illinois, Urbana-Champaign
*/

/// @file MetalVectorOps.h
/// @brief Pure C++ API for Metal-accelerated vector operations.
///
/// Provides a process-wide singleton MetalVectorContext that holds the Metal
/// device, command queue, and pipeline states for all vector algebra kernels.
/// All functions accept raw float* pointers (host-side, on unified memory).

#ifndef POWER_GRID_MetalVectorOps_h
#define POWER_GRID_MetalVectorOps_h

#ifdef METAL_COMPUTE

#include <cstddef>

struct MetalVectorContext;

/// Get or create the process-wide singleton context.
/// Thread-safe (uses std::call_once internally).
/// Returns nullptr if no Apple GPU is available.
MetalVectorContext* metal_vector_get_context();

// ============================================================================
// Element-wise real operations.  n = number of float elements.
// A, B are inputs; C is output.  In-place (C == A) is supported.
// ============================================================================

void metal_vec_add(MetalVectorContext* ctx,
                   const float* A, const float* B, float* C, size_t n);
void metal_vec_sub(MetalVectorContext* ctx,
                   const float* A, const float* B, float* C, size_t n);
void metal_vec_mul(MetalVectorContext* ctx,
                   const float* A, const float* B, float* C, size_t n);
void metal_vec_div(MetalVectorContext* ctx,
                   const float* A, const float* B, float* C, size_t n);
void metal_vec_add_scalar(MetalVectorContext* ctx,
                          const float* A, float scalar, float* C, size_t n);
void metal_vec_mul_scalar(MetalVectorContext* ctx,
                          const float* A, float scalar, float* C, size_t n);

// ============================================================================
// Element-wise complex operations.  n = number of complex elements.
// Buffers have 2*n floats (interleaved real/imag).
// ============================================================================

void metal_cvec_add(MetalVectorContext* ctx,
                    const float* A, const float* B, float* C, size_t n);
void metal_cvec_sub(MetalVectorContext* ctx,
                    const float* A, const float* B, float* C, size_t n);
void metal_cvec_mul(MetalVectorContext* ctx,
                    const float* A, const float* B, float* C, size_t n);
void metal_cvec_div(MetalVectorContext* ctx,
                    const float* A, const float* B, float* C, size_t n);

/// Real weight * complex: C[i] = W[i] * X[i].
/// W has n floats; X and C have 2*n floats.
void metal_rvec_cmul(MetalVectorContext* ctx,
                     const float* W, const float* X, float* C, size_t n);

/// Complex AXPY: C[i] = A[i] + (alphaRe + i*alphaIm) * B[i].
void metal_cvec_axpy(MetalVectorContext* ctx,
                     const float* A, const float* B, float* C,
                     float alphaRe, float alphaIm, size_t n);

/// Complex scalar multiply: C[i] = (alphaRe + i*alphaIm) * A[i].
void metal_cvec_mul_scalar(MetalVectorContext* ctx,
                           const float* A, float alphaRe, float alphaIm,
                           float* C, size_t n);

// ============================================================================
// Reductions.  n = number of complex elements (buffers have 2*n floats).
// ============================================================================

/// Complex dot product: sum(conj(A[i]) * B[i]).
void metal_cvec_cdot(MetalVectorContext* ctx,
                     const float* A, const float* B,
                     float* outRe, float* outIm, size_t n);

/// L2 norm squared of complex vector: sum(|A[i]|^2).
float metal_cvec_norm2sq(MetalVectorContext* ctx, const float* A, size_t n);

/// Sum of real float vector.
float metal_vec_sum(MetalVectorContext* ctx, const float* A, size_t n);

// ============================================================================
// Zero-copy variants.  Caller memory MUST be page-aligned (e.g. pgCol's
// aligned_alloc(16384, ...)).  These create ephemeral MTLBuffers that wrap the
// caller's memory directly — no memcpy.  Used for benchmarking and for
// production paths where pgCol data is already page-aligned.
// ============================================================================

// Element-wise real
void metal_vec_add_zc(MetalVectorContext* ctx,
                      float* A, float* B, float* C, size_t n);
void metal_vec_sub_zc(MetalVectorContext* ctx,
                      float* A, float* B, float* C, size_t n);
void metal_vec_mul_zc(MetalVectorContext* ctx,
                      float* A, float* B, float* C, size_t n);
void metal_vec_div_zc(MetalVectorContext* ctx,
                      float* A, float* B, float* C, size_t n);
void metal_vec_add_scalar_zc(MetalVectorContext* ctx,
                              float* A, float scalar, float* C, size_t n);
void metal_vec_mul_scalar_zc(MetalVectorContext* ctx,
                              float* A, float scalar, float* C, size_t n);

// Element-wise complex
void metal_cvec_add_zc(MetalVectorContext* ctx,
                       float* A, float* B, float* C, size_t n);
void metal_cvec_sub_zc(MetalVectorContext* ctx,
                       float* A, float* B, float* C, size_t n);
void metal_cvec_mul_zc(MetalVectorContext* ctx,
                       float* A, float* B, float* C, size_t n);
void metal_cvec_div_zc(MetalVectorContext* ctx,
                       float* A, float* B, float* C, size_t n);
void metal_rvec_cmul_zc(MetalVectorContext* ctx,
                        float* W, float* X, float* C, size_t n);
void metal_cvec_axpy_zc(MetalVectorContext* ctx,
                        float* A, float* B, float* C,
                        float alphaRe, float alphaIm, size_t n);
void metal_cvec_mul_scalar_zc(MetalVectorContext* ctx,
                               float* A, float alphaRe, float alphaIm,
                               float* C, size_t n);

// Reductions
void metal_cvec_cdot_zc(MetalVectorContext* ctx,
                        float* A, float* B,
                        float* outRe, float* outIm, size_t n);
float metal_cvec_norm2sq_zc(MetalVectorContext* ctx, float* A, size_t n);
float metal_vec_sum_zc(MetalVectorContext* ctx, float* A, size_t n);

// ============================================================================
// Dispatch statistics — counts command buffer commits and cumulative GPU wait.
// ============================================================================

/// Total number of Metal command buffer commits (vector ops only, not gridding).
uint64_t metal_vecops_dispatch_count();

/// Cumulative time in seconds spent in waitUntilCompleted (vector ops only).
double metal_vecops_wait_seconds();

/// Reset both counters to zero.
void metal_vecops_reset_stats();

#endif // METAL_COMPUTE
#endif // POWER_GRID_MetalVectorOps_h
