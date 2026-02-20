/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                    MRFIL Research Groups
               University of Illinois, Urbana-Champaign
*/

/// @file vectorops_metal.metal
/// @brief Metal compute shaders for pgCol/pgMat vector algebra.
///
/// Element-wise kernels operate on raw float* buffers.
/// Complex values are interleaved: [re0, im0, re1, im1, ...].
/// Reduction kernels use threadgroup shared memory for per-group partial sums;
/// the CPU finishes the small final reduction.

#include <metal_stdlib>
using namespace metal;

// ============================================================================
// Element-wise real kernels — one thread per float element
// ============================================================================

kernel void pg_vec_add(
    device const float* A [[buffer(0)]],
    device const float* B [[buffer(1)]],
    device float*       C [[buffer(2)]],
    uint gid [[thread_position_in_grid]])
{
    C[gid] = A[gid] + B[gid];
}

kernel void pg_vec_sub(
    device const float* A [[buffer(0)]],
    device const float* B [[buffer(1)]],
    device float*       C [[buffer(2)]],
    uint gid [[thread_position_in_grid]])
{
    C[gid] = A[gid] - B[gid];
}

kernel void pg_vec_mul(
    device const float* A [[buffer(0)]],
    device const float* B [[buffer(1)]],
    device float*       C [[buffer(2)]],
    uint gid [[thread_position_in_grid]])
{
    C[gid] = A[gid] * B[gid];
}

kernel void pg_vec_div(
    device const float* A [[buffer(0)]],
    device const float* B [[buffer(1)]],
    device float*       C [[buffer(2)]],
    uint gid [[thread_position_in_grid]])
{
    C[gid] = A[gid] / B[gid];
}

kernel void pg_vec_add_scalar(
    device const float* A [[buffer(0)]],
    constant float&     s [[buffer(1)]],
    device float*       C [[buffer(2)]],
    uint gid [[thread_position_in_grid]])
{
    C[gid] = A[gid] + s;
}

kernel void pg_vec_mul_scalar(
    device const float* A [[buffer(0)]],
    constant float&     s [[buffer(1)]],
    device float*       C [[buffer(2)]],
    uint gid [[thread_position_in_grid]])
{
    C[gid] = A[gid] * s;
}

// ============================================================================
// Complex element-wise kernels — one thread per complex element
// gid ranges [0, n_complex).  Buffers have 2*n_complex floats.
// ============================================================================

kernel void pg_cvec_add(
    device const float* A [[buffer(0)]],
    device const float* B [[buffer(1)]],
    device float*       C [[buffer(2)]],
    uint gid [[thread_position_in_grid]])
{
    uint idx = gid * 2;
    C[idx]     = A[idx]     + B[idx];
    C[idx + 1] = A[idx + 1] + B[idx + 1];
}

kernel void pg_cvec_sub(
    device const float* A [[buffer(0)]],
    device const float* B [[buffer(1)]],
    device float*       C [[buffer(2)]],
    uint gid [[thread_position_in_grid]])
{
    uint idx = gid * 2;
    C[idx]     = A[idx]     - B[idx];
    C[idx + 1] = A[idx + 1] - B[idx + 1];
}

kernel void pg_cvec_mul(
    device const float* A [[buffer(0)]],
    device const float* B [[buffer(1)]],
    device float*       C [[buffer(2)]],
    uint gid [[thread_position_in_grid]])
{
    uint idx = gid * 2;
    float ar = A[idx], ai = A[idx + 1];
    float br = B[idx], bi = B[idx + 1];
    C[idx]     = ar * br - ai * bi;
    C[idx + 1] = ar * bi + ai * br;
}

kernel void pg_cvec_div(
    device const float* A [[buffer(0)]],
    device const float* B [[buffer(1)]],
    device float*       C [[buffer(2)]],
    uint gid [[thread_position_in_grid]])
{
    uint idx = gid * 2;
    float ar = A[idx], ai = A[idx + 1];
    float br = B[idx], bi = B[idx + 1];
    float denom = br * br + bi * bi;
    C[idx]     = (ar * br + ai * bi) / denom;
    C[idx + 1] = (ai * br - ar * bi) / denom;
}

/// Real weight * Complex element-wise: C[i] = W[i] * X[i]
/// W has n floats (real weights), X and C have 2*n floats (complex).
kernel void pg_rvec_cmul(
    device const float* W [[buffer(0)]],
    device const float* X [[buffer(1)]],
    device float*       C [[buffer(2)]],
    uint gid [[thread_position_in_grid]])
{
    float w = W[gid];
    uint idx = gid * 2;
    C[idx]     = w * X[idx];
    C[idx + 1] = w * X[idx + 1];
}

/// Complex AXPY: C[i] = A[i] + alpha * B[i]
/// alpha is passed as constant float2 (re, im).
kernel void pg_cvec_axpy(
    device const float*  A     [[buffer(0)]],
    device const float*  B     [[buffer(1)]],
    device float*        C     [[buffer(2)]],
    constant float2&     alpha [[buffer(3)]],
    uint gid [[thread_position_in_grid]])
{
    uint idx = gid * 2;
    float br = B[idx], bi = B[idx + 1];
    float pr = alpha.x * br - alpha.y * bi;
    float pi = alpha.x * bi + alpha.y * br;
    C[idx]     = A[idx]     + pr;
    C[idx + 1] = A[idx + 1] + pi;
}

/// Complex scalar multiply: C[i] = alpha * A[i]
kernel void pg_cvec_mul_scalar(
    device const float*  A     [[buffer(0)]],
    constant float2&     alpha [[buffer(1)]],
    device float*        C     [[buffer(2)]],
    uint gid [[thread_position_in_grid]])
{
    uint idx = gid * 2;
    float ar = A[idx], ai = A[idx + 1];
    C[idx]     = alpha.x * ar - alpha.y * ai;
    C[idx + 1] = alpha.x * ai + alpha.y * ar;
}

// ============================================================================
// Reduction kernels — per-threadgroup partial sums, CPU finishes.
// Threadgroup size = 256.
// ============================================================================

/// Sum of float vector.  Writes one partial sum per threadgroup.
kernel void pg_vec_reduce_sum(
    device const float* A       [[buffer(0)]],
    device float*       partial [[buffer(1)]],
    constant uint&      N       [[buffer(2)]],
    uint gid    [[thread_position_in_grid]],
    uint tid    [[thread_index_in_threadgroup]],
    uint tgSize [[threads_per_threadgroup]],
    uint tgID   [[threadgroup_position_in_grid]])
{
    threadgroup float sdata[256];
    sdata[tid] = (gid < N) ? A[gid] : 0.0f;
    threadgroup_barrier(mem_flags::mem_threadgroup);

    for (uint s = tgSize / 2; s > 0; s >>= 1) {
        if (tid < s) sdata[tid] += sdata[tid + s];
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    if (tid == 0) partial[tgID] = sdata[0];
}

/// Complex dot product: partial sums of conj(A) * B.
/// Writes partial_re[tgID] and partial_im[tgID].
/// gid ranges [0, N) where N = number of complex elements.
kernel void pg_cvec_reduce_cdot(
    device const float* A          [[buffer(0)]],
    device const float* B          [[buffer(1)]],
    device float*       partial_re [[buffer(2)]],
    device float*       partial_im [[buffer(3)]],
    constant uint&      N          [[buffer(4)]],
    uint gid    [[thread_position_in_grid]],
    uint tid    [[thread_index_in_threadgroup]],
    uint tgSize [[threads_per_threadgroup]],
    uint tgID   [[threadgroup_position_in_grid]])
{
    threadgroup float sre[256];
    threadgroup float sim[256];

    if (gid < N) {
        uint idx = gid * 2;
        float ar = A[idx], ai = A[idx + 1];
        float br = B[idx], bi = B[idx + 1];
        // conj(A)*B = (ar - i*ai)(br + i*bi) = (ar*br + ai*bi) + i(ar*bi - ai*br)
        sre[tid] = ar * br + ai * bi;
        sim[tid] = ar * bi - ai * br;
    } else {
        sre[tid] = 0.0f;
        sim[tid] = 0.0f;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    for (uint s = tgSize / 2; s > 0; s >>= 1) {
        if (tid < s) {
            sre[tid] += sre[tid + s];
            sim[tid] += sim[tid + s];
        }
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    if (tid == 0) {
        partial_re[tgID] = sre[0];
        partial_im[tgID] = sim[0];
    }
}

/// L2 norm squared of complex vector: partial sums of |A[i]|^2.
/// gid ranges [0, N) where N = number of complex elements.
kernel void pg_cvec_reduce_norm2sq(
    device const float* A       [[buffer(0)]],
    device float*       partial [[buffer(1)]],
    constant uint&      N       [[buffer(2)]],
    uint gid    [[thread_position_in_grid]],
    uint tid    [[thread_index_in_threadgroup]],
    uint tgSize [[threads_per_threadgroup]],
    uint tgID   [[threadgroup_position_in_grid]])
{
    threadgroup float sdata[256];
    if (gid < N) {
        uint idx = gid * 2;
        sdata[tid] = A[idx] * A[idx] + A[idx + 1] * A[idx + 1];
    } else {
        sdata[tid] = 0.0f;
    }
    threadgroup_barrier(mem_flags::mem_threadgroup);

    for (uint s = tgSize / 2; s > 0; s >>= 1) {
        if (tid < s) sdata[tid] += sdata[tid + s];
        threadgroup_barrier(mem_flags::mem_threadgroup);
    }
    if (tid == 0) partial[tgID] = sdata[0];
}
