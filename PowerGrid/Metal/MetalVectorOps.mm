/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file MetalVectorOps.mm
/// @brief Objective-C++ implementation of Metal-accelerated vector operations.
///
/// Process-wide singleton MetalVectorContext holds the Metal device, command
/// queue, and pipeline states for all vector algebra kernels.  Scratch buffers
/// are grow-only and reused across calls to avoid per-call allocation overhead.
///
/// The compiled Metal library (vectorops.metallib) is embedded as a C byte
/// array via the generated header `vectorops_metallib.h`.

#ifdef METAL_COMPUTE

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>

#include "MetalVectorOps.h"

// Embedded metallib (generated at build time)
#include "vectorops_metallib.h"

#include <cstdio>
#include <cstring>
#include <mutex>

static constexpr NSUInteger kTGSize = 256;

// ---------------------------------------------------------------------------
// MetalVectorContext — process-wide singleton
// ---------------------------------------------------------------------------
struct MetalVectorContext {
    id<MTLDevice>       device;
    id<MTLCommandQueue> queue;
    id<MTLLibrary>      library;

    // Element-wise real pipelines
    id<MTLComputePipelineState> ps_vec_add;
    id<MTLComputePipelineState> ps_vec_sub;
    id<MTLComputePipelineState> ps_vec_mul;
    id<MTLComputePipelineState> ps_vec_div;
    id<MTLComputePipelineState> ps_vec_add_scalar;
    id<MTLComputePipelineState> ps_vec_mul_scalar;

    // Element-wise complex pipelines
    id<MTLComputePipelineState> ps_cvec_add;
    id<MTLComputePipelineState> ps_cvec_sub;
    id<MTLComputePipelineState> ps_cvec_mul;
    id<MTLComputePipelineState> ps_cvec_div;
    id<MTLComputePipelineState> ps_rvec_cmul;
    id<MTLComputePipelineState> ps_cvec_axpy;
    id<MTLComputePipelineState> ps_cvec_mul_scalar;

    // Reduction pipelines
    id<MTLComputePipelineState> ps_reduce_sum;
    id<MTLComputePipelineState> ps_reduce_cdot;
    id<MTLComputePipelineState> ps_reduce_norm2sq;

    // Reusable scratch buffers (grow-only, shared memory)
    id<MTLBuffer> bufA;       size_t bufABytes   = 0;
    id<MTLBuffer> bufB;       size_t bufBBytes   = 0;
    id<MTLBuffer> bufC;       size_t bufCBytes   = 0;
    id<MTLBuffer> bufPartRe;  size_t partReCount = 0;  // current capacity in floats
    id<MTLBuffer> bufPartIm;  size_t partImCount = 0;
};

// ---------------------------------------------------------------------------
// Singleton creation
// ---------------------------------------------------------------------------
static MetalVectorContext* g_vecCtx  = nullptr;
static std::once_flag      g_vecOnce;

static id<MTLComputePipelineState> buildPS(id<MTLDevice> dev,
                                            id<MTLLibrary> lib,
                                            const char* name)
{
    NSError* err = nil;
    id<MTLFunction> fn = [lib newFunctionWithName:[NSString stringWithUTF8String:name]];
    if (!fn) {
        fprintf(stderr, "[MetalVectorOps] kernel not found: %s\n", name);
        return nil;
    }
    id<MTLComputePipelineState> ps =
        [dev newComputePipelineStateWithFunction:fn error:&err];
    if (!ps) {
        fprintf(stderr, "[MetalVectorOps] pipeline error for %s: %s\n",
                name, [[err localizedDescription] UTF8String]);
    }
    return ps;
}

static void createVecContext()
{
    id<MTLDevice> dev = MTLCreateSystemDefaultDevice();
    if (!dev) {
        fprintf(stderr, "[MetalVectorOps] No Metal device.\n");
        return;
    }

    // Load embedded metallib via dispatch_data_create (NOT NSData cast).
    dispatch_data_t dd = dispatch_data_create(
        vectorops_metallib, vectorops_metallib_len,
        NULL, DISPATCH_DATA_DESTRUCTOR_DEFAULT);
    NSError* err = nil;
    id<MTLLibrary> lib = [dev newLibraryWithData:dd error:&err];
    if (!lib) {
        fprintf(stderr, "[MetalVectorOps] metallib load failed: %s\n",
                [[err localizedDescription] UTF8String]);
        return;
    }

    auto* ctx = new MetalVectorContext();
    ctx->device  = dev;
    ctx->queue   = [dev newCommandQueue];
    ctx->library = lib;

    // Build all 16 pipeline states
    #define PS(field, kname) ctx->field = buildPS(dev, lib, kname)
    PS(ps_vec_add,        "pg_vec_add");
    PS(ps_vec_sub,        "pg_vec_sub");
    PS(ps_vec_mul,        "pg_vec_mul");
    PS(ps_vec_div,        "pg_vec_div");
    PS(ps_vec_add_scalar, "pg_vec_add_scalar");
    PS(ps_vec_mul_scalar, "pg_vec_mul_scalar");
    PS(ps_cvec_add,        "pg_cvec_add");
    PS(ps_cvec_sub,        "pg_cvec_sub");
    PS(ps_cvec_mul,        "pg_cvec_mul");
    PS(ps_cvec_div,        "pg_cvec_div");
    PS(ps_rvec_cmul,       "pg_rvec_cmul");
    PS(ps_cvec_axpy,       "pg_cvec_axpy");
    PS(ps_cvec_mul_scalar, "pg_cvec_mul_scalar");
    PS(ps_reduce_sum,      "pg_vec_reduce_sum");
    PS(ps_reduce_cdot,     "pg_cvec_reduce_cdot");
    PS(ps_reduce_norm2sq,  "pg_cvec_reduce_norm2sq");
    #undef PS

    g_vecCtx = ctx;
}

MetalVectorContext* metal_vector_get_context()
{
    std::call_once(g_vecOnce, createVecContext);
    return g_vecCtx;
}

// ---------------------------------------------------------------------------
// Buffer management — grow-only scratch buffers
// ---------------------------------------------------------------------------
static inline void growBuf(id<MTLDevice> dev,
                            id<MTLBuffer> __strong& buf,
                            size_t& curBytes,
                            size_t needed)
{
    if (needed > curBytes) {
        buf      = [dev newBufferWithLength:needed
                        options:MTLResourceStorageModeShared];
        curBytes = needed;
    }
}

static inline void growPartial(MetalVectorContext* ctx, size_t numGroups)
{
    size_t bytes = numGroups * sizeof(float);
    if (numGroups > ctx->partReCount) {
        ctx->bufPartRe  = [ctx->device newBufferWithLength:bytes
                                       options:MTLResourceStorageModeShared];
        ctx->partReCount = numGroups;
    }
    if (numGroups > ctx->partImCount) {
        ctx->bufPartIm  = [ctx->device newBufferWithLength:bytes
                                       options:MTLResourceStorageModeShared];
        ctx->partImCount = numGroups;
    }
}

// ---------------------------------------------------------------------------
// Generic element-wise dispatch: A op B → C  (3-buffer binding)
// bytesA/B/C may differ (e.g. rvec_cmul: A is real, B/C are complex).
// gridWidth = number of threads to launch (1 per element or complex pair).
// ---------------------------------------------------------------------------
static void ewise3(MetalVectorContext* ctx,
                    id<MTLComputePipelineState> ps,
                    const float* A, size_t bytesA,
                    const float* B, size_t bytesB,
                    float*       C, size_t bytesC,
                    size_t gridWidth)
{
    growBuf(ctx->device, ctx->bufA, ctx->bufABytes, bytesA);
    growBuf(ctx->device, ctx->bufB, ctx->bufBBytes, bytesB);
    growBuf(ctx->device, ctx->bufC, ctx->bufCBytes, bytesC);

    memcpy([ctx->bufA contents], A, bytesA);
    memcpy([ctx->bufB contents], B, bytesB);

    id<MTLCommandBuffer>         cmd = [ctx->queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:ps];
    [enc setBuffer:ctx->bufA offset:0 atIndex:0];
    [enc setBuffer:ctx->bufB offset:0 atIndex:1];
    [enc setBuffer:ctx->bufC offset:0 atIndex:2];

    NSUInteger tg = MIN((NSUInteger)ps.maxTotalThreadsPerThreadgroup, kTGSize);
    [enc dispatchThreads:MTLSizeMake(gridWidth, 1, 1)
         threadsPerThreadgroup:MTLSizeMake(tg, 1, 1)];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];

    memcpy(C, [ctx->bufC contents], bytesC);
}

// ===========================================================================
//  Element-wise REAL operations  (n = number of float elements)
// ===========================================================================

void metal_vec_add(MetalVectorContext* ctx,
                   const float* A, const float* B, float* C, size_t n)
{
    size_t bytes = n * sizeof(float);
    ewise3(ctx, ctx->ps_vec_add, A, bytes, B, bytes, C, bytes, n);
}

void metal_vec_sub(MetalVectorContext* ctx,
                   const float* A, const float* B, float* C, size_t n)
{
    size_t bytes = n * sizeof(float);
    ewise3(ctx, ctx->ps_vec_sub, A, bytes, B, bytes, C, bytes, n);
}

void metal_vec_mul(MetalVectorContext* ctx,
                   const float* A, const float* B, float* C, size_t n)
{
    size_t bytes = n * sizeof(float);
    ewise3(ctx, ctx->ps_vec_mul, A, bytes, B, bytes, C, bytes, n);
}

void metal_vec_div(MetalVectorContext* ctx,
                   const float* A, const float* B, float* C, size_t n)
{
    size_t bytes = n * sizeof(float);
    ewise3(ctx, ctx->ps_vec_div, A, bytes, B, bytes, C, bytes, n);
}

void metal_vec_add_scalar(MetalVectorContext* ctx,
                           const float* A, float scalar, float* C, size_t n)
{
    size_t bytes = n * sizeof(float);
    growBuf(ctx->device, ctx->bufA, ctx->bufABytes, bytes);
    growBuf(ctx->device, ctx->bufC, ctx->bufCBytes, bytes);
    memcpy([ctx->bufA contents], A, bytes);

    id<MTLCommandBuffer>         cmd = [ctx->queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:ctx->ps_vec_add_scalar];
    [enc setBuffer:ctx->bufA offset:0 atIndex:0];
    [enc setBytes:&scalar length:sizeof(float) atIndex:1];
    [enc setBuffer:ctx->bufC offset:0 atIndex:2];

    NSUInteger tg = MIN((NSUInteger)ctx->ps_vec_add_scalar.maxTotalThreadsPerThreadgroup, kTGSize);
    [enc dispatchThreads:MTLSizeMake(n, 1, 1)
         threadsPerThreadgroup:MTLSizeMake(tg, 1, 1)];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];

    memcpy(C, [ctx->bufC contents], bytes);
}

void metal_vec_mul_scalar(MetalVectorContext* ctx,
                           const float* A, float scalar, float* C, size_t n)
{
    size_t bytes = n * sizeof(float);
    growBuf(ctx->device, ctx->bufA, ctx->bufABytes, bytes);
    growBuf(ctx->device, ctx->bufC, ctx->bufCBytes, bytes);
    memcpy([ctx->bufA contents], A, bytes);

    id<MTLCommandBuffer>         cmd = [ctx->queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:ctx->ps_vec_mul_scalar];
    [enc setBuffer:ctx->bufA offset:0 atIndex:0];
    [enc setBytes:&scalar length:sizeof(float) atIndex:1];
    [enc setBuffer:ctx->bufC offset:0 atIndex:2];

    NSUInteger tg = MIN((NSUInteger)ctx->ps_vec_mul_scalar.maxTotalThreadsPerThreadgroup, kTGSize);
    [enc dispatchThreads:MTLSizeMake(n, 1, 1)
         threadsPerThreadgroup:MTLSizeMake(tg, 1, 1)];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];

    memcpy(C, [ctx->bufC contents], bytes);
}

// ===========================================================================
//  Element-wise COMPLEX operations  (n = complex element count; 2*n floats)
// ===========================================================================

void metal_cvec_add(MetalVectorContext* ctx,
                    const float* A, const float* B, float* C, size_t n)
{
    size_t bytes = 2 * n * sizeof(float);
    ewise3(ctx, ctx->ps_cvec_add, A, bytes, B, bytes, C, bytes, n);
}

void metal_cvec_sub(MetalVectorContext* ctx,
                    const float* A, const float* B, float* C, size_t n)
{
    size_t bytes = 2 * n * sizeof(float);
    ewise3(ctx, ctx->ps_cvec_sub, A, bytes, B, bytes, C, bytes, n);
}

void metal_cvec_mul(MetalVectorContext* ctx,
                    const float* A, const float* B, float* C, size_t n)
{
    size_t bytes = 2 * n * sizeof(float);
    ewise3(ctx, ctx->ps_cvec_mul, A, bytes, B, bytes, C, bytes, n);
}

void metal_cvec_div(MetalVectorContext* ctx,
                    const float* A, const float* B, float* C, size_t n)
{
    size_t bytes = 2 * n * sizeof(float);
    ewise3(ctx, ctx->ps_cvec_div, A, bytes, B, bytes, C, bytes, n);
}

void metal_rvec_cmul(MetalVectorContext* ctx,
                     const float* W, const float* X, float* C, size_t n)
{
    size_t realBytes = n * sizeof(float);
    size_t cplxBytes = 2 * n * sizeof(float);
    ewise3(ctx, ctx->ps_rvec_cmul, W, realBytes, X, cplxBytes, C, cplxBytes, n);
}

void metal_cvec_axpy(MetalVectorContext* ctx,
                     const float* A, const float* B, float* C,
                     float alphaRe, float alphaIm, size_t n)
{
    size_t bytes = 2 * n * sizeof(float);
    growBuf(ctx->device, ctx->bufA, ctx->bufABytes, bytes);
    growBuf(ctx->device, ctx->bufB, ctx->bufBBytes, bytes);
    growBuf(ctx->device, ctx->bufC, ctx->bufCBytes, bytes);

    memcpy([ctx->bufA contents], A, bytes);
    memcpy([ctx->bufB contents], B, bytes);

    float alpha[2] = {alphaRe, alphaIm};

    id<MTLCommandBuffer>         cmd = [ctx->queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:ctx->ps_cvec_axpy];
    [enc setBuffer:ctx->bufA offset:0 atIndex:0];
    [enc setBuffer:ctx->bufB offset:0 atIndex:1];
    [enc setBuffer:ctx->bufC offset:0 atIndex:2];
    [enc setBytes:alpha length:sizeof(float) * 2 atIndex:3];

    NSUInteger tg = MIN((NSUInteger)ctx->ps_cvec_axpy.maxTotalThreadsPerThreadgroup, kTGSize);
    [enc dispatchThreads:MTLSizeMake(n, 1, 1)
         threadsPerThreadgroup:MTLSizeMake(tg, 1, 1)];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];

    memcpy(C, [ctx->bufC contents], bytes);
}

void metal_cvec_mul_scalar(MetalVectorContext* ctx,
                           const float* A, float alphaRe, float alphaIm,
                           float* C, size_t n)
{
    size_t bytes = 2 * n * sizeof(float);
    growBuf(ctx->device, ctx->bufA, ctx->bufABytes, bytes);
    growBuf(ctx->device, ctx->bufC, ctx->bufCBytes, bytes);

    memcpy([ctx->bufA contents], A, bytes);

    float alpha[2] = {alphaRe, alphaIm};

    id<MTLCommandBuffer>         cmd = [ctx->queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:ctx->ps_cvec_mul_scalar];
    [enc setBuffer:ctx->bufA offset:0 atIndex:0];
    [enc setBytes:alpha length:sizeof(float) * 2 atIndex:1];
    [enc setBuffer:ctx->bufC offset:0 atIndex:2];

    NSUInteger tg = MIN((NSUInteger)ctx->ps_cvec_mul_scalar.maxTotalThreadsPerThreadgroup, kTGSize);
    [enc dispatchThreads:MTLSizeMake(n, 1, 1)
         threadsPerThreadgroup:MTLSizeMake(tg, 1, 1)];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];

    memcpy(C, [ctx->bufC contents], bytes);
}

// ===========================================================================
//  Reductions — GPU partial sums + CPU finish
//  Threadgroup size is fixed at kTGSize (256).
// ===========================================================================

float metal_vec_sum(MetalVectorContext* ctx, const float* A, size_t n)
{
    size_t bytesA    = n * sizeof(float);
    size_t numGroups = (n + kTGSize - 1) / kTGSize;

    growBuf(ctx->device, ctx->bufA, ctx->bufABytes, bytesA);
    growPartial(ctx, numGroups);
    memcpy([ctx->bufA contents], A, bytesA);

    uint32_t N = (uint32_t)n;

    id<MTLCommandBuffer>         cmd = [ctx->queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:ctx->ps_reduce_sum];
    [enc setBuffer:ctx->bufA      offset:0 atIndex:0];
    [enc setBuffer:ctx->bufPartRe offset:0 atIndex:1];
    [enc setBytes:&N length:sizeof(uint32_t) atIndex:2];

    [enc dispatchThreadgroups:MTLSizeMake(numGroups, 1, 1)
         threadsPerThreadgroup:MTLSizeMake(kTGSize, 1, 1)];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];

    const float* p = (const float*)[ctx->bufPartRe contents];
    float sum = 0.0f;
    for (size_t i = 0; i < numGroups; i++) sum += p[i];
    return sum;
}

void metal_cvec_cdot(MetalVectorContext* ctx,
                     const float* A, const float* B,
                     float* outRe, float* outIm, size_t n)
{
    size_t bytes     = 2 * n * sizeof(float);
    size_t numGroups = (n + kTGSize - 1) / kTGSize;

    growBuf(ctx->device, ctx->bufA, ctx->bufABytes, bytes);
    growBuf(ctx->device, ctx->bufB, ctx->bufBBytes, bytes);
    growPartial(ctx, numGroups);

    memcpy([ctx->bufA contents], A, bytes);
    memcpy([ctx->bufB contents], B, bytes);

    uint32_t N = (uint32_t)n;

    id<MTLCommandBuffer>         cmd = [ctx->queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:ctx->ps_reduce_cdot];
    [enc setBuffer:ctx->bufA      offset:0 atIndex:0];
    [enc setBuffer:ctx->bufB      offset:0 atIndex:1];
    [enc setBuffer:ctx->bufPartRe offset:0 atIndex:2];
    [enc setBuffer:ctx->bufPartIm offset:0 atIndex:3];
    [enc setBytes:&N length:sizeof(uint32_t) atIndex:4];

    [enc dispatchThreadgroups:MTLSizeMake(numGroups, 1, 1)
         threadsPerThreadgroup:MTLSizeMake(kTGSize, 1, 1)];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];

    const float* pre = (const float*)[ctx->bufPartRe contents];
    const float* pim = (const float*)[ctx->bufPartIm contents];
    float re = 0.0f, im = 0.0f;
    for (size_t i = 0; i < numGroups; i++) {
        re += pre[i];
        im += pim[i];
    }
    *outRe = re;
    *outIm = im;
}

float metal_cvec_norm2sq(MetalVectorContext* ctx, const float* A, size_t n)
{
    size_t bytes     = 2 * n * sizeof(float);
    size_t numGroups = (n + kTGSize - 1) / kTGSize;

    growBuf(ctx->device, ctx->bufA, ctx->bufABytes, bytes);
    growPartial(ctx, numGroups);
    memcpy([ctx->bufA contents], A, bytes);

    uint32_t N = (uint32_t)n;

    id<MTLCommandBuffer>         cmd = [ctx->queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:ctx->ps_reduce_norm2sq];
    [enc setBuffer:ctx->bufA      offset:0 atIndex:0];
    [enc setBuffer:ctx->bufPartRe offset:0 atIndex:1];
    [enc setBytes:&N length:sizeof(uint32_t) atIndex:2];

    [enc dispatchThreadgroups:MTLSizeMake(numGroups, 1, 1)
         threadsPerThreadgroup:MTLSizeMake(kTGSize, 1, 1)];
    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];

    const float* p = (const float*)[ctx->bufPartRe contents];
    float sum = 0.0f;
    for (size_t i = 0; i < numGroups; i++) sum += p[i];
    return sum;
}

#endif // METAL_COMPUTE
