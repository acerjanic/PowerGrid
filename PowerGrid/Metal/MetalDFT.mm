/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file MetalDFT.mm
/// @brief Objective-C++ bridge between the C++ DFT layer and Apple Metal.
///
/// All Objective-C types (id<MTLDevice>, etc.) are confined to this file.
/// The public API exposed to C++ callers is declared in MetalDFT.h.
///
/// The compiled Metal library (dft_metallib.h) is embedded as a C byte array
/// by the build system and included via the generated header.

#ifdef METAL_COMPUTE

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>

#include "MetalDFT.h"

// The xxd-embedded metallib
#include "dft_metallib.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>

// ---------------------------------------------------------------------------
// DFTParamsMSL — must match the Metal struct DFTParams in dft_metal.metal
// ---------------------------------------------------------------------------
struct DFTParamsMSL {
    uint32_t num_k;
    uint32_t num_i;
    uint32_t num_x;
    uint32_t num_y;
    uint32_t num_z;
};

// ---------------------------------------------------------------------------
// MetalDFTContext — opaque context struct
// ---------------------------------------------------------------------------
struct MetalDFTContext {
    id<MTLDevice>              device;
    id<MTLCommandQueue>        queue;
    id<MTLLibrary>             library;

    // Pipeline states
    id<MTLComputePipelineState> psForward;
    id<MTLComputePipelineState> psAdjoint;

    // Persistent GPU buffers for trajectory/coordinates (uploaded once)
    id<MTLBuffer> bufKx;
    id<MTLBuffer> bufKy;
    id<MTLBuffer> bufKz;
    id<MTLBuffer> bufIx;
    id<MTLBuffer> bufIy;
    id<MTLBuffer> bufIz;
    id<MTLBuffer> bufFM;
    id<MTLBuffer> bufT;

    // Optional gradient buffers (GdftR2 only)
    id<MTLBuffer> bufGx;
    id<MTLBuffer> bufGy;
    id<MTLBuffer> bufGz;
    bool hasGrads;

    // Working buffers for input/output (allocated once, reused)
    id<MTLBuffer> bufInR;
    id<MTLBuffer> bufInI;
    id<MTLBuffer> bufOutR;
    id<MTLBuffer> bufOutI;

    // Cached parameters
    DFTParamsMSL params;
};

// ---------------------------------------------------------------------------
// Helper: create pipeline state from kernel name
// ---------------------------------------------------------------------------
static id<MTLComputePipelineState> makePipeline(id<MTLDevice> dev,
                                                 id<MTLLibrary> lib,
                                                 const char* name) {
    NSError* err = nil;
    NSString* ns = [NSString stringWithUTF8String:name];
    id<MTLFunction> fn = [lib newFunctionWithName:ns];
    if (!fn) {
        fprintf(stderr, "[MetalDFT] kernel not found: %s\n", name);
        return nil;
    }
    id<MTLComputePipelineState> ps = [dev newComputePipelineStateWithFunction:fn error:&err];
    if (!ps) {
        fprintf(stderr, "[MetalDFT] pipeline creation failed for %s: %s\n",
                name, [[err localizedDescription] UTF8String]);
    }
    return ps;
}

// ---------------------------------------------------------------------------
// Helper: create shared buffer and upload host data
// ---------------------------------------------------------------------------
static id<MTLBuffer> uploadBuffer(id<MTLDevice> dev, const float* src, size_t count) {
    NSUInteger bytes = (NSUInteger)count * sizeof(float);
    id<MTLBuffer> buf = [dev newBufferWithLength:bytes options:MTLResourceStorageModeShared];
    if (buf && src) {
        memcpy([buf contents], src, bytes);
    }
    return buf;
}

// ---------------------------------------------------------------------------
// Common initialization
// ---------------------------------------------------------------------------
static MetalDFTContext* createCommon(
    const float* kx, const float* ky, const float* kz,
    const float* ix, const float* iy, const float* iz,
    const float* FM, const float* t,
    unsigned int num_k, unsigned int num_i,
    bool withGrads,
    const float* Gx, const float* Gy, const float* Gz,
    unsigned int num_x, unsigned int num_y, unsigned int num_z)
{
    id<MTLDevice> device = MTLCreateSystemDefaultDevice();
    if (!device) {
        fprintf(stderr, "[MetalDFT] No Metal device available.\n");
        return nullptr;
    }

    // Load pre-compiled metallib from embedded byte array
    NSError* err = nil;
    dispatch_data_t metallibDispData = dispatch_data_create(
        (const void*)dft_metallib,
        (size_t)dft_metallib_len,
        NULL,
        DISPATCH_DATA_DESTRUCTOR_DEFAULT);
    id<MTLLibrary> library = [device newLibraryWithData:metallibDispData error:&err];
    if (!library) {
        fprintf(stderr, "[MetalDFT] Failed to load metallib: %s\n",
                [[err localizedDescription] UTF8String]);
        return nullptr;
    }

    // Build pipeline states for the appropriate kernel variants
    const char* fwdName = withGrads ? "pg_dft_forward_grads" : "pg_dft_forward";
    const char* adjName = withGrads ? "pg_dft_adjoint_grads" : "pg_dft_adjoint";

    id<MTLComputePipelineState> psForward = makePipeline(device, library, fwdName);
    id<MTLComputePipelineState> psAdjoint = makePipeline(device, library, adjName);
    if (!psForward || !psAdjoint) {
        return nullptr;
    }

    id<MTLCommandQueue> queue = [device newCommandQueue];
    if (!queue) {
        fprintf(stderr, "[MetalDFT] Failed to create command queue.\n");
        return nullptr;
    }

    // Allocate context
    MetalDFTContext* ctx = new MetalDFTContext();
    ctx->device = device;
    ctx->queue = queue;
    ctx->library = library;
    ctx->psForward = psForward;
    ctx->psAdjoint = psAdjoint;

    // Upload constant arrays
    ctx->bufKx = uploadBuffer(device, kx, num_k);
    ctx->bufKy = uploadBuffer(device, ky, num_k);
    ctx->bufKz = uploadBuffer(device, kz, num_k);
    ctx->bufIx = uploadBuffer(device, ix, num_i);
    ctx->bufIy = uploadBuffer(device, iy, num_i);
    ctx->bufIz = uploadBuffer(device, iz, num_i);
    ctx->bufFM = uploadBuffer(device, FM, num_i);
    ctx->bufT  = uploadBuffer(device, t,  num_k);

    // Gradient buffers (GdftR2 only)
    ctx->hasGrads = withGrads;
    if (withGrads) {
        ctx->bufGx = uploadBuffer(device, Gx, num_i);
        ctx->bufGy = uploadBuffer(device, Gy, num_i);
        ctx->bufGz = uploadBuffer(device, Gz, num_i);
    }

    // Working buffers — max of num_k and num_i for input/output
    size_t maxN = (num_k > num_i) ? num_k : num_i;
    ctx->bufInR  = uploadBuffer(device, nullptr, maxN);
    ctx->bufInI  = uploadBuffer(device, nullptr, maxN);
    ctx->bufOutR = uploadBuffer(device, nullptr, maxN);
    ctx->bufOutI = uploadBuffer(device, nullptr, maxN);

    // Parameters
    ctx->params.num_k = num_k;
    ctx->params.num_i = num_i;
    ctx->params.num_x = num_x;
    ctx->params.num_y = num_y;
    ctx->params.num_z = num_z;

    return ctx;
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

MetalDFTContext* metal_dft_create(
    const float* kx, const float* ky, const float* kz,
    const float* ix, const float* iy, const float* iz,
    const float* FM, const float* t,
    unsigned int num_k, unsigned int num_i)
{
    return createCommon(kx, ky, kz, ix, iy, iz, FM, t,
                        num_k, num_i, false,
                        nullptr, nullptr, nullptr, 0, 0, 0);
}

MetalDFTContext* metal_dft_create_with_grads(
    const float* kx, const float* ky, const float* kz,
    const float* ix, const float* iy, const float* iz,
    const float* FM, const float* t,
    const float* Gx, const float* Gy, const float* Gz,
    unsigned int num_k, unsigned int num_i,
    unsigned int num_x, unsigned int num_y, unsigned int num_z)
{
    return createCommon(kx, ky, kz, ix, iy, iz, FM, t,
                        num_k, num_i, true,
                        Gx, Gy, Gz, num_x, num_y, num_z);
}

void metal_dft_forward(MetalDFTContext* ctx,
    const float* idata_r, const float* idata_i,
    float* kdata_r, float* kdata_i)
{
    if (!ctx) return;

    uint32_t num_k = ctx->params.num_k;
    uint32_t num_i = ctx->params.num_i;

    // Upload input data (image-space)
    memcpy([ctx->bufInR contents], idata_r, num_i * sizeof(float));
    memcpy([ctx->bufInI contents], idata_i, num_i * sizeof(float));

    id<MTLCommandBuffer> cmd = [ctx->queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:ctx->psForward];

    // Buffer bindings matching dft_metal.metal
    [enc setBuffer:ctx->bufInR  offset:0 atIndex:0];  // idata_r
    [enc setBuffer:ctx->bufInI  offset:0 atIndex:1];  // idata_i
    [enc setBuffer:ctx->bufOutR offset:0 atIndex:2];  // kdata_r
    [enc setBuffer:ctx->bufOutI offset:0 atIndex:3];  // kdata_i
    [enc setBuffer:ctx->bufKx   offset:0 atIndex:4];
    [enc setBuffer:ctx->bufKy   offset:0 atIndex:5];
    [enc setBuffer:ctx->bufKz   offset:0 atIndex:6];
    [enc setBuffer:ctx->bufIx   offset:0 atIndex:7];
    [enc setBuffer:ctx->bufIy   offset:0 atIndex:8];
    [enc setBuffer:ctx->bufIz   offset:0 atIndex:9];
    [enc setBuffer:ctx->bufFM   offset:0 atIndex:10];
    [enc setBuffer:ctx->bufT    offset:0 atIndex:11];
    [enc setBytes:&ctx->params length:sizeof(DFTParamsMSL) atIndex:12];

    if (ctx->hasGrads) {
        [enc setBuffer:ctx->bufGx offset:0 atIndex:13];
        [enc setBuffer:ctx->bufGy offset:0 atIndex:14];
        [enc setBuffer:ctx->bufGz offset:0 atIndex:15];
    }

    // Dispatch: one thread per k-space point
    NSUInteger tgSize = MIN(ctx->psForward.maxTotalThreadsPerThreadgroup, (NSUInteger)256);
    MTLSize grid = MTLSizeMake(num_k, 1, 1);
    MTLSize tg   = MTLSizeMake(tgSize, 1, 1);
    [enc dispatchThreads:grid threadsPerThreadgroup:tg];

    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];

    // Read back output
    memcpy(kdata_r, [ctx->bufOutR contents], num_k * sizeof(float));
    memcpy(kdata_i, [ctx->bufOutI contents], num_k * sizeof(float));
}

void metal_dft_adjoint(MetalDFTContext* ctx,
    const float* kdata_r, const float* kdata_i,
    float* idata_r, float* idata_i)
{
    if (!ctx) return;

    uint32_t num_k = ctx->params.num_k;
    uint32_t num_i = ctx->params.num_i;

    // Upload input data (k-space)
    memcpy([ctx->bufInR contents], kdata_r, num_k * sizeof(float));
    memcpy([ctx->bufInI contents], kdata_i, num_k * sizeof(float));

    id<MTLCommandBuffer> cmd = [ctx->queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:ctx->psAdjoint];

    // Buffer bindings matching dft_metal.metal
    [enc setBuffer:ctx->bufInR  offset:0 atIndex:0];  // kdata_r
    [enc setBuffer:ctx->bufInI  offset:0 atIndex:1];  // kdata_i
    [enc setBuffer:ctx->bufOutR offset:0 atIndex:2];  // idata_r
    [enc setBuffer:ctx->bufOutI offset:0 atIndex:3];  // idata_i
    [enc setBuffer:ctx->bufKx   offset:0 atIndex:4];
    [enc setBuffer:ctx->bufKy   offset:0 atIndex:5];
    [enc setBuffer:ctx->bufKz   offset:0 atIndex:6];
    [enc setBuffer:ctx->bufIx   offset:0 atIndex:7];
    [enc setBuffer:ctx->bufIy   offset:0 atIndex:8];
    [enc setBuffer:ctx->bufIz   offset:0 atIndex:9];
    [enc setBuffer:ctx->bufFM   offset:0 atIndex:10];
    [enc setBuffer:ctx->bufT    offset:0 atIndex:11];
    [enc setBytes:&ctx->params length:sizeof(DFTParamsMSL) atIndex:12];

    if (ctx->hasGrads) {
        [enc setBuffer:ctx->bufGx offset:0 atIndex:13];
        [enc setBuffer:ctx->bufGy offset:0 atIndex:14];
        [enc setBuffer:ctx->bufGz offset:0 atIndex:15];
    }

    // Dispatch: one thread per image pixel
    NSUInteger tgSize = MIN(ctx->psAdjoint.maxTotalThreadsPerThreadgroup, (NSUInteger)256);
    MTLSize grid = MTLSizeMake(num_i, 1, 1);
    MTLSize tg   = MTLSizeMake(tgSize, 1, 1);
    [enc dispatchThreads:grid threadsPerThreadgroup:tg];

    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];

    // Read back output
    memcpy(idata_r, [ctx->bufOutR contents], num_i * sizeof(float));
    memcpy(idata_i, [ctx->bufOutI contents], num_i * sizeof(float));
}

void metal_dft_destroy(MetalDFTContext* ctx) {
    if (ctx) {
        delete ctx;
    }
}

#endif // METAL_COMPUTE
