/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file MetalGridding.mm
/// @brief Objective-C++ bridge between the C++ gridding layer and Apple Metal.
///
/// All Objective-C types (id<MTLDevice>, etc.) are confined to this file.
/// The public API exposed to C++ callers is declared in MetalGridding.h.
///
/// The compiled Metal library (gridding_metal.metallib) is embedded as a C
/// byte array by the build system (xxd -i) and included via the generated
/// header `gridding_metallib.h`.  This avoids any runtime file-path dependency.

#ifdef METAL_COMPUTE

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>

#include "MetalGridding.h"

// The xxd-embedded metallib: defines gridding_metallib[] and gridding_metallib_len
// Generated at build time into ${CMAKE_BINARY_DIR}/PowerGrid/Metal/
#include "gridding_metallib.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <atomic>
#include <mach/mach_time.h>

// ---------------------------------------------------------------------------
// Gridding dispatch statistics
// ---------------------------------------------------------------------------
static std::atomic<uint64_t> g_gridDispatchCount{0};
static std::atomic<uint64_t> g_gridWaitTicks{0};

static double gridTicksToSeconds(uint64_t ticks) {
    static mach_timebase_info_data_t tb = [] {
        mach_timebase_info_data_t info;
        mach_timebase_info(&info);
        return info;
    }();
    return (double)ticks * tb.numer / tb.denom / 1e9;
}

static void gridCommitAndWait(id<MTLCommandBuffer> cmd) {
    uint64_t t0 = mach_absolute_time();
    [cmd commit];
    [cmd waitUntilCompleted];
    uint64_t t1 = mach_absolute_time();
    g_gridDispatchCount.fetch_add(1, std::memory_order_relaxed);
    g_gridWaitTicks.fetch_add(t1 - t0, std::memory_order_relaxed);
}

uint64_t metal_gridding_dispatch_count() {
    return g_gridDispatchCount.load(std::memory_order_relaxed);
}

double metal_gridding_wait_seconds() {
    return gridTicksToSeconds(g_gridWaitTicks.load(std::memory_order_relaxed));
}

void metal_gridding_reset_stats() {
    g_gridDispatchCount.store(0, std::memory_order_relaxed);
    g_gridWaitTicks.store(0, std::memory_order_relaxed);
}

// ---------------------------------------------------------------------------
// GridParamsMSL — must match the Metal struct in gridding_metal.metal exactly
// ---------------------------------------------------------------------------
struct GridParamsMSL {
    int numSamples;
    int imageSize[3];
    int gridSize[3];
    float gridOS;
    float kernelWidth;
    int sizeLUT;
    int pad;
};

// ---------------------------------------------------------------------------
// SampleFMSL — matches ReconstructionSample<float> / SampleF in .metal
// ---------------------------------------------------------------------------
struct SampleFMSL {
    float real, imag, kX, kY, kZ, sdc, t, dummy;
};

// ---------------------------------------------------------------------------
// MetalGriddingContext — opaque context struct
// ---------------------------------------------------------------------------
struct MetalGriddingContext {
    // Metal objects (retain counted via ARC)
    id<MTLDevice>              device;
    id<MTLCommandQueue>        queue;
    id<MTLLibrary>             library;

    // Pipeline states for the five kernels
    id<MTLComputePipelineState> psZeroInit;
    id<MTLComputePipelineState> psAdjoint2D;
    id<MTLComputePipelineState> psAdjoint3D;
    id<MTLComputePipelineState> psForward2D;
    id<MTLComputePipelineState> psForward3D;

    // Persistent GPU buffers (MTLResourceStorageModeShared — zero-copy on AS)
    id<MTLBuffer> bufLUT;       ///< Kaiser-Bessel LUT (float[sizeLUT])
    id<MTLBuffer> bufKx;        ///< k-space x coords (float[numSamples])
    id<MTLBuffer> bufKy;        ///< k-space y coords (float[numSamples])
    id<MTLBuffer> bufKz;        ///< k-space z coords (float[numSamples])

    // Per-call working buffers: allocated once at creation, reused
    id<MTLBuffer> bufSamples;   ///< SampleFMSL[numSamples] (adjoint input)
    id<MTLBuffer> bufGrid;      ///< float[2 * gridNx * gridNy * gridNz]
    id<MTLBuffer> bufSamplesOut;///< float[2 * numSamples] (forward output)

    // Cached params
    GridParamsMSL params;
    int numSamples;
    int gridNumElems; ///< gNx * gNy * gNz
};

// ---------------------------------------------------------------------------
// Helper: compile a pipeline state for a named kernel function
// ---------------------------------------------------------------------------
static id<MTLComputePipelineState> makePipeline(id<MTLDevice> dev,
                                                 id<MTLLibrary> lib,
                                                 const char* name) {
    NSError* err = nil;
    NSString* ns = [NSString stringWithUTF8String:name];
    id<MTLFunction> fn = [lib newFunctionWithName:ns];
    if (!fn) {
        fprintf(stderr, "[MetalGridding] kernel not found: %s\n", name);
        return nil;
    }
    id<MTLComputePipelineState> ps = [dev newComputePipelineStateWithFunction:fn error:&err];
    if (!ps) {
        fprintf(stderr, "[MetalGridding] pipeline creation failed for %s: %s\n",
                name, [[err localizedDescription] UTF8String]);
    }
    return ps;
}

// ---------------------------------------------------------------------------
// metal_gridding_create
// ---------------------------------------------------------------------------
MetalGriddingContext* metal_gridding_create(
    int gridNx, int gridNy, int gridNz,
    int imageNx, int imageNy, int imageNz,
    float gridOS, float kernelWidth,
    const float* LUT, int sizeLUT,
    const float* kx, const float* ky, const float* kz,
    int numSamples)
{
    // Pick the default GPU
    id<MTLDevice> device = MTLCreateSystemDefaultDevice();
    if (!device) {
        fprintf(stderr, "[MetalGridding] No Metal device available.\n");
        return nullptr;
    }

    // Require Apple6+ family for atomic<float>
    if (![device supportsFamily:MTLGPUFamilyApple6]) {
        fprintf(stderr, "[MetalGridding] Device does not support MTLGPUFamilyApple6 "
                        "(atomic<float>); falling back to CPU.\n");
        return nullptr;
    }

    // Load the pre-compiled metallib from the embedded byte array.
    // dispatch_data_create (not NSData cast) is required — NSData does not
    // conform to the OS_dispatch_data protocol that Metal expects.
    NSError* err = nil;
    dispatch_data_t metallibDispData = dispatch_data_create(
        (const void*)gridding_metallib,
        (size_t)gridding_metallib_len,
        NULL,
        DISPATCH_DATA_DESTRUCTOR_DEFAULT);
    id<MTLLibrary> library = [device newLibraryWithData:metallibDispData error:&err];
    if (!library) {
        fprintf(stderr, "[MetalGridding] Failed to load metallib: %s\n",
                [[err localizedDescription] UTF8String]);
        return nullptr;
    }

    // Build pipeline states
    id<MTLComputePipelineState> psZeroInit  = makePipeline(device, library, "pg_zero_init");
    id<MTLComputePipelineState> psAdj2D    = makePipeline(device, library, "pg_gridding_adjoint_2D");
    id<MTLComputePipelineState> psAdj3D    = makePipeline(device, library, "pg_gridding_adjoint_3D");
    id<MTLComputePipelineState> psFwd2D    = makePipeline(device, library, "pg_gridding_forward_2D");
    id<MTLComputePipelineState> psFwd3D    = makePipeline(device, library, "pg_gridding_forward_3D");

    if (!psZeroInit || !psAdj2D || !psAdj3D || !psFwd2D || !psFwd3D) {
        return nullptr;
    }

    id<MTLCommandQueue> queue = [device newCommandQueue];
    if (!queue) {
        fprintf(stderr, "[MetalGridding] Failed to create command queue.\n");
        return nullptr;
    }

    // Allocate persistent shared-memory buffers
    const NSUInteger lutBytes     = (NSUInteger)sizeLUT   * sizeof(float);
    const NSUInteger kBytes       = (NSUInteger)numSamples * sizeof(float);
    const NSUInteger samplesBytes = (NSUInteger)numSamples * sizeof(SampleFMSL);
    const int gridNumElems        = gridNx * gridNy * gridNz;
    const NSUInteger gridBytes    = (NSUInteger)gridNumElems * 2 * sizeof(float);
    const NSUInteger outBytes     = (NSUInteger)numSamples * 2 * sizeof(float);

    auto mkBuf = [&](NSUInteger sz) -> id<MTLBuffer> {
        return [device newBufferWithLength:sz options:MTLResourceStorageModeShared];
    };

    id<MTLBuffer> bufLUT     = mkBuf(lutBytes);
    id<MTLBuffer> bufKx      = mkBuf(kBytes);
    id<MTLBuffer> bufKy      = mkBuf(kBytes);
    id<MTLBuffer> bufKz      = mkBuf(kBytes);
    id<MTLBuffer> bufSamples = mkBuf(samplesBytes);
    id<MTLBuffer> bufGrid    = mkBuf(gridBytes);
    id<MTLBuffer> bufSampOut = mkBuf(outBytes);

    if (!bufLUT || !bufKx || !bufKy || !bufKz || !bufSamples || !bufGrid || !bufSampOut) {
        fprintf(stderr, "[MetalGridding] Buffer allocation failed.\n");
        return nullptr;
    }

    // Upload LUT and k-space coords (persistent, uploaded once)
    memcpy([bufLUT contents], LUT, lutBytes);
    memcpy([bufKx  contents], kx,  kBytes);
    memcpy([bufKy  contents], ky,  kBytes);
    memcpy([bufKz  contents], kz,  kBytes);

    // Build the params struct
    GridParamsMSL params;
    params.numSamples   = numSamples;
    params.imageSize[0] = imageNx;
    params.imageSize[1] = imageNy;
    params.imageSize[2] = imageNz;
    params.gridSize[0]  = gridNx;
    params.gridSize[1]  = gridNy;
    params.gridSize[2]  = gridNz;
    params.gridOS       = gridOS;
    params.kernelWidth  = kernelWidth;
    params.sizeLUT      = sizeLUT;
    params.pad          = 0;

    // Allocate and populate context
    MetalGriddingContext* ctx = new MetalGriddingContext();
    ctx->device        = device;
    ctx->queue         = queue;
    ctx->library       = library;
    ctx->psZeroInit    = psZeroInit;
    ctx->psAdjoint2D   = psAdj2D;
    ctx->psAdjoint3D   = psAdj3D;
    ctx->psForward2D   = psFwd2D;
    ctx->psForward3D   = psFwd3D;
    ctx->bufLUT        = bufLUT;
    ctx->bufKx         = bufKx;
    ctx->bufKy         = bufKy;
    ctx->bufKz         = bufKz;
    ctx->bufSamples    = bufSamples;
    ctx->bufGrid       = bufGrid;
    ctx->bufSamplesOut = bufSampOut;
    ctx->params        = params;
    ctx->numSamples    = numSamples;
    ctx->gridNumElems  = gridNumElems;

    return ctx;
}

// ---------------------------------------------------------------------------
// metal_gridding_destroy
// ---------------------------------------------------------------------------
void metal_gridding_destroy(MetalGriddingContext* ctx) {
    if (ctx) {
        // ARC releases all Obj-C members when ctx is deleted
        delete ctx;
    }
}

// ---------------------------------------------------------------------------
// Helper: run pg_zero_init on a buffer
// ---------------------------------------------------------------------------
static void zeroBuffer(MetalGriddingContext* ctx, id<MTLBuffer> buf, NSUInteger numFloats) {
    id<MTLCommandBuffer>      cmd = [ctx->queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:ctx->psZeroInit];
    [enc setBuffer:buf offset:0 atIndex:0];
    NSUInteger tgSize = MIN((NSUInteger)ctx->psZeroInit.maxTotalThreadsPerThreadgroup, 256);
    MTLSize grid = MTLSizeMake(numFloats, 1, 1);
    MTLSize tg   = MTLSizeMake(tgSize, 1, 1);
    [enc dispatchThreads:grid threadsPerThreadgroup:tg];
    [enc endEncoding];
    gridCommitAndWait(cmd);
}

// ---------------------------------------------------------------------------
// Helper: pack dIn (interleaved float*) into SampleFMSL buffer, filling kxyz from context
// ---------------------------------------------------------------------------
static void packSamples(MetalGriddingContext* ctx, const float* dIn) {
    SampleFMSL* dst = (SampleFMSL*)[ctx->bufSamples contents];
    const float* kxHost = (const float*)[ctx->bufKx contents];
    const float* kyHost = (const float*)[ctx->bufKy contents];
    const float* kzHost = (const float*)[ctx->bufKz contents];
    for (int i = 0; i < ctx->numSamples; ++i) {
        dst[i].real  = dIn[2 * i];
        dst[i].imag  = dIn[2 * i + 1];
        dst[i].kX    = kxHost[i];
        dst[i].kY    = kyHost[i];
        dst[i].kZ    = kzHost[i];
        dst[i].sdc   = 1.0f;
        dst[i].t     = 0.0f;
        dst[i].dummy = 0.0f;
    }
}

// ---------------------------------------------------------------------------
// metal_gridding_adjoint_2D
// ---------------------------------------------------------------------------
void metal_gridding_adjoint_2D(MetalGriddingContext* ctx,
                               const float* dIn, float* pGridOut)
{
    const NSUInteger gridFloats = (NSUInteger)ctx->gridNumElems * 2;

    // Zero grid buffer on GPU
    zeroBuffer(ctx, ctx->bufGrid, gridFloats);

    // Pack samples on CPU (shared memory — no copy needed)
    packSamples(ctx, dIn);

    // Encode adjoint kernel
    id<MTLCommandBuffer>         cmd = [ctx->queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:ctx->psAdjoint2D];
    [enc setBuffer:ctx->bufSamples offset:0 atIndex:0];
    [enc setBuffer:ctx->bufGrid    offset:0 atIndex:1];
    [enc setBuffer:ctx->bufLUT     offset:0 atIndex:2];
    [enc setBytes:&ctx->params length:sizeof(GridParamsMSL) atIndex:3];

    NSUInteger tgSize = MIN((NSUInteger)ctx->psAdjoint2D.maxTotalThreadsPerThreadgroup, 256);
    MTLSize grid = MTLSizeMake((NSUInteger)ctx->numSamples, 1, 1);
    MTLSize tg   = MTLSizeMake(tgSize, 1, 1);
    [enc dispatchThreads:grid threadsPerThreadgroup:tg];
    [enc endEncoding];
    gridCommitAndWait(cmd);

    // Copy result (shared mem — already coherent on Apple Silicon)
    memcpy(pGridOut, [ctx->bufGrid contents], gridFloats * sizeof(float));
}

// ---------------------------------------------------------------------------
// metal_gridding_adjoint_3D
// ---------------------------------------------------------------------------
void metal_gridding_adjoint_3D(MetalGriddingContext* ctx,
                               const float* dIn, float* pGridOut)
{
    const NSUInteger gridFloats = (NSUInteger)ctx->gridNumElems * 2;

    zeroBuffer(ctx, ctx->bufGrid, gridFloats);
    packSamples(ctx, dIn);

    id<MTLCommandBuffer>         cmd = [ctx->queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:ctx->psAdjoint3D];
    [enc setBuffer:ctx->bufSamples offset:0 atIndex:0];
    [enc setBuffer:ctx->bufGrid    offset:0 atIndex:1];
    [enc setBuffer:ctx->bufLUT     offset:0 atIndex:2];
    [enc setBytes:&ctx->params length:sizeof(GridParamsMSL) atIndex:3];

    NSUInteger tgSize = MIN((NSUInteger)ctx->psAdjoint3D.maxTotalThreadsPerThreadgroup, 256);
    MTLSize grid = MTLSizeMake((NSUInteger)ctx->numSamples, 1, 1);
    MTLSize tg   = MTLSizeMake(tgSize, 1, 1);
    [enc dispatchThreads:grid threadsPerThreadgroup:tg];
    [enc endEncoding];
    gridCommitAndWait(cmd);

    memcpy(pGridOut, [ctx->bufGrid contents], gridFloats * sizeof(float));
}

// ---------------------------------------------------------------------------
// Helper: run forward kernel and copy results
// ---------------------------------------------------------------------------
static void runForward(MetalGriddingContext* ctx,
                       id<MTLComputePipelineState> ps,
                       const float* pGridIn,
                       float* pSamplesOut,
                       bool is3D)
{
    const NSUInteger gridFloats = (NSUInteger)ctx->gridNumElems * 2;
    const NSUInteger sampFloats = (NSUInteger)ctx->numSamples   * 2;

    // Copy grid data into shared buffer
    memcpy([ctx->bufGrid    contents], pGridIn, gridFloats * sizeof(float));

    // Zero the output samples buffer
    zeroBuffer(ctx, ctx->bufSamplesOut, sampFloats);

    id<MTLCommandBuffer>         cmd = [ctx->queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];
    [enc setComputePipelineState:ps];
    [enc setBuffer:ctx->bufGrid       offset:0 atIndex:0];
    [enc setBuffer:ctx->bufSamplesOut offset:0 atIndex:1];
    [enc setBuffer:ctx->bufLUT        offset:0 atIndex:2];
    [enc setBytes:&ctx->params length:sizeof(GridParamsMSL) atIndex:3];
    [enc setBuffer:ctx->bufKx         offset:0 atIndex:4];
    [enc setBuffer:ctx->bufKy         offset:0 atIndex:5];
    if (is3D) {
        [enc setBuffer:ctx->bufKz     offset:0 atIndex:6];
    }

    NSUInteger tgSize = MIN((NSUInteger)ps.maxTotalThreadsPerThreadgroup, 256);
    MTLSize grid = MTLSizeMake((NSUInteger)ctx->numSamples, 1, 1);
    MTLSize tg   = MTLSizeMake(tgSize, 1, 1);
    [enc dispatchThreads:grid threadsPerThreadgroup:tg];
    [enc endEncoding];
    gridCommitAndWait(cmd);

    memcpy(pSamplesOut, [ctx->bufSamplesOut contents], sampFloats * sizeof(float));
}

// ---------------------------------------------------------------------------
// metal_gridding_forward_2D / _3D
// ---------------------------------------------------------------------------
void metal_gridding_forward_2D(MetalGriddingContext* ctx,
                               const float* pGridIn, float* pSamplesOut)
{
    runForward(ctx, ctx->psForward2D, pGridIn, pSamplesOut, false);
}

void metal_gridding_forward_3D(MetalGriddingContext* ctx,
                               const float* pGridIn, float* pSamplesOut)
{
    runForward(ctx, ctx->psForward3D, pGridIn, pSamplesOut, true);
}

#endif // METAL_COMPUTE
