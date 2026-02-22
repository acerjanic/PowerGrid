/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file MetalNufftPipeline.mm
/// @brief Full GPU NUFFT pipeline using Metal compute + MPSGraph FFT.
///
/// Keeps all intermediate buffers on GPU.  Only two host↔GPU memcpy per call
/// (input and output).  FFT via MPSGraph supports non-power-of-2 sizes.

#ifdef METAL_COMPUTE

#import <Metal/Metal.h>
#import <Foundation/Foundation.h>
#import <MetalPerformanceShadersGraph/MetalPerformanceShadersGraph.h>

#include "MetalNufftPipeline.h"

// Embedded metallib byte arrays (generated at build time).
// gridding_metallib is already defined in MetalGridding.mm — use extern to avoid duplicate symbols.
extern const unsigned char gridding_metallib[];
extern const unsigned int gridding_metallib_len;
#include "nufft_support_metallib.h"

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cmath>

// ---------------------------------------------------------------------------
// MSL struct mirrors — must match Metal shader structs exactly
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

struct NufftSupportParamsMSL {
    int imageSize[3];
    int gridSize[3];
    float gridOS;
    float kernelWidth;
    float beta;
    int shift[3];
};

struct SampleFMSL {
    float real, imag, kX, kY, kZ, sdc, t, dummy;
};

// ---------------------------------------------------------------------------
// MetalNufftPipelineContext
// ---------------------------------------------------------------------------

struct MetalNufftPipelineContext {
    id<MTLDevice>              device;
    id<MTLCommandQueue>        queue;
    id<MTLLibrary>             griddingLib;
    id<MTLLibrary>             supportLib;

    // Gridding pipeline states
    id<MTLComputePipelineState> psZeroInit;
    id<MTLComputePipelineState> psAdjoint2D, psAdjoint3D;
    id<MTLComputePipelineState> psForward2D, psForward3D;

    // Support pipeline states
    id<MTLComputePipelineState> psDeapod2D, psDeapod3D;
    id<MTLComputePipelineState> psZeroPad2D, psZeroPad3D;
    id<MTLComputePipelineState> psCircshift2D, psCircshift3D;
    id<MTLComputePipelineState> psCrop2D, psCrop3D;

    // Pre-compiled MPSGraph FFT executables
    MPSGraphExecutable* fftFwdExec;
    MPSGraphExecutable* fftInvExec;

    // Persistent GPU buffers
    id<MTLBuffer> bufLUT, bufKx, bufKy, bufKz;

    // Working buffers
    id<MTLBuffer> bufImage;      // 2 * imageNumElems floats
    id<MTLBuffer> bufImageD;     // 2 * imageNumElems floats
    id<MTLBuffer> bufGridA;      // 2 * gridNumElems floats
    id<MTLBuffer> bufGridB;      // 2 * gridNumElems floats
    id<MTLBuffer> bufSamples;    // SampleFMSL[numSamples]
    id<MTLBuffer> bufSamplesOut; // 2 * numSamples floats

    // Dimensions
    int Nx, Ny, Nz, gNx, gNy, gNz;
    float gridOS, kernelWidth, beta;
    int numSamples, sizeLUT;
    int imageNumElems, gridNumElems;
    bool is3D;

    // Cached params
    GridParamsMSL griddingParams;
};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static id<MTLComputePipelineState> makePipeline(id<MTLDevice> dev,
                                                 id<MTLLibrary> lib,
                                                 const char* name) {
    NSError* err = nil;
    id<MTLFunction> fn = [lib newFunctionWithName:[NSString stringWithUTF8String:name]];
    if (!fn) { fprintf(stderr, "[NufftPipeline] kernel not found: %s\n", name); return nil; }
    id<MTLComputePipelineState> ps = [dev newComputePipelineStateWithFunction:fn error:&err];
    if (!ps) fprintf(stderr, "[NufftPipeline] pipeline failed: %s: %s\n",
                     name, [[err localizedDescription] UTF8String]);
    return ps;
}

static id<MTLLibrary> loadMetallib(id<MTLDevice> dev,
                                    const unsigned char* data, unsigned int len) {
    NSError* err = nil;
    dispatch_data_t dd = dispatch_data_create(data, (size_t)len, NULL,
                                               DISPATCH_DATA_DESTRUCTOR_DEFAULT);
    id<MTLLibrary> lib = [dev newLibraryWithData:dd error:&err];
    if (!lib) fprintf(stderr, "[NufftPipeline] metallib load failed: %s\n",
                      [[err localizedDescription] UTF8String]);
    return lib;
}

/// Encode one support kernel dispatch into an existing encoder.
static void encodeSupportKernel(id<MTLComputeCommandEncoder> enc,
                                id<MTLComputePipelineState> ps,
                                id<MTLBuffer> src, id<MTLBuffer> dst,
                                const NufftSupportParamsMSL& params,
                                NSUInteger numThreads) {
    [enc setComputePipelineState:ps];
    [enc setBuffer:src offset:0 atIndex:0];
    [enc setBuffer:dst offset:0 atIndex:1];
    [enc setBytes:&params length:sizeof(params) atIndex:2];
    NSUInteger tg = MIN((NSUInteger)ps.maxTotalThreadsPerThreadgroup, 256);
    [enc dispatchThreads:MTLSizeMake(numThreads, 1, 1)
   threadsPerThreadgroup:MTLSizeMake(tg, 1, 1)];
}

/// Insert a buffer-scope memory barrier.
static void encodeBarrier(id<MTLComputeCommandEncoder> enc) {
    [enc memoryBarrierWithScope:MTLBarrierScopeBuffers];
}

// ---------------------------------------------------------------------------
// MPSGraph FFT compilation (macOS 14.0+)
// ---------------------------------------------------------------------------

static MPSGraphExecutable* compileFFTExecutable(id<MTLDevice> device,
                                                 NSArray<NSNumber*>* shape,
                                                 BOOL inverse)
    API_AVAILABLE(macos(14.0))
{
    MPSGraph* graph = [[MPSGraph alloc] init];
    MPSGraphTensor* input = [graph placeholderWithShape:shape
                                              dataType:MPSDataTypeComplexFloat32
                                                  name:@"fft_in"];

    MPSGraphFFTDescriptor* desc = [MPSGraphFFTDescriptor descriptor];
    desc.inverse = inverse;
    desc.scalingMode = MPSGraphFFTScalingModeNone; // match FFTW unnormalized

    NSMutableArray<NSNumber*>* axes = [NSMutableArray array];
    for (NSUInteger i = 0; i < shape.count; i++)
        [axes addObject:@(i)];

    MPSGraphTensor* output = [graph fastFourierTransformWithTensor:input
                                                             axes:axes
                                                       descriptor:desc
                                                             name:@"fft_out"];

    MPSGraphDevice* mpsDevice = [MPSGraphDevice deviceWithMTLDevice:device];
    MPSGraphShapedType* inType = [[MPSGraphShapedType alloc]
        initWithShape:shape dataType:MPSDataTypeComplexFloat32];

    MPSGraphExecutable* exec = [graph compileWithDevice:mpsDevice
                                                  feeds:@{input: inType}
                                          targetTensors:@[output]
                                       targetOperations:nil
                                  compilationDescriptor:nil];
    return exec;
}

// ---------------------------------------------------------------------------
// metal_nufft_pipeline_create
// ---------------------------------------------------------------------------

MetalNufftPipelineContext* metal_nufft_pipeline_create(
    int gridNx, int gridNy, int gridNz,
    int imageNx, int imageNy, int imageNz,
    float gridOS, float kernelWidth, float beta,
    const float* LUT, int sizeLUT,
    const float* kx, const float* ky, const float* kz,
    int numSamples)
{
    // Require macOS 14.0+ for MPSGraph FFT
    if (@available(macOS 14.0, *)) {
        // OK
    } else {
        fprintf(stderr, "[NufftPipeline] Requires macOS 14.0+ for MPSGraph FFT.\n");
        return nullptr;
    }

    id<MTLDevice> device = MTLCreateSystemDefaultDevice();
    if (!device || ![device supportsFamily:MTLGPUFamilyApple6]) {
        fprintf(stderr, "[NufftPipeline] No suitable Metal device.\n");
        return nullptr;
    }

    // Load metallibs
    id<MTLLibrary> gLib = loadMetallib(device, gridding_metallib, gridding_metallib_len);
    id<MTLLibrary> sLib = loadMetallib(device, nufft_support_metallib, nufft_support_metallib_len);
    if (!gLib || !sLib) return nullptr;

    // Create pipeline states
    auto ps = [&](id<MTLLibrary> lib, const char* name) { return makePipeline(device, lib, name); };

    id<MTLComputePipelineState> psZeroInit  = ps(gLib, "pg_zero_init");
    id<MTLComputePipelineState> psAdj2D     = ps(gLib, "pg_gridding_adjoint_2D");
    id<MTLComputePipelineState> psAdj3D     = ps(gLib, "pg_gridding_adjoint_3D");
    id<MTLComputePipelineState> psFwd2D     = ps(gLib, "pg_gridding_forward_2D");
    id<MTLComputePipelineState> psFwd3D     = ps(gLib, "pg_gridding_forward_3D");
    id<MTLComputePipelineState> psDeapod2D  = ps(sLib, "pg_deapodize_2d");
    id<MTLComputePipelineState> psDeapod3D  = ps(sLib, "pg_deapodize_3d");
    id<MTLComputePipelineState> psZeroPad2D = ps(sLib, "pg_zero_pad_2d");
    id<MTLComputePipelineState> psZeroPad3D = ps(sLib, "pg_zero_pad_3d");
    id<MTLComputePipelineState> psCirc2D    = ps(sLib, "pg_circshift_2d");
    id<MTLComputePipelineState> psCirc3D    = ps(sLib, "pg_circshift_3d");
    id<MTLComputePipelineState> psCrop2D    = ps(sLib, "pg_crop_center_2d");
    id<MTLComputePipelineState> psCrop3D    = ps(sLib, "pg_crop_center_3d");

    if (!psZeroInit || !psAdj2D || !psAdj3D || !psFwd2D || !psFwd3D ||
        !psDeapod2D || !psDeapod3D || !psZeroPad2D || !psZeroPad3D ||
        !psCirc2D || !psCirc3D || !psCrop2D || !psCrop3D)
        return nullptr;

    id<MTLCommandQueue> queue = [device newCommandQueue];
    if (!queue) return nullptr;

    // Compile MPSGraph FFT executables
    bool is3D = (imageNz > 1);
    NSArray<NSNumber*>* fftShape = is3D
        ? @[@(gridNz), @(gridNy), @(gridNx)]
        : @[@(gridNy), @(gridNx)];

    MPSGraphExecutable* fftFwd = nil;
    MPSGraphExecutable* fftInv = nil;
    if (@available(macOS 14.0, *)) {
        fftFwd = compileFFTExecutable(device, fftShape, NO);
        fftInv = compileFFTExecutable(device, fftShape, YES);
    }
    if (!fftFwd || !fftInv) {
        fprintf(stderr, "[NufftPipeline] MPSGraph FFT compilation failed.\n");
        return nullptr;
    }

    // Allocate buffers
    int imageN = imageNx * imageNy * (is3D ? imageNz : 1);
    int gridN  = gridNx * gridNy * (is3D ? gridNz : 1);
    auto mkBuf = [&](NSUInteger bytes) {
        return [device newBufferWithLength:bytes options:MTLResourceStorageModeShared];
    };

    id<MTLBuffer> bufLUT  = mkBuf(sizeLUT * sizeof(float));
    id<MTLBuffer> bufKx   = mkBuf(numSamples * sizeof(float));
    id<MTLBuffer> bufKy   = mkBuf(numSamples * sizeof(float));
    id<MTLBuffer> bufKz   = mkBuf(numSamples * sizeof(float));
    id<MTLBuffer> bufImg  = mkBuf(2 * imageN * sizeof(float));
    id<MTLBuffer> bufImgD = mkBuf(2 * imageN * sizeof(float));
    id<MTLBuffer> bufGA   = mkBuf(2 * gridN * sizeof(float));
    id<MTLBuffer> bufGB   = mkBuf(2 * gridN * sizeof(float));
    id<MTLBuffer> bufSamp = mkBuf(numSamples * sizeof(SampleFMSL));
    id<MTLBuffer> bufSOut = mkBuf(2 * numSamples * sizeof(float));

    if (!bufLUT || !bufKx || !bufKy || !bufKz ||
        !bufImg || !bufImgD || !bufGA || !bufGB || !bufSamp || !bufSOut)
        return nullptr;

    // Upload persistent data
    memcpy([bufLUT contents], LUT, sizeLUT * sizeof(float));
    memcpy([bufKx contents], kx, numSamples * sizeof(float));
    memcpy([bufKy contents], ky, numSamples * sizeof(float));
    memcpy([bufKz contents], kz, numSamples * sizeof(float));

    // Build gridding params
    GridParamsMSL gp = {};
    gp.numSamples   = numSamples;
    gp.imageSize[0] = imageNx; gp.imageSize[1] = imageNy; gp.imageSize[2] = imageNz;
    gp.gridSize[0]  = gridNx;  gp.gridSize[1]  = gridNy;  gp.gridSize[2]  = gridNz;
    gp.gridOS       = gridOS;
    gp.kernelWidth  = kernelWidth;
    gp.sizeLUT      = sizeLUT;

    // Populate context
    MetalNufftPipelineContext* ctx = new MetalNufftPipelineContext();
    ctx->device = device;  ctx->queue = queue;
    ctx->griddingLib = gLib;  ctx->supportLib = sLib;
    ctx->psZeroInit = psZeroInit;
    ctx->psAdjoint2D = psAdj2D;  ctx->psAdjoint3D = psAdj3D;
    ctx->psForward2D = psFwd2D;  ctx->psForward3D = psFwd3D;
    ctx->psDeapod2D = psDeapod2D;  ctx->psDeapod3D = psDeapod3D;
    ctx->psZeroPad2D = psZeroPad2D;  ctx->psZeroPad3D = psZeroPad3D;
    ctx->psCircshift2D = psCirc2D;  ctx->psCircshift3D = psCirc3D;
    ctx->psCrop2D = psCrop2D;  ctx->psCrop3D = psCrop3D;
    ctx->fftFwdExec = fftFwd;  ctx->fftInvExec = fftInv;
    ctx->bufLUT = bufLUT;  ctx->bufKx = bufKx;  ctx->bufKy = bufKy;  ctx->bufKz = bufKz;
    ctx->bufImage = bufImg;  ctx->bufImageD = bufImgD;
    ctx->bufGridA = bufGA;  ctx->bufGridB = bufGB;
    ctx->bufSamples = bufSamp;  ctx->bufSamplesOut = bufSOut;
    ctx->Nx = imageNx;  ctx->Ny = imageNy;  ctx->Nz = imageNz;
    ctx->gNx = gridNx;  ctx->gNy = gridNy;  ctx->gNz = gridNz;
    ctx->gridOS = gridOS;  ctx->kernelWidth = kernelWidth;  ctx->beta = beta;
    ctx->numSamples = numSamples;  ctx->sizeLUT = sizeLUT;
    ctx->imageNumElems = imageN;  ctx->gridNumElems = gridN;
    ctx->is3D = is3D;
    ctx->griddingParams = gp;

    fprintf(stderr, "[NufftPipeline] Created: image=%dx%dx%d grid=%dx%dx%d samples=%d\n",
            imageNx, imageNy, imageNz, gridNx, gridNy, gridNz, numSamples);
    return ctx;
}

// ---------------------------------------------------------------------------
// metal_nufft_pipeline_destroy
// ---------------------------------------------------------------------------

void metal_nufft_pipeline_destroy(MetalNufftPipelineContext* ctx) {
    if (ctx) delete ctx; // ARC releases Obj-C members
}

// ---------------------------------------------------------------------------
// Helper: run MPSGraph FFT
// ---------------------------------------------------------------------------

static void runFFT(MetalNufftPipelineContext* ctx,
                   MPSGraphExecutable* exec,
                   id<MTLBuffer> inputBuf,
                   id<MTLBuffer> outputBuf)
    API_AVAILABLE(macos(14.0))
{
    NSArray<NSNumber*>* shape = ctx->is3D
        ? @[@(ctx->gNz), @(ctx->gNy), @(ctx->gNx)]
        : @[@(ctx->gNy), @(ctx->gNx)];

    MPSGraphTensorData* inData = [[MPSGraphTensorData alloc]
        initWithMTLBuffer:inputBuf shape:shape dataType:MPSDataTypeComplexFloat32];
    MPSGraphTensorData* outData = [[MPSGraphTensorData alloc]
        initWithMTLBuffer:outputBuf shape:shape dataType:MPSDataTypeComplexFloat32];

    [exec runWithMTLCommandQueue:ctx->queue
                     inputsArray:@[inData]
                    resultsArray:@[outData]
              executionDescriptor:nil];
}

// ---------------------------------------------------------------------------
// Helper: pack k-space samples into SampleFMSL struct array
// ---------------------------------------------------------------------------

static void packSamples(MetalNufftPipelineContext* ctx, const float* dIn) {
    SampleFMSL* dst = (SampleFMSL*)[ctx->bufSamples contents];
    const float* kxH = (const float*)[ctx->bufKx contents];
    const float* kyH = (const float*)[ctx->bufKy contents];
    const float* kzH = (const float*)[ctx->bufKz contents];
    for (int i = 0; i < ctx->numSamples; ++i) {
        dst[i] = { dIn[2*i], dIn[2*i+1], kxH[i], kyH[i], kzH[i], 1.0f, 0.0f, 0.0f };
    }
}

// ---------------------------------------------------------------------------
// Helper: build NufftSupportParamsMSL with specific shift values
// ---------------------------------------------------------------------------

static NufftSupportParamsMSL makeSupportParams(MetalNufftPipelineContext* ctx,
                                                int sx, int sy, int sz) {
    NufftSupportParamsMSL p = {};
    p.imageSize[0] = ctx->Nx;  p.imageSize[1] = ctx->Ny;  p.imageSize[2] = ctx->Nz;
    p.gridSize[0]  = ctx->gNx; p.gridSize[1]  = ctx->gNy; p.gridSize[2]  = ctx->gNz;
    p.gridOS = ctx->gridOS;  p.kernelWidth = ctx->kernelWidth;  p.beta = ctx->beta;
    p.shift[0] = sx;  p.shift[1] = sy;  p.shift[2] = sz;
    return p;
}

// ---------------------------------------------------------------------------
// metal_nufft_forward
//
// Pipeline: deapodize → zero_pad → fftshift → FFT → ifftshift → gridding
// ---------------------------------------------------------------------------

void metal_nufft_forward(MetalNufftPipelineContext* ctx,
                         const float* imageIn, float* samplesOut)
{
    if (@available(macOS 14.0, *)) {} else return;

    const bool is3D = ctx->is3D;

    // 1. Copy image to GPU shared buffer
    memcpy([ctx->bufImage contents], imageIn, 2 * ctx->imageNumElems * sizeof(float));

    // 2. Encode: deapodize → zero_pad → fftshift
    NufftSupportParamsMSL sp = makeSupportParams(ctx, 0, 0, 0);

    id<MTLCommandBuffer> cmd = [ctx->queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];

    // Deapodize: bufImage → bufImageD
    encodeSupportKernel(enc,
        is3D ? ctx->psDeapod3D : ctx->psDeapod2D,
        ctx->bufImage, ctx->bufImageD, sp,
        (NSUInteger)ctx->imageNumElems);
    encodeBarrier(enc);

    // Zero-pad: bufImageD → bufGridA
    encodeSupportKernel(enc,
        is3D ? ctx->psZeroPad3D : ctx->psZeroPad2D,
        ctx->bufImageD, ctx->bufGridA, sp,
        (NSUInteger)ctx->gridNumElems);
    encodeBarrier(enc);

    // fftshift: bufGridA → bufGridB
    // Shifts: floor(dim/2) for all dimensions
    sp.shift[0] = ctx->gNx / 2;
    sp.shift[1] = ctx->gNy / 2;
    sp.shift[2] = is3D ? ctx->gNz / 2 : 0;
    encodeSupportKernel(enc,
        is3D ? ctx->psCircshift3D : ctx->psCircshift2D,
        ctx->bufGridA, ctx->bufGridB, sp,
        (NSUInteger)ctx->gridNumElems);

    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];

    // 3. Forward FFT: bufGridB → bufGridA
    runFFT(ctx, ctx->fftFwdExec, ctx->bufGridB, ctx->bufGridA);

    // 4. Encode: ifftshift → forward gridding
    cmd = [ctx->queue commandBuffer];
    enc = [cmd computeCommandEncoder];

    // ifftshift: bufGridA → bufGridB
    // Shifts: ceil(dim/2) for x,y; floor for z (ifftshift3 asymmetry)
    sp.shift[0] = (ctx->gNx + 1) / 2;
    sp.shift[1] = (ctx->gNy + 1) / 2;
    sp.shift[2] = is3D ? ctx->gNz / 2 : 0;
    encodeSupportKernel(enc,
        is3D ? ctx->psCircshift3D : ctx->psCircshift2D,
        ctx->bufGridA, ctx->bufGridB, sp,
        (NSUInteger)ctx->gridNumElems);
    encodeBarrier(enc);

    // Forward gridding: bufGridB → bufSamplesOut
    auto& gp = ctx->griddingParams;
    [enc setComputePipelineState:is3D ? ctx->psForward3D : ctx->psForward2D];
    [enc setBuffer:ctx->bufGridB     offset:0 atIndex:0];
    [enc setBuffer:ctx->bufSamplesOut offset:0 atIndex:1];
    [enc setBuffer:ctx->bufLUT       offset:0 atIndex:2];
    [enc setBytes:&gp length:sizeof(gp) atIndex:3];
    [enc setBuffer:ctx->bufKx        offset:0 atIndex:4];
    [enc setBuffer:ctx->bufKy        offset:0 atIndex:5];
    if (is3D) [enc setBuffer:ctx->bufKz offset:0 atIndex:6];

    NSUInteger tg = MIN((NSUInteger)(is3D ? ctx->psForward3D : ctx->psForward2D).maxTotalThreadsPerThreadgroup, 256);
    [enc dispatchThreads:MTLSizeMake((NSUInteger)ctx->numSamples, 1, 1)
   threadsPerThreadgroup:MTLSizeMake(tg, 1, 1)];

    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];

    // 5. Copy result to host
    memcpy(samplesOut, [ctx->bufSamplesOut contents], 2 * ctx->numSamples * sizeof(float));
}

// ---------------------------------------------------------------------------
// metal_nufft_adjoint
//
// Pipeline: gridding → ifftshift → IFFT → fftshift → crop → deapodize
// ---------------------------------------------------------------------------

void metal_nufft_adjoint(MetalNufftPipelineContext* ctx,
                         const float* samplesIn, float* imageOut)
{
    if (@available(macOS 14.0, *)) {} else return;

    const bool is3D = ctx->is3D;

    // 1. Pack samples into SampleFMSL structs (CPU, shared memory)
    packSamples(ctx, samplesIn);

    // 2. Encode: zero grid → adjoint gridding → ifftshift
    id<MTLCommandBuffer> cmd = [ctx->queue commandBuffer];
    id<MTLComputeCommandEncoder> enc = [cmd computeCommandEncoder];

    // Zero the grid buffer (adjoint gridding uses atomics)
    NSUInteger gridFloats = (NSUInteger)ctx->gridNumElems * 2;
    [enc setComputePipelineState:ctx->psZeroInit];
    [enc setBuffer:ctx->bufGridA offset:0 atIndex:0];
    NSUInteger tg = MIN((NSUInteger)ctx->psZeroInit.maxTotalThreadsPerThreadgroup, 256);
    [enc dispatchThreads:MTLSizeMake(gridFloats, 1, 1)
   threadsPerThreadgroup:MTLSizeMake(tg, 1, 1)];
    encodeBarrier(enc);

    // Adjoint gridding: bufSamples → bufGridA (atomic scatter)
    auto& gp = ctx->griddingParams;
    [enc setComputePipelineState:is3D ? ctx->psAdjoint3D : ctx->psAdjoint2D];
    [enc setBuffer:ctx->bufSamples offset:0 atIndex:0];
    [enc setBuffer:ctx->bufGridA   offset:0 atIndex:1];
    [enc setBuffer:ctx->bufLUT     offset:0 atIndex:2];
    [enc setBytes:&gp length:sizeof(gp) atIndex:3];
    tg = MIN((NSUInteger)(is3D ? ctx->psAdjoint3D : ctx->psAdjoint2D).maxTotalThreadsPerThreadgroup, 256);
    [enc dispatchThreads:MTLSizeMake((NSUInteger)ctx->numSamples, 1, 1)
   threadsPerThreadgroup:MTLSizeMake(tg, 1, 1)];
    encodeBarrier(enc);

    // ifftshift: bufGridA → bufGridB
    NufftSupportParamsMSL sp = makeSupportParams(ctx, 0, 0, 0);
    if (is3D) {
        // ifftshift3: ceil(x/2), ceil(y/2), floor(z/2)
        sp.shift[0] = (ctx->gNx + 1) / 2;
        sp.shift[1] = (ctx->gNy + 1) / 2;
        sp.shift[2] = ctx->gNz / 2;
    } else {
        // ifftshift2: ceil(x/2), ceil(y/2)
        sp.shift[0] = (ctx->gNx + 1) / 2;
        sp.shift[1] = (ctx->gNy + 1) / 2;
    }
    encodeSupportKernel(enc,
        is3D ? ctx->psCircshift3D : ctx->psCircshift2D,
        ctx->bufGridA, ctx->bufGridB, sp,
        (NSUInteger)ctx->gridNumElems);

    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];

    // 3. Inverse FFT: bufGridB → bufGridA (unnormalized, matching FFTW_BACKWARD)
    runFFT(ctx, ctx->fftInvExec, ctx->bufGridB, ctx->bufGridA);

    // 4. Encode: fftshift → crop → deapodize
    cmd = [ctx->queue commandBuffer];
    enc = [cmd computeCommandEncoder];

    // fftshift: bufGridA → bufGridB
    sp.shift[0] = ctx->gNx / 2;
    sp.shift[1] = ctx->gNy / 2;
    sp.shift[2] = is3D ? ctx->gNz / 2 : 0;
    encodeSupportKernel(enc,
        is3D ? ctx->psCircshift3D : ctx->psCircshift2D,
        ctx->bufGridA, ctx->bufGridB, sp,
        (NSUInteger)ctx->gridNumElems);
    encodeBarrier(enc);

    // crop_center: bufGridB → bufImageD
    encodeSupportKernel(enc,
        is3D ? ctx->psCrop3D : ctx->psCrop2D,
        ctx->bufGridB, ctx->bufImageD, sp,
        (NSUInteger)ctx->imageNumElems);
    encodeBarrier(enc);

    // deapodize: bufImageD → bufImage
    encodeSupportKernel(enc,
        is3D ? ctx->psDeapod3D : ctx->psDeapod2D,
        ctx->bufImageD, ctx->bufImage, sp,
        (NSUInteger)ctx->imageNumElems);

    [enc endEncoding];
    [cmd commit];
    [cmd waitUntilCompleted];

    // 5. Copy result to host
    memcpy(imageOut, [ctx->bufImage contents], 2 * ctx->imageNumElems * sizeof(float));
}

#endif // METAL_COMPUTE
