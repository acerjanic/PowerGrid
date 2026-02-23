/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file MetalNufftPipeline.h
/// @brief Full GPU NUFFT pipeline — pure C++ interface.
///
/// Replaces the split CPU+GPU approach (MetalGridding.h + CPU deapodize/
/// zero_pad/fftshift/FFT) with a single pipeline that keeps all intermediate
/// buffers on the GPU.  FFT uses MPSGraph (handles non-power-of-2 sizes).
///
/// Pipeline stages (all on GPU):
///   Forward (image → k-space):
///     deapodize → zero_pad → fftshift → FFT → ifftshift → gridding
///   Adjoint (k-space → image):
///     gridding → ifftshift → IFFT → fftshift → crop_center → deapodize

#pragma once
#ifdef METAL_COMPUTE

#include <cstddef>

/// @brief Opaque context holding Metal device, command queue, pipeline states,
///        MPSGraph FFT executables, and all GPU working buffers.
struct MetalNufftPipelineContext;

/// @brief Create a full GPU NUFFT pipeline context.
///
/// Uploads the Kaiser-Bessel LUT and k-space trajectory to Metal shared
/// buffers, allocates working buffers for all intermediate stages, and
/// pre-compiles MPSGraph FFT executables for the oversampled grid dimensions.
///
/// @param gridNx      Oversampled grid size in X (= ceil(gridOS * Nx), made even)
/// @param gridNy      Oversampled grid size in Y
/// @param gridNz      Oversampled grid size in Z (1 for 2D)
/// @param imageNx     Image size in X
/// @param imageNy     Image size in Y
/// @param imageNz     Image size in Z (1 for 2D)
/// @param gridOS      Grid oversampling factor
/// @param kernelWidth Kaiser-Bessel kernel width in grid cells
/// @param beta        Kaiser-Bessel beta parameter
/// @param LUT         Kaiser-Bessel lookup table (host pointer)
/// @param sizeLUT     Number of entries in @p LUT
/// @param kx          k-space X coordinates (host, numSamples elements)
/// @param ky          k-space Y coordinates (host, numSamples elements)
/// @param kz          k-space Z coordinates (host, numSamples elements)
/// @param numSamples  Number of k-space samples per coil
/// @returns Pointer to a new context, or nullptr on failure.
MetalNufftPipelineContext* metal_nufft_pipeline_create(
    int gridNx, int gridNy, int gridNz,
    int imageNx, int imageNy, int imageNz,
    float gridOS, float kernelWidth, float beta,
    const float* LUT, int sizeLUT,
    const float* kx, const float* ky, const float* kz,
    int numSamples);

/// @brief Destroy a pipeline context and release all GPU resources.
void metal_nufft_pipeline_destroy(MetalNufftPipelineContext* ctx);

/// @brief Full forward NUFFT: image → k-space (all on GPU).
///
/// @param ctx         Active pipeline context.
/// @param imageIn     Input image, interleaved complex (2 * Nx*Ny*Nz floats).
/// @param samplesOut  Output k-space, interleaved complex (2 * numSamples floats).
void metal_nufft_forward(MetalNufftPipelineContext* ctx,
                         const float* imageIn, float* samplesOut);

/// @brief Full adjoint NUFFT: k-space → image (all on GPU).
///
/// @param ctx         Active pipeline context.
/// @param samplesIn   Input k-space, interleaved complex (2 * numSamples floats).
/// @param imageOut    Output image, interleaved complex (2 * Nx*Ny*Nz floats).
void metal_nufft_adjoint(MetalNufftPipelineContext* ctx,
                         const float* samplesIn, float* imageOut);

// Dispatch statistics
uint64_t metal_nufft_dispatch_count();
double metal_nufft_wait_seconds();
void metal_nufft_reset_stats();

#endif // METAL_COMPUTE
