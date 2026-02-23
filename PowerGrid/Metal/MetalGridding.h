/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file MetalGridding.h
/// @brief Pure C++ interface to the Objective-C++ Metal gridding bridge.
///
/// All types in this header are plain C++ — no Objective-C or Metal API
/// types appear here, so this header can be included from regular C++ files.
///
/// The opaque `MetalGriddingContext` struct is defined in MetalGridding.mm.
/// Callers create a context with metal_gridding_create(), use it for adjoint
/// and forward gridding calls, then destroy it with metal_gridding_destroy().

#pragma once
#ifdef METAL_COMPUTE

#include <cstddef>

/// @brief Opaque context that holds Metal device, command queue, pipeline
///        states, and pre-uploaded buffers (LUT, kx/ky/kz).
struct MetalGriddingContext;

/// @brief Allocate and initialise a Metal gridding context.
///
/// Upload the Kaiser-Bessel LUT and k-space trajectory to Metal buffers once
/// here so they need not be transferred on every operator call.
///
/// @param gridNx      Oversampled grid size in X (= ceil(gridOS * Nx))
/// @param gridNy      Oversampled grid size in Y (= ceil(gridOS * Ny))
/// @param gridNz      Oversampled grid size in Z (1 for 2D)
/// @param imageNx     Image size in X
/// @param imageNy     Image size in Y
/// @param imageNz     Image size in Z (1 for 2D)
/// @param gridOS      Grid oversampling factor
/// @param kernelWidth Kaiser-Bessel kernel width in grid cells
/// @param LUT         Kaiser-Bessel lookup table (host pointer)
/// @param sizeLUT     Number of entries in @p LUT
/// @param kx          k-space coordinates in X (host pointer, @p numSamples elements)
/// @param ky          k-space coordinates in Y (host pointer, @p numSamples elements)
/// @param kz          k-space coordinates in Z (host pointer, @p numSamples elements; may be all-zero for 2D)
/// @param numSamples  Number of k-space samples per coil
/// @returns Pointer to a new context, or nullptr if Metal is unavailable
///          (e.g. device does not support MTLGPUFamilyApple6).
MetalGriddingContext* metal_gridding_create(
    int gridNx, int gridNy, int gridNz,
    int imageNx, int imageNy, int imageNz,
    float gridOS, float kernelWidth,
    const float* LUT, int sizeLUT,
    const float* kx, const float* ky, const float* kz,
    int numSamples);

/// @brief Destroy a Metal gridding context and release all Metal resources.
/// @param ctx  Context previously created by metal_gridding_create().
void metal_gridding_destroy(MetalGriddingContext* ctx);

/// @brief Adjoint (k-space → image) gridding for a 2D problem.
///
/// Zeroes the grid buffer, runs the adjoint scatter kernel on the GPU,
/// and leaves the result in @p pGridOut (interleaved real/imag floats,
/// gridNx * gridNy complex values).
///
/// @param ctx     Active Metal context.
/// @param dIn     Input k-space data, interleaved real/imag (2*numSamples floats).
/// @param pGridOut Output oversampled grid buffer (2 * gridNx * gridNy floats).
void metal_gridding_adjoint_2D(MetalGriddingContext* ctx,
                               const float* dIn, float* pGridOut);

/// @brief Adjoint (k-space → image) gridding for a 3D problem.
/// @param ctx     Active Metal context.
/// @param dIn     Input k-space data, interleaved real/imag (2*numSamples floats).
/// @param pGridOut Output oversampled grid buffer (2 * gridNx * gridNy * gridNz floats).
void metal_gridding_adjoint_3D(MetalGriddingContext* ctx,
                               const float* dIn, float* pGridOut);

/// @brief Forward (image → k-space) gridding for a 2D problem.
///
/// Reads from the oversampled grid @p pGridIn and accumulates into
/// @p pSamplesOut (one complex value per k-space sample).
///
/// @param ctx        Active Metal context.
/// @param pGridIn    Input oversampled grid (2 * gridNx * gridNy floats).
/// @param pSamplesOut Output k-space samples, interleaved real/imag (2*numSamples floats).
void metal_gridding_forward_2D(MetalGriddingContext* ctx,
                               const float* pGridIn, float* pSamplesOut);

/// @brief Forward (image → k-space) gridding for a 3D problem.
/// @param ctx        Active Metal context.
/// @param pGridIn    Input oversampled grid (2 * gridNx * gridNy * gridNz floats).
/// @param pSamplesOut Output k-space samples, interleaved real/imag (2*numSamples floats).
void metal_gridding_forward_3D(MetalGriddingContext* ctx,
                               const float* pGridIn, float* pSamplesOut);

// ============================================================================
// Dispatch statistics — counts command buffer commits and cumulative GPU wait.
// ============================================================================

uint64_t metal_gridding_dispatch_count();
double metal_gridding_wait_seconds();
void metal_gridding_reset_stats();

#endif // METAL_COMPUTE
