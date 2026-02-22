/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file MetalDFT.h
/// @brief Pure C++ interface to the Objective-C++ Metal DFT bridge.
///
/// The opaque `MetalDFTContext` struct is defined in MetalDFT.mm.
/// Callers create a context with metal_dft_create(), use it for forward
/// and adjoint DFT calls, then destroy it with metal_dft_destroy().

#pragma once
#ifdef METAL_COMPUTE

#include <cstddef>

struct MetalDFTContext;

/// @brief Create a Metal DFT context for Gdft (no gradient correction).
///
/// Uploads k-space trajectory, image coordinates, field map, and timing
/// vector to persistent Metal buffers.
///
/// @param kx,ky,kz  k-space coordinates (num_k elements each)
/// @param ix,iy,iz  image-space coordinates (num_i elements each)
/// @param FM        field map in radians/sec (num_i elements)
/// @param t         timing vector in seconds (num_k elements)
/// @param num_k     number of k-space samples
/// @param num_i     number of image pixels
/// @returns Pointer to new context, or nullptr on failure.
MetalDFTContext* metal_dft_create(
    const float* kx, const float* ky, const float* kz,
    const float* ix, const float* iy, const float* iz,
    const float* FM, const float* t,
    unsigned int num_k, unsigned int num_i);

/// @brief Create a Metal DFT context for GdftR2 (with gradient correction).
///
/// Same as metal_dft_create() but also uploads gradient maps and grid dims.
MetalDFTContext* metal_dft_create_with_grads(
    const float* kx, const float* ky, const float* kz,
    const float* ix, const float* iy, const float* iz,
    const float* FM, const float* t,
    const float* Gx, const float* Gy, const float* Gz,
    unsigned int num_k, unsigned int num_i,
    unsigned int num_x, unsigned int num_y, unsigned int num_z);

/// @brief Forward DFT: image → k-space.
///
/// @param ctx       Active Metal DFT context.
/// @param idata_r   Input image real part (num_i floats).
/// @param idata_i   Input image imag part (num_i floats).
/// @param kdata_r   Output k-space real part (num_k floats).
/// @param kdata_i   Output k-space imag part (num_k floats).
void metal_dft_forward(MetalDFTContext* ctx,
    const float* idata_r, const float* idata_i,
    float* kdata_r, float* kdata_i);

/// @brief Adjoint DFT: k-space → image.
///
/// @param ctx       Active Metal DFT context.
/// @param kdata_r   Input k-space real part (num_k floats).
/// @param kdata_i   Input k-space imag part (num_k floats).
/// @param idata_r   Output image real part (num_i floats).
/// @param idata_i   Output image imag part (num_i floats).
void metal_dft_adjoint(MetalDFTContext* ctx,
    const float* kdata_r, const float* kdata_i,
    float* idata_r, float* idata_i);

/// @brief Destroy a Metal DFT context and release resources.
void metal_dft_destroy(MetalDFTContext* ctx);

#endif // METAL_COMPUTE
