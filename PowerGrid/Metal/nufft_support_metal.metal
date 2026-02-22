/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file nufft_support_metal.metal
/// @brief Metal compute shaders for NUFFT support operations.
///
/// 8 compute kernels for the operations surrounding the core gridding/FFT:
///   pg_deapodize_2d/3d    — Kaiser-Bessel deapodization
///   pg_zero_pad_2d/3d     — Zero-pad image into oversampled grid (fused)
///   pg_circshift_2d/3d    — Circular shift (for fftshift/ifftshift)
///   pg_crop_center_2d/3d  — Extract center region from oversampled grid
///
/// All data uses interleaved complex format: [re0, im0, re1, im1, ...]
/// Index conventions exactly replicate griddingSupport.cpp.

#include <metal_stdlib>
using namespace metal;

// ---------------------------------------------------------------------------
// Shared parameter block — host fills whichever fields the kernel needs
// ---------------------------------------------------------------------------

struct NufftSupportParams {
    int imageSize[3];   ///< [Nx, Ny, Nz]  — original image dimensions
    int gridSize[3];    ///< [gNx, gNy, gNz] — oversampled grid dimensions
    float gridOS;       ///< grid oversampling factor
    float kernelWidth;  ///< Kaiser-Bessel kernel width
    float beta;         ///< Kaiser-Bessel beta parameter
    int shift[3];       ///< [xshift, yshift, zshift] for circshift
};

// =========================================================================
// pg_deapodize_2d — 1 thread per image pixel
// Index: Y + X * imageY  (column-major, matching CPU deapodization2d)
// Dispatch: imageX * imageY threads
// =========================================================================

kernel void pg_deapodize_2d(
    device const float*           pSrc   [[buffer(0)]],
    device float*                 pDst   [[buffer(1)]],
    constant NufftSupportParams&  params [[buffer(2)]],
    uint gid [[thread_position_in_grid]])
{
    const int imageX = params.imageSize[0];
    const int imageY = params.imageSize[1];
    if ((int)gid >= imageX * imageY) return;

    // Column-major decomposition: index = Y + X * imageY
    const int X = (int)gid / imageY;
    const int Y = (int)gid % imageY;

    // CPU 2D uses integer division: imageY/2, imageX/2
    const float gkY = float(Y - imageY / 2) / float(imageY);
    const float gkX = float(X - imageX / 2) / float(imageX);

    const float pi2 = M_PI_F * M_PI_F;
    const float kw  = params.kernelWidth;
    const float kw2 = kw * kw;
    const float b2  = params.beta * params.beta;

    const float exprX = pi2 * kw2 * gkX * gkX - b2;
    const float exprY = pi2 * kw2 * gkY * gkY - b2;

    float vX, vY;
    if (exprX >= 0.0f) { float s = sqrt(exprX); vX = sin(s) / s; }
    else               { float s = sqrt(-exprX); vX = sinh(s) / s; }
    if (exprY >= 0.0f) { float s = sqrt(exprY); vY = sin(s) / s; }
    else               { float s = sqrt(-exprY); vY = sinh(s) / s; }

    const float gk  = vX * vY;
    const int   idx = 2 * (int)gid;

    if (isnan(gk)) {
        pDst[idx]     = 0.0f;
        pDst[idx + 1] = 0.0f;
    } else {
        const float scale = 1.0f / (gk * params.gridOS * params.gridOS);
        pDst[idx]     = pSrc[idx]     * scale;
        pDst[idx + 1] = pSrc[idx + 1] * scale;
    }
}

// =========================================================================
// pg_deapodize_3d — 1 thread per image voxel
// Index: Z*imageX*imageY + X*imageY + Y  (matching CPU deapodization3d)
// Dispatch: imageX * imageY * imageZ threads
// =========================================================================

kernel void pg_deapodize_3d(
    device const float*           pSrc   [[buffer(0)]],
    device float*                 pDst   [[buffer(1)]],
    constant NufftSupportParams&  params [[buffer(2)]],
    uint gid [[thread_position_in_grid]])
{
    const int imageX = params.imageSize[0];
    const int imageY = params.imageSize[1];
    const int imageZ = params.imageSize[2];
    if ((int)gid >= imageX * imageY * imageZ) return;

    // Decompose: index = Z*imageX*imageY + X*imageY + Y
    const int Y = (int)gid % imageY;
    const int X = ((int)gid / imageY) % imageX;
    const int Z = (int)gid / (imageX * imageY);

    // CPU 3D uses float division: (T1)imageZ / 2.0
    const float gkZ = (float(Z) - float(imageZ) / 2.0f) / float(imageZ);
    const float gkY = (float(Y) - float(imageY) / 2.0f) / float(imageY);
    const float gkX = (float(X) - float(imageX) / 2.0f) / float(imageX);

    const float pi2 = M_PI_F * M_PI_F;
    const float kw  = params.kernelWidth;
    const float kw2 = kw * kw;
    const float b2  = params.beta * params.beta;

    const float exprX = pi2 * kw2 * gkX * gkX - b2;
    const float exprY = pi2 * kw2 * gkY * gkY - b2;
    const float exprZ = pi2 * kw2 * gkZ * gkZ - b2;

    float vX, vY, vZ;
    if (exprX >= 0.0f) { float s = sqrt(exprX); vX = sin(s) / s; }
    else               { float s = sqrt(-exprX); vX = sinh(s) / s; }
    if (exprY >= 0.0f) { float s = sqrt(exprY); vY = sin(s) / s; }
    else               { float s = sqrt(-exprY); vY = sinh(s) / s; }
    if (exprZ >= 0.0f) { float s = sqrt(exprZ); vZ = sin(s) / s; }
    else               { float s = sqrt(-exprZ); vZ = sinh(s) / s; }

    const float gk  = vX * vY * vZ;
    const int   idx = 2 * (int)gid;

    if (isnan(gk)) {
        pDst[idx]     = 0.0f;
        pDst[idx + 1] = 0.0f;
    } else {
        const float scale = 1.0f / (gk * params.gridOS * params.gridOS * params.gridOS);
        pDst[idx]     = pSrc[idx]     * scale;
        pDst[idx + 1] = pSrc[idx + 1] * scale;
    }
}

// =========================================================================
// pg_zero_pad_2d — fused zero-init + copy, 1 thread per dest element
// Dest index: dY_dst * gridSizeX + dX_dst  (row-major, matching CPU)
// Dispatch: gridSizeX * gridSizeY threads
// =========================================================================

kernel void pg_zero_pad_2d(
    device const float*           pSrc   [[buffer(0)]],
    device float*                 pDst   [[buffer(1)]],
    constant NufftSupportParams&  params [[buffer(2)]],
    uint gid [[thread_position_in_grid]])
{
    const int gNx    = params.gridSize[0];
    const int gNy    = params.gridSize[1];
    const int imageX = params.imageSize[0];
    const int imageY = params.imageSize[1];
    if ((int)gid >= gNx * gNy) return;

    // Row-major decomposition: dst_idx = dY_dst * gNx + dX_dst
    const int dX_dst = (int)gid % gNx;
    const int dY_dst = (int)gid / gNx;

    const int offsetX = (gNx - imageX) / 2;
    const int offsetY = (gNy - imageY) / 2;

    const int dX_src = dX_dst - offsetX;
    const int dY_src = dY_dst - offsetY;

    const int idx = 2 * (int)gid;
    if (dX_src >= 0 && dX_src < imageX && dY_src >= 0 && dY_src < imageY) {
        // CPU src index: dY_src * imageSizeX + dX_src
        const int src_idx = dY_src * imageX + dX_src;
        pDst[idx]     = pSrc[2 * src_idx];
        pDst[idx + 1] = pSrc[2 * src_idx + 1];
    } else {
        pDst[idx]     = 0.0f;
        pDst[idx + 1] = 0.0f;
    }
}

// =========================================================================
// pg_zero_pad_3d — fused zero-init + copy, 1 thread per dest element
// Dest index: dZ_dst*gNx*gNy + dX_dst*gNy + dY_dst  (matching CPU)
// Dispatch: gridSizeX * gridSizeY * gridSizeZ threads
// =========================================================================

kernel void pg_zero_pad_3d(
    device const float*           pSrc   [[buffer(0)]],
    device float*                 pDst   [[buffer(1)]],
    constant NufftSupportParams&  params [[buffer(2)]],
    uint gid [[thread_position_in_grid]])
{
    const int gNx    = params.gridSize[0];
    const int gNy    = params.gridSize[1];
    const int gNz    = params.gridSize[2];
    const int imageX = params.imageSize[0];
    const int imageY = params.imageSize[1];
    const int imageZ = params.imageSize[2];
    if ((int)gid >= gNx * gNy * gNz) return;

    // Decompose: dst_idx = dZ*gNx*gNy + dX*gNy + dY
    const int dY_dst = (int)gid % gNy;
    const int dX_dst = ((int)gid / gNy) % gNx;
    const int dZ_dst = (int)gid / (gNx * gNy);

    const int offsetX = (gNx - imageX) / 2;
    const int offsetY = (gNy - imageY) / 2;
    const int offsetZ = (gNz - imageZ) / 2;

    const int dX_src = dX_dst - offsetX;
    const int dY_src = dY_dst - offsetY;
    const int dZ_src = dZ_dst - offsetZ;

    const int idx = 2 * (int)gid;
    if (dX_src >= 0 && dX_src < imageX &&
        dY_src >= 0 && dY_src < imageY &&
        dZ_src >= 0 && dZ_src < imageZ) {
        // CPU src index: dZ*imageX*imageY + dX*imageY + dY
        const int src_idx = dZ_src * imageX * imageY + dX_src * imageY + dY_src;
        pDst[idx]     = pSrc[2 * src_idx];
        pDst[idx + 1] = pSrc[2 * src_idx + 1];
    } else {
        pDst[idx]     = 0.0f;
        pDst[idx + 1] = 0.0f;
    }
}

// =========================================================================
// pg_circshift_2d — 1 thread per element (source-side)
// Index: x + y * xdim  (matching CPU circshift2, xdim = gridSize[0])
// Dispatch: gridSize[0] * gridSize[1] threads
// =========================================================================

kernel void pg_circshift_2d(
    device const float*           pSrc   [[buffer(0)]],
    device float*                 pDst   [[buffer(1)]],
    constant NufftSupportParams&  params [[buffer(2)]],
    uint gid [[thread_position_in_grid]])
{
    const int xdim = params.gridSize[0];
    const int ydim = params.gridSize[1];
    if ((int)gid >= xdim * ydim) return;

    // Decompose: src_idx = x + y * xdim
    const int x = (int)gid % xdim;
    const int y = (int)gid / xdim;

    const int ii = (x + params.shift[0]) % xdim;
    const int jj = (y + params.shift[1]) % ydim;

    const int dst_idx = ii + jj * xdim;
    pDst[2 * dst_idx]     = pSrc[2 * (int)gid];
    pDst[2 * dst_idx + 1] = pSrc[2 * (int)gid + 1];
}

// =========================================================================
// pg_circshift_3d — 1 thread per element (source-side)
// Index: x + y*xdim + z*xdim*ydim  (matching CPU circshift3)
// Dispatch: gridSize[0] * gridSize[1] * gridSize[2] threads
// =========================================================================

kernel void pg_circshift_3d(
    device const float*           pSrc   [[buffer(0)]],
    device float*                 pDst   [[buffer(1)]],
    constant NufftSupportParams&  params [[buffer(2)]],
    uint gid [[thread_position_in_grid]])
{
    const int xdim = params.gridSize[0];
    const int ydim = params.gridSize[1];
    const int zdim = params.gridSize[2];
    if ((int)gid >= xdim * ydim * zdim) return;

    // Decompose: src_idx = x + y*xdim + z*xdim*ydim
    const int x = (int)gid % xdim;
    const int y = ((int)gid / xdim) % ydim;
    const int z = (int)gid / (xdim * ydim);

    const int ii = (x + params.shift[0]) % xdim;
    const int jj = (y + params.shift[1]) % ydim;
    const int kk = (z + params.shift[2]) % zdim;

    const int dst_idx = ii + jj * xdim + kk * xdim * ydim;
    pDst[2 * dst_idx]     = pSrc[2 * (int)gid];
    pDst[2 * dst_idx + 1] = pSrc[2 * (int)gid + 1];
}

// =========================================================================
// pg_crop_center_2d — 1 thread per destination pixel
// Dst index: dX_dst * imageSizeY + dY_dst  (column-major, matching CPU)
// Dispatch: imageSizeX * imageSizeY threads
// =========================================================================

kernel void pg_crop_center_2d(
    device const float*           pSrc   [[buffer(0)]],
    device float*                 pDst   [[buffer(1)]],
    constant NufftSupportParams&  params [[buffer(2)]],
    uint gid [[thread_position_in_grid]])
{
    const int imageX = params.imageSize[0];
    const int imageY = params.imageSize[1];
    const int gNx    = params.gridSize[0];
    const int gNy    = params.gridSize[1];
    if ((int)gid >= imageX * imageY) return;

    // Column-major decomposition: dst_idx = dX_dst * imageY + dY_dst
    const int dY_dst = (int)gid % imageY;
    const int dX_dst = (int)gid / imageY;

    const int offsetY = (gNy - imageY) / 2;
    const int offsetX = (gNx - imageX) / 2;

    const int dY_src = dY_dst + offsetY;
    const int dX_src = dX_dst + offsetX;

    // CPU src index: dX_src * gridSizeY + dY_src
    const int src_idx = dX_src * gNy + dY_src;
    pDst[2 * (int)gid]     = pSrc[2 * src_idx];
    pDst[2 * (int)gid + 1] = pSrc[2 * src_idx + 1];
}

// =========================================================================
// pg_crop_center_3d — 1 thread per destination voxel
// Dst index: dZ*imageX*imageY + dX*imageY + dY  (matching CPU)
// Dispatch: imageSizeX * imageSizeY * imageSizeZ threads
// =========================================================================

kernel void pg_crop_center_3d(
    device const float*           pSrc   [[buffer(0)]],
    device float*                 pDst   [[buffer(1)]],
    constant NufftSupportParams&  params [[buffer(2)]],
    uint gid [[thread_position_in_grid]])
{
    const int imageX = params.imageSize[0];
    const int imageY = params.imageSize[1];
    const int imageZ = params.imageSize[2];
    const int gNx    = params.gridSize[0];
    const int gNy    = params.gridSize[1];
    const int gNz    = params.gridSize[2];
    if ((int)gid >= imageX * imageY * imageZ) return;

    // Decompose: dst_idx = dZ*imageX*imageY + dX*imageY + dY
    const int dY_dst = (int)gid % imageY;
    const int dX_dst = ((int)gid / imageY) % imageX;
    const int dZ_dst = (int)gid / (imageX * imageY);

    const int offsetY = (gNy - imageY) / 2;
    const int offsetX = (gNx - imageX) / 2;
    const int offsetZ = (gNz - imageZ) / 2;

    const int dY_src = dY_dst + offsetY;
    const int dX_src = dX_dst + offsetX;
    const int dZ_src = dZ_dst + offsetZ;

    // CPU src index: dZ*gridSizeX*gridSizeY + dX*gridSizeY + dY
    const int src_idx = dZ_src * gNx * gNy + dX_src * gNy + dY_src;
    pDst[2 * (int)gid]     = pSrc[2 * src_idx];
    pDst[2 * (int)gid + 1] = pSrc[2 * src_idx + 1];
}
