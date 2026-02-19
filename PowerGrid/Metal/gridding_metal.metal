/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file gridding_metal.metal
/// @brief Metal compute shaders for NUFFT gridding.
///
/// Kernels (dispatched from MetalGridding.mm):
///   pg_zero_init           — zero-initialise a float buffer (one thread per element)
///   pg_gridding_adjoint_2D — scatter k-space → 2D oversampled grid (atomic<float>)
///   pg_gridding_adjoint_3D — scatter k-space → 3D oversampled grid (atomic<float>)
///   pg_gridding_forward_2D — gather 2D grid → k-space (no atomics)
///   pg_gridding_forward_3D — gather 3D grid → k-space (no atomics)
///
/// Buffer index conventions:
///   Adjoint:  [0]=SampleF[], [1]=atomic_float grid[], [2]=LUT[], [3]=GridParams
///   Forward:  [0]=float grid[], [1]=float pSamples[], [2]=LUT[], [3]=GridParams,
///             [4]=kx[], [5]=ky[], ([6]=kz[] for 3D)
///
/// Requires MTLGPUFamilyApple6 (M1+) for `atomic<float>` support.

#include <metal_stdlib>
using namespace metal;

// ---------------------------------------------------------------------------
// Structs — must match griddingTypes.h for float instantiation and
// MetalGridding.mm GridParamsMSL.
// ---------------------------------------------------------------------------

/// @brief Gridding parameters (constant buffer 3).
struct GridParams {
    int numSamples;
    int imageSize[3];   ///< [Nx, Ny, Nz]
    int gridSize[3];    ///< [gNx, gNy, gNz] (oversampled dimensions)
    float gridOS;
    float kernelWidth;
    int sizeLUT;
    int pad;            ///< struct alignment padding
};

/// @brief One k-space sample — matches ReconstructionSample<float> layout.
struct SampleF {
    float real;
    float imag;
    float kX;
    float kY;
    float kZ;
    float sdc;
    float t;
    float dummy;
};

// ---------------------------------------------------------------------------
// pg_zero_init
// ---------------------------------------------------------------------------

/// @brief Set every float element of @p buf to 0.0.  One thread per element.
kernel void pg_zero_init(
    device float* buf  [[buffer(0)]],
    uint gid           [[thread_position_in_grid]])
{
    buf[gid] = 0.0f;
}

// ---------------------------------------------------------------------------
// pg_gridding_adjoint_2D
// ---------------------------------------------------------------------------

/// @brief 2D adjoint NUFFT gridding: scatter one k-space sample onto the grid.
///
/// Grid index layout (matches gridding.cpp): idx = ny + nx * gNy
kernel void pg_gridding_adjoint_2D(
    device const SampleF*  samples [[buffer(0)]],
    device atomic_float*   pGData  [[buffer(1)]],  ///< Re/Im interleaved, 2*gNx*gNy
    constant const float*  LUT     [[buffer(2)]],
    constant GridParams&   params  [[buffer(3)]],
    uint gid [[thread_position_in_grid]])
{
    if ((int)gid >= params.numSamples) return;

    const SampleF pt        = samples[gid];
    const float   gridOS    = params.gridOS;
    const float   kw        = params.kernelWidth;
    const float   kw2       = kw * kw;
    const int     sizeLUT   = params.sizeLUT;
    const int     Nx        = params.imageSize[0];
    const int     Ny        = params.imageSize[1];
    const int     gNx       = params.gridSize[0];
    const int     gNy       = params.gridSize[1];

    const float shKx = gridOS * (pt.kX + (float)Nx * 0.5f);
    const float shKy = gridOS * (pt.kY + (float)Ny * 0.5f);

    const int NxL = max(0,     (int)ceil (shKx - kw * gridOS * 0.5f));
    const int NxH = min(gNx-1, (int)floor(shKx + kw * gridOS * 0.5f));
    const int NyL = max(0,     (int)ceil (shKy - kw * gridOS * 0.5f));
    const int NyH = min(gNy-1, (int)floor(shKy + kw * gridOS * 0.5f));

    for (int nx = NxL; nx <= NxH; ++nx) {
        const float dx  = fabs(shKx - (float)nx) / gridOS;
        const int   k0x = (int)(dx * dx * 4.0f / kw2 * (float)sizeLUT);
        const float kbX = (k0x < sizeLUT) ? LUT[k0x] : 0.0f;
        if (kbX == 0.0f) continue;

        for (int ny = NyL; ny <= NyH; ++ny) {
            const float dy  = fabs(shKy - (float)ny) / gridOS;
            const int   k0y = (int)(dy * dy * 4.0f / kw2 * (float)sizeLUT);
            const float w   = kbX * ((k0y < sizeLUT) ? LUT[k0y] : 0.0f);
            if (w == 0.0f) continue;

            const int idx = 2 * (ny + nx * gNy);
            atomic_fetch_add_explicit(&pGData[idx],     w * pt.real, memory_order_relaxed);
            atomic_fetch_add_explicit(&pGData[idx + 1], w * pt.imag, memory_order_relaxed);
        }
    }
}

// ---------------------------------------------------------------------------
// pg_gridding_adjoint_3D
// ---------------------------------------------------------------------------

/// @brief 3D adjoint NUFFT gridding: scatter one k-space sample onto the 3D grid.
///
/// Grid index layout (matches gridding.cpp): idx = ny + nx*gNy + nz*gNx*gNy
kernel void pg_gridding_adjoint_3D(
    device const SampleF*  samples [[buffer(0)]],
    device atomic_float*   pGData  [[buffer(1)]],  ///< Re/Im interleaved, 2*gNx*gNy*gNz
    constant const float*  LUT     [[buffer(2)]],
    constant GridParams&   params  [[buffer(3)]],
    uint gid [[thread_position_in_grid]])
{
    if ((int)gid >= params.numSamples) return;

    const SampleF pt        = samples[gid];
    const float   gridOS    = params.gridOS;
    const float   kw        = params.kernelWidth;
    const float   kw2       = kw * kw;
    const int     sizeLUT   = params.sizeLUT;
    const int     Nx        = params.imageSize[0];
    const int     Ny        = params.imageSize[1];
    const int     Nz        = params.imageSize[2];
    const int     gNx       = params.gridSize[0];
    const int     gNy       = params.gridSize[1];
    const int     gNz       = params.gridSize[2];

    const float shKx = gridOS * (pt.kX + (float)Nx * 0.5f);
    const float shKy = gridOS * (pt.kY + (float)Ny * 0.5f);
    const float shKz = gridOS * (pt.kZ + (float)Nz * 0.5f);

    const int NxL = max(0,     (int)ceil (shKx - kw * gridOS * 0.5f));
    const int NxH = min(gNx-1, (int)floor(shKx + kw * gridOS * 0.5f));
    const int NyL = max(0,     (int)ceil (shKy - kw * gridOS * 0.5f));
    const int NyH = min(gNy-1, (int)floor(shKy + kw * gridOS * 0.5f));
    const int NzL = max(0,     (int)ceil (shKz - kw * gridOS * 0.5f));
    const int NzH = min(gNz-1, (int)floor(shKz + kw * gridOS * 0.5f));

    for (int nz = NzL; nz <= NzH; ++nz) {
        const float dz  = fabs(shKz - (float)nz) / gridOS;
        const int   k0z = (int)(dz * dz * 4.0f / kw2 * (float)sizeLUT);
        const float kbZ = (k0z < sizeLUT) ? LUT[k0z] : 0.0f;
        if (kbZ == 0.0f) continue;

        for (int nx = NxL; nx <= NxH; ++nx) {
            const float dx  = fabs(shKx - (float)nx) / gridOS;
            const int   k0x = (int)(dx * dx * 4.0f / kw2 * (float)sizeLUT);
            const float kbX = (k0x < sizeLUT) ? LUT[k0x] : 0.0f;
            if (kbX == 0.0f) continue;

            for (int ny = NyL; ny <= NyH; ++ny) {
                const float dy  = fabs(shKy - (float)ny) / gridOS;
                const int   k0y = (int)(dy * dy * 4.0f / kw2 * (float)sizeLUT);
                const float w   = kbX * kbZ * ((k0y < sizeLUT) ? LUT[k0y] : 0.0f);
                if (w == 0.0f) continue;

                const int idx = 2 * (ny + nx * gNy + nz * gNx * gNy);
                atomic_fetch_add_explicit(&pGData[idx],     w * pt.real, memory_order_relaxed);
                atomic_fetch_add_explicit(&pGData[idx + 1], w * pt.imag, memory_order_relaxed);
            }
        }
    }
}

// ---------------------------------------------------------------------------
// pg_gridding_forward_2D
// ---------------------------------------------------------------------------

/// @brief 2D forward NUFFT gridding: gather from grid into one k-space sample.
///
/// No atomics — each thread writes exclusively to its own output index.
/// k-space coordinates supplied via separate float buffers (buffers 4, 5).
kernel void pg_gridding_forward_2D(
    device const float*    pGridData [[buffer(0)]],   ///< 2*gNx*gNy (input grid)
    device float*          pSamples  [[buffer(1)]],   ///< 2*numSamples (output, pre-zeroed)
    constant const float*  LUT       [[buffer(2)]],
    constant GridParams&   params    [[buffer(3)]],
    device const float*    kx        [[buffer(4)]],
    device const float*    ky        [[buffer(5)]],
    uint gid [[thread_position_in_grid]])
{
    if ((int)gid >= params.numSamples) return;

    const float gridOS  = params.gridOS;
    const float kw      = params.kernelWidth;
    const float kw2     = kw * kw;
    const int   sizeLUT = params.sizeLUT;
    const int   Nx      = params.imageSize[0];
    const int   Ny      = params.imageSize[1];
    const int   gNx     = params.gridSize[0];
    const int   gNy     = params.gridSize[1];

    const float shKx = gridOS * (kx[gid] + (float)Nx * 0.5f);
    const float shKy = gridOS * (ky[gid] + (float)Ny * 0.5f);

    const int NxL = max(0,     (int)ceil (shKx - kw * gridOS * 0.5f));
    const int NxH = min(gNx-1, (int)floor(shKx + kw * gridOS * 0.5f));
    const int NyL = max(0,     (int)ceil (shKy - kw * gridOS * 0.5f));
    const int NyH = min(gNy-1, (int)floor(shKy + kw * gridOS * 0.5f));

    float accReal = 0.0f;
    float accImag = 0.0f;

    for (int nx = NxL; nx <= NxH; ++nx) {
        const float dx  = fabs(shKx - (float)nx) / gridOS;
        const int   k0x = (int)(dx * dx * 4.0f / kw2 * (float)sizeLUT);
        const float kbX = (k0x < sizeLUT) ? LUT[k0x] : 0.0f;
        if (kbX == 0.0f) continue;

        for (int ny = NyL; ny <= NyH; ++ny) {
            const float dy  = fabs(shKy - (float)ny) / gridOS;
            const int   k0y = (int)(dy * dy * 4.0f / kw2 * (float)sizeLUT);
            const float w   = kbX * ((k0y < sizeLUT) ? LUT[k0y] : 0.0f);

            const int idx = 2 * (ny + nx * gNy);
            accReal += w * pGridData[idx];
            accImag += w * pGridData[idx + 1];
        }
    }

    pSamples[2 * (int)gid]     = accReal;
    pSamples[2 * (int)gid + 1] = accImag;
}

// ---------------------------------------------------------------------------
// pg_gridding_forward_3D
// ---------------------------------------------------------------------------

/// @brief 3D forward NUFFT gridding: gather from 3D grid into one k-space sample.
kernel void pg_gridding_forward_3D(
    device const float*    pGridData [[buffer(0)]],   ///< 2*gNx*gNy*gNz (input grid)
    device float*          pSamples  [[buffer(1)]],   ///< 2*numSamples (output, pre-zeroed)
    constant const float*  LUT       [[buffer(2)]],
    constant GridParams&   params    [[buffer(3)]],
    device const float*    kx        [[buffer(4)]],
    device const float*    ky        [[buffer(5)]],
    device const float*    kz        [[buffer(6)]],
    uint gid [[thread_position_in_grid]])
{
    if ((int)gid >= params.numSamples) return;

    const float gridOS  = params.gridOS;
    const float kw      = params.kernelWidth;
    const float kw2     = kw * kw;
    const int   sizeLUT = params.sizeLUT;
    const int   Nx      = params.imageSize[0];
    const int   Ny      = params.imageSize[1];
    const int   Nz      = params.imageSize[2];
    const int   gNx     = params.gridSize[0];
    const int   gNy     = params.gridSize[1];
    const int   gNz     = params.gridSize[2];

    const float shKx = gridOS * (kx[gid] + (float)Nx * 0.5f);
    const float shKy = gridOS * (ky[gid] + (float)Ny * 0.5f);
    const float shKz = gridOS * (kz[gid] + (float)Nz * 0.5f);

    const int NxL = max(0,     (int)ceil (shKx - kw * gridOS * 0.5f));
    const int NxH = min(gNx-1, (int)floor(shKx + kw * gridOS * 0.5f));
    const int NyL = max(0,     (int)ceil (shKy - kw * gridOS * 0.5f));
    const int NyH = min(gNy-1, (int)floor(shKy + kw * gridOS * 0.5f));
    const int NzL = max(0,     (int)ceil (shKz - kw * gridOS * 0.5f));
    const int NzH = min(gNz-1, (int)floor(shKz + kw * gridOS * 0.5f));

    float accReal = 0.0f;
    float accImag = 0.0f;

    for (int nz = NzL; nz <= NzH; ++nz) {
        const float dz  = fabs(shKz - (float)nz) / gridOS;
        const int   k0z = (int)(dz * dz * 4.0f / kw2 * (float)sizeLUT);
        const float kbZ = (k0z < sizeLUT) ? LUT[k0z] : 0.0f;
        if (kbZ == 0.0f) continue;

        for (int nx = NxL; nx <= NxH; ++nx) {
            const float dx  = fabs(shKx - (float)nx) / gridOS;
            const int   k0x = (int)(dx * dx * 4.0f / kw2 * (float)sizeLUT);
            const float kbX = (k0x < sizeLUT) ? LUT[k0x] : 0.0f;
            if (kbX == 0.0f) continue;

            for (int ny = NyL; ny <= NyH; ++ny) {
                const float dy  = fabs(shKy - (float)ny) / gridOS;
                const int   k0y = (int)(dy * dy * 4.0f / kw2 * (float)sizeLUT);
                const float w   = kbX * kbZ * ((k0y < sizeLUT) ? LUT[k0y] : 0.0f);

                const int idx = 2 * (ny + nx * gNy + nz * gNx * gNy);
                accReal += w * pGridData[idx];
                accImag += w * pGridData[idx + 1];
            }
        }
    }

    pSamples[2 * (int)gid]     = accReal;
    pSamples[2 * (int)gid + 1] = accImag;
}
