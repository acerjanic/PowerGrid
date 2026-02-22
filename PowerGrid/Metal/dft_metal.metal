/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file dft_metal.metal
/// @brief Metal compute shaders for brute-force DFT (Gdft) and DFT with
///        gradient correction (GdftR2).
///
/// Kernels:
///   pg_dft_forward   — forward DFT: image → k-space (one thread per k-point)
///   pg_dft_adjoint   — adjoint DFT: k-space → image (one thread per image pixel)
///   pg_dft_forward_grads — forward DFT with sinc gradient correction (GdftR2)
///   pg_dft_adjoint_grads — adjoint DFT with sinc gradient correction (GdftR2)
///
/// Buffer layout:
///   [0] input real,  [1] input imag
///   [2] output real, [3] output imag
///   [4] kx, [5] ky, [6] kz
///   [7] ix, [8] iy, [9] iz
///   [10] FM, [11] t
///   [12] DFTParams (constant)
///   For WithGrads variants:
///   [13] Gx, [14] Gy, [15] Gz

#include <metal_stdlib>
using namespace metal;

constant float MRI_PI_F = 3.14159265358979323846f;

struct DFTParams {
    uint num_k;
    uint num_i;
    uint num_x;  // only used by WithGrads variants
    uint num_y;
    uint num_z;
};

// sincPG: sinc(x) = sin(pi*x)/(pi*x), returns 1.0 for |x| < 0.0001
inline float sincPG(float x) {
    if (abs(x) < 0.0001f) return 1.0f;
    float pix = MRI_PI_F * x;
    return sin(pix) / pix;
}

// ---- Forward DFT: image → k-space ----
// One thread per k-space point. Inner loop over all image points.
kernel void pg_dft_forward(
    device const float* idata_r [[buffer(0)]],
    device const float* idata_i [[buffer(1)]],
    device float* kdata_r       [[buffer(2)]],
    device float* kdata_i       [[buffer(3)]],
    device const float* kx      [[buffer(4)]],
    device const float* ky      [[buffer(5)]],
    device const float* kz      [[buffer(6)]],
    device const float* ix      [[buffer(7)]],
    device const float* iy      [[buffer(8)]],
    device const float* iz      [[buffer(9)]],
    device const float* FM      [[buffer(10)]],
    device const float* t       [[buffer(11)]],
    constant DFTParams& params  [[buffer(12)]],
    uint gid [[thread_position_in_grid]])
{
    if (gid >= params.num_k) return;

    float sumr = 0.0f;
    float sumi = 0.0f;

    float kxtpi = kx[gid] * 2.0f * MRI_PI_F;
    float kytpi = ky[gid] * 2.0f * MRI_PI_F;
    float kztpi = kz[gid] * 2.0f * MRI_PI_F;
    float myti  = t[gid];

    for (uint j = 0; j < params.num_i; j++) {
        float expr = kxtpi * ix[j] + kytpi * iy[j] + kztpi * iz[j] + FM[j] * myti;
        float c = cos(expr);
        float s = sin(expr);
        sumr += c * idata_r[j] + s * idata_i[j];
        sumi += -s * idata_r[j] + c * idata_i[j];
    }

    kdata_r[gid] = sumr;
    kdata_i[gid] = sumi;
}

// ---- Adjoint DFT: k-space → image ----
// One thread per image pixel. Inner loop over all k-space points.
kernel void pg_dft_adjoint(
    device const float* kdata_r [[buffer(0)]],
    device const float* kdata_i [[buffer(1)]],
    device float* idata_r       [[buffer(2)]],
    device float* idata_i       [[buffer(3)]],
    device const float* kx      [[buffer(4)]],
    device const float* ky      [[buffer(5)]],
    device const float* kz      [[buffer(6)]],
    device const float* ix      [[buffer(7)]],
    device const float* iy      [[buffer(8)]],
    device const float* iz      [[buffer(9)]],
    device const float* FM      [[buffer(10)]],
    device const float* t       [[buffer(11)]],
    constant DFTParams& params  [[buffer(12)]],
    uint gid [[thread_position_in_grid]])
{
    if (gid >= params.num_i) return;

    float sumr = 0.0f;
    float sumi = 0.0f;

    float tpi = 2.0f * MRI_PI_F;
    float ix_tpi = ix[gid] * tpi;
    float iy_tpi = iy[gid] * tpi;
    float iz_tpi = iz[gid] * tpi;
    float myfm   = FM[gid];

    for (uint i = 0; i < params.num_k; i++) {
        float expr = kx[i] * ix_tpi + ky[i] * iy_tpi + kz[i] * iz_tpi + myfm * t[i];
        float c = cos(expr);
        float s = sin(expr);
        sumr += c * kdata_r[i] - s * kdata_i[i];
        sumi += s * kdata_r[i] + c * kdata_i[i];
    }

    idata_r[gid] = sumr;
    idata_i[gid] = sumi;
}

// ---- Forward DFT with gradient correction (GdftR2) ----
kernel void pg_dft_forward_grads(
    device const float* idata_r [[buffer(0)]],
    device const float* idata_i [[buffer(1)]],
    device float* kdata_r       [[buffer(2)]],
    device float* kdata_i       [[buffer(3)]],
    device const float* kx      [[buffer(4)]],
    device const float* ky      [[buffer(5)]],
    device const float* kz      [[buffer(6)]],
    device const float* ix      [[buffer(7)]],
    device const float* iy      [[buffer(8)]],
    device const float* iz      [[buffer(9)]],
    device const float* FM      [[buffer(10)]],
    device const float* t       [[buffer(11)]],
    constant DFTParams& params  [[buffer(12)]],
    device const float* Gx      [[buffer(13)]],
    device const float* Gy      [[buffer(14)]],
    device const float* Gz      [[buffer(15)]],
    uint gid [[thread_position_in_grid]])
{
    if (gid >= params.num_k) return;

    float sumr = 0.0f;
    float sumi = 0.0f;

    float kxtpi = kx[gid] * 2.0f * MRI_PI_F;
    float kytpi = ky[gid] * 2.0f * MRI_PI_F;
    float kztpi = kz[gid] * 2.0f * MRI_PI_F;
    float myti  = t[gid];

    float kx_nx = kx[gid] / float(params.num_x);
    float ky_ny = ky[gid] / float(params.num_y);
    float kz_nz = kz[gid] / float(params.num_z);

    for (uint j = 0; j < params.num_i; j++) {
        float expr = kxtpi * ix[j] + kytpi * iy[j] + kztpi * iz[j] + FM[j] * myti;
        float c = cos(expr);
        float s = sin(expr);
        float bfunc = sincPG(kx_nx + Gx[j] * myti)
                    * sincPG(ky_ny + Gy[j] * myti)
                    * sincPG(kz_nz + Gz[j] * myti);
        sumr += bfunc * (c * idata_r[j] + s * idata_i[j]);
        sumi += bfunc * (-s * idata_r[j] + c * idata_i[j]);
    }

    kdata_r[gid] = sumr;
    kdata_i[gid] = sumi;
}

// ---- Adjoint DFT with gradient correction (GdftR2) ----
kernel void pg_dft_adjoint_grads(
    device const float* kdata_r [[buffer(0)]],
    device const float* kdata_i [[buffer(1)]],
    device float* idata_r       [[buffer(2)]],
    device float* idata_i       [[buffer(3)]],
    device const float* kx      [[buffer(4)]],
    device const float* ky      [[buffer(5)]],
    device const float* kz      [[buffer(6)]],
    device const float* ix      [[buffer(7)]],
    device const float* iy      [[buffer(8)]],
    device const float* iz      [[buffer(9)]],
    device const float* FM      [[buffer(10)]],
    device const float* t       [[buffer(11)]],
    constant DFTParams& params  [[buffer(12)]],
    device const float* Gx      [[buffer(13)]],
    device const float* Gy      [[buffer(14)]],
    device const float* Gz      [[buffer(15)]],
    uint gid [[thread_position_in_grid]])
{
    if (gid >= params.num_i) return;

    float sumr = 0.0f;
    float sumi = 0.0f;

    float tpi = 2.0f * MRI_PI_F;
    float ix_tpi = ix[gid] * tpi;
    float iy_tpi = iy[gid] * tpi;
    float iz_tpi = iz[gid] * tpi;
    float myfm   = FM[gid];
    float myGx   = Gx[gid];
    float myGy   = Gy[gid];
    float myGz   = Gz[gid];

    for (uint i = 0; i < params.num_k; i++) {
        float expr = kx[i] * ix_tpi + ky[i] * iy_tpi + kz[i] * iz_tpi + myfm * t[i];
        float c = cos(expr);
        float s = sin(expr);
        float bfunc = sincPG(kx[i] / float(params.num_x) + myGx * t[i])
                    * sincPG(ky[i] / float(params.num_y) + myGy * t[i])
                    * sincPG(kz[i] / float(params.num_z) + myGz * t[i]);
        sumr += bfunc * (c * kdata_r[i] - s * kdata_i[i]);
        sumi += bfunc * (s * kdata_r[i] + c * kdata_i[i]);
    }

    idata_r[gid] = sumr;
    idata_i[gid] = sumi;
}
