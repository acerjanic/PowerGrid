/*
(C) Copyright 2015-2024 The Board of Trustees of the University of Illinois.
All rights reserved.

See LICENSE.txt for the University of Illinois/NCSA Open Source license.

Developed by:
                     MRFIL Research Groups
                University of Illinois, Urbana-Champaign
*/

/// @file MetalVectorOpsBenchmarks.cpp
/// @brief Standalone benchmark: Accelerate vs Metal-memcpy vs Metal-zerocopy.
///
/// Tests each vector operation at multiple sizes to find the crossover point
/// where Metal GPU dispatch outperforms CPU/Accelerate.
///
/// Build with METAL_COMPUTE=ON, then run:
///   ./metal_vecops_bench
///
/// Override Metal dispatch threshold:
///   PG_METAL_MIN_SIZE=256 ./metal_vecops_bench

#ifdef METAL_COMPUTE

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <vector>

#include "Core/PGIncludes.h"
#include "Core/pgCol.hpp"
#include "Core/pgComplex.hpp"
#include "Metal/MetalVectorOps.h"
#include "AccelerateDispatch.hpp"

#include <Accelerate/Accelerate.h>

using Clock = std::chrono::high_resolution_clock;
using Ms    = std::chrono::duration<double, std::milli>;

template<typename T>
static void doNotOptimize(const T& val) {
    asm volatile("" : : "r,m"(&val) : "memory");
}

template <typename Fn>
static double median_ms(Fn fn, int nwarmup = 3, int nruns = 15) {
    for (int i = 0; i < nwarmup; i++) fn();
    std::vector<double> t(nruns);
    for (int i = 0; i < nruns; i++) {
        auto t0 = Clock::now();
        fn();
        t[i] = Ms(Clock::now() - t0).count();
    }
    std::sort(t.begin(), t.end());
    return t[nruns / 2];
}

/// Allocate page-aligned buffer (16384 = Apple Silicon page size).
static float* pageAlloc(size_t count) {
    size_t bytes = count * sizeof(float);
    size_t pageSize = 16384;
    size_t aligned = (bytes + pageSize - 1) & ~(pageSize - 1);
    return (float*)std::aligned_alloc(pageSize, aligned);
}

/// Throughput in GB/s given bytes transferred and time in ms.
static double gbps(size_t bytes, double ms) {
    return (ms > 0.0) ? (double)bytes / (ms * 1e6) : 0.0;
}

// ---------------------------------------------------------------------------
// Benchmark runner: prints one line per size comparing 3 backends.
// ---------------------------------------------------------------------------
struct BenchResult {
    double accel_ms;
    double metal_ms;
    double metal_zc_ms;
};

static void printRow(size_t N, const char* name, size_t bytesPerElem,
                     const BenchResult& r) {
    size_t totalBytes = N * bytesPerElem * 3; // rough: 2 inputs + 1 output
    printf("  N=%-8zu  %-20s  Accel: %7.3fms (%5.1f GB/s)  "
           "Metal: %7.3fms (%5.1f GB/s)  "
           "Metal-zc: %7.3fms (%5.1f GB/s)  "
           "zc vs Accel: %.1fx\n",
           N, name,
           r.accel_ms, gbps(totalBytes, r.accel_ms),
           r.metal_ms, gbps(totalBytes, r.metal_ms),
           r.metal_zc_ms, gbps(totalBytes, r.metal_zc_ms),
           (r.metal_zc_ms > 0.0) ? r.accel_ms / r.metal_zc_ms : 0.0);
}

// ---------------------------------------------------------------------------
int main()
// ---------------------------------------------------------------------------
{
    auto* mctx = metal_vector_get_context();
    if (!mctx) {
        printf("ERROR: Metal GPU not available.\n");
        return 1;
    }
    printf("=== Metal Vector Ops Benchmark: Accelerate vs Metal (memcpy) vs Metal (zero-copy) ===\n");
    printf("    PG_METAL_MIN_SIZE = %zu (set env var to override)\n\n",
           pg_metal::metalMinSize());

    const std::vector<size_t> sizes = {1024, 4096, 16384, 65536, 262144, 864000, 1048576, 4194304};

    // =====================================================================
    // 1. Complex element-wise multiply: C = A .* B
    // =====================================================================
    printf("--- Complex element-wise multiply (C = A .* B) ---\n");
    for (size_t N : sizes) {
        float* A = pageAlloc(2 * N);
        float* B = pageAlloc(2 * N);
        float* C = pageAlloc(2 * N);
        for (size_t i = 0; i < 2 * N; i++) {
            A[i] = (float)(i % 1000) * 0.001f - 0.5f;
            B[i] = (float)((i + 333) % 1000) * 0.001f - 0.5f;
        }

        BenchResult r;
        // Accelerate: use our optimized try_accel_mul
        r.accel_ms = median_ms([&]{
            pg_accel::try_accel_mul(
                reinterpret_cast<const pgComplex<float>*>(A),
                reinterpret_cast<const pgComplex<float>*>(B),
                reinterpret_cast<pgComplex<float>*>(C), N);
            doNotOptimize(C[0]);
        });
        // Metal memcpy
        r.metal_ms = median_ms([&]{
            metal_cvec_mul(mctx, A, B, C, N);
            doNotOptimize(C[0]);
        });
        // Metal zero-copy
        r.metal_zc_ms = median_ms([&]{
            metal_cvec_mul_zc(mctx, A, B, C, N);
            doNotOptimize(C[0]);
        });

        printRow(N, "cvec_mul", 2 * sizeof(float), r);
        std::free(A); std::free(B); std::free(C);
    }
    printf("\n");

    // =====================================================================
    // 2. Complex element-wise divide: C = A ./ B
    // =====================================================================
    printf("--- Complex element-wise divide (C = A ./ B) ---\n");
    for (size_t N : sizes) {
        float* A = pageAlloc(2 * N);
        float* B = pageAlloc(2 * N);
        float* C = pageAlloc(2 * N);
        for (size_t i = 0; i < 2 * N; i++) {
            A[i] = (float)(i % 1000) * 0.001f - 0.5f;
            B[i] = (float)((i + 333) % 1000) * 0.001f + 0.1f; // avoid div-by-zero
        }

        BenchResult r;
        r.accel_ms = median_ms([&]{
            pg_accel::try_accel_div(
                reinterpret_cast<const pgComplex<float>*>(A),
                reinterpret_cast<const pgComplex<float>*>(B),
                reinterpret_cast<pgComplex<float>*>(C), N);
            doNotOptimize(C[0]);
        });
        r.metal_ms = median_ms([&]{
            metal_cvec_div(mctx, A, B, C, N);
            doNotOptimize(C[0]);
        });
        r.metal_zc_ms = median_ms([&]{
            metal_cvec_div_zc(mctx, A, B, C, N);
            doNotOptimize(C[0]);
        });

        printRow(N, "cvec_div", 2 * sizeof(float), r);
        std::free(A); std::free(B); std::free(C);
    }
    printf("\n");

    // =====================================================================
    // 3. Complex scalar multiply: C = alpha * A
    // =====================================================================
    printf("--- Complex scalar multiply (C = alpha * A) ---\n");
    for (size_t N : sizes) {
        float* A = pageAlloc(2 * N);
        float* C = pageAlloc(2 * N);
        for (size_t i = 0; i < 2 * N; i++)
            A[i] = (float)(i % 1000) * 0.001f - 0.5f;
        pgComplex<float> alpha(0.7f, -0.3f);

        BenchResult r;
        r.accel_ms = median_ms([&]{
            pg_accel::try_accel_mul_scalar(
                reinterpret_cast<const pgComplex<float>*>(A),
                alpha,
                reinterpret_cast<pgComplex<float>*>(C), N);
            doNotOptimize(C[0]);
        });
        r.metal_ms = median_ms([&]{
            metal_cvec_mul_scalar(mctx, A, alpha.real(), alpha.imag(), C, N);
            doNotOptimize(C[0]);
        });
        r.metal_zc_ms = median_ms([&]{
            metal_cvec_mul_scalar_zc(mctx, A, alpha.real(), alpha.imag(), C, N);
            doNotOptimize(C[0]);
        });

        printRow(N, "cvec_mul_scalar", 2 * sizeof(float), r);
        std::free(A); std::free(C);
    }
    printf("\n");

    // =====================================================================
    // 4. Complex AXPY: C = A + alpha * B
    // =====================================================================
    printf("--- Complex AXPY (C = A + alpha * B) ---\n");
    for (size_t N : sizes) {
        float* A = pageAlloc(2 * N);
        float* B = pageAlloc(2 * N);
        float* C = pageAlloc(2 * N);
        for (size_t i = 0; i < 2 * N; i++) {
            A[i] = (float)(i % 1000) * 0.001f - 0.5f;
            B[i] = (float)((i + 333) % 1000) * 0.001f - 0.5f;
        }
        pgComplex<float> alpha(0.7f, -0.3f);

        BenchResult r;
        r.accel_ms = median_ms([&]{
            pg_accel::try_accel_axpy(
                reinterpret_cast<const pgComplex<float>*>(A),
                reinterpret_cast<const pgComplex<float>*>(B),
                reinterpret_cast<pgComplex<float>*>(C),
                alpha, N);
            doNotOptimize(C[0]);
        });
        r.metal_ms = median_ms([&]{
            metal_cvec_axpy(mctx, A, B, C, alpha.real(), alpha.imag(), N);
            doNotOptimize(C[0]);
        });
        r.metal_zc_ms = median_ms([&]{
            metal_cvec_axpy_zc(mctx, A, B, C, alpha.real(), alpha.imag(), N);
            doNotOptimize(C[0]);
        });

        printRow(N, "cvec_axpy", 2 * sizeof(float), r);
        std::free(A); std::free(B); std::free(C);
    }
    printf("\n");

    // =====================================================================
    // 5. Complex add: C = A + B
    // =====================================================================
    printf("--- Complex element-wise add (C = A + B) ---\n");
    for (size_t N : sizes) {
        float* A = pageAlloc(2 * N);
        float* B = pageAlloc(2 * N);
        float* C = pageAlloc(2 * N);
        for (size_t i = 0; i < 2 * N; i++) {
            A[i] = (float)(i % 1000) * 0.001f;
            B[i] = (float)((i + 500) % 1000) * 0.001f;
        }

        BenchResult r;
        r.accel_ms = median_ms([&]{
            pg_accel::try_accel_add(
                reinterpret_cast<const pgComplex<float>*>(A),
                reinterpret_cast<const pgComplex<float>*>(B),
                reinterpret_cast<pgComplex<float>*>(C), N);
            doNotOptimize(C[0]);
        });
        r.metal_ms = median_ms([&]{
            metal_cvec_add(mctx, A, B, C, N);
            doNotOptimize(C[0]);
        });
        r.metal_zc_ms = median_ms([&]{
            metal_cvec_add_zc(mctx, A, B, C, N);
            doNotOptimize(C[0]);
        });

        printRow(N, "cvec_add", 2 * sizeof(float), r);
        std::free(A); std::free(B); std::free(C);
    }
    printf("\n");

    // =====================================================================
    // 6. Complex sub: C = A - B
    // =====================================================================
    printf("--- Complex element-wise sub (C = A - B) ---\n");
    for (size_t N : sizes) {
        float* A = pageAlloc(2 * N);
        float* B = pageAlloc(2 * N);
        float* C = pageAlloc(2 * N);
        for (size_t i = 0; i < 2 * N; i++) {
            A[i] = (float)(i % 1000) * 0.001f;
            B[i] = (float)((i + 500) % 1000) * 0.001f;
        }

        BenchResult r;
        r.accel_ms = median_ms([&]{
            pg_accel::try_accel_sub(
                reinterpret_cast<const pgComplex<float>*>(A),
                reinterpret_cast<const pgComplex<float>*>(B),
                reinterpret_cast<pgComplex<float>*>(C), N);
            doNotOptimize(C[0]);
        });
        r.metal_ms = median_ms([&]{
            metal_cvec_sub(mctx, A, B, C, N);
            doNotOptimize(C[0]);
        });
        r.metal_zc_ms = median_ms([&]{
            metal_cvec_sub_zc(mctx, A, B, C, N);
            doNotOptimize(C[0]);
        });

        printRow(N, "cvec_sub", 2 * sizeof(float), r);
        std::free(A); std::free(B); std::free(C);
    }
    printf("\n");

    // =====================================================================
    // 7. Real-weight * complex: C = W .* X
    // =====================================================================
    printf("--- Real-weight * complex (C = W .* X) ---\n");
    for (size_t N : sizes) {
        float* W = pageAlloc(N);
        float* X = pageAlloc(2 * N);
        float* C = pageAlloc(2 * N);
        for (size_t i = 0; i < N; i++)
            W[i] = (float)(i % 1000) * 0.001f;
        for (size_t i = 0; i < 2 * N; i++)
            X[i] = (float)((i + 333) % 1000) * 0.001f - 0.5f;

        BenchResult r;
        r.accel_ms = median_ms([&]{
            pg_accel::try_accel_rvec_cmul(
                W,
                reinterpret_cast<const pgComplex<float>*>(X),
                reinterpret_cast<pgComplex<float>*>(C), N);
            doNotOptimize(C[0]);
        });
        r.metal_ms = median_ms([&]{
            metal_rvec_cmul(mctx, W, X, C, N);
            doNotOptimize(C[0]);
        });
        r.metal_zc_ms = median_ms([&]{
            metal_rvec_cmul_zc(mctx, W, X, C, N);
            doNotOptimize(C[0]);
        });

        // bytes: N*4 (W) + N*8 (X) + N*8 (C) = N*20
        size_t totalBytes = N * 20;
        printf("  N=%-8zu  %-20s  Accel: %7.3fms (%5.1f GB/s)  "
               "Metal: %7.3fms (%5.1f GB/s)  "
               "Metal-zc: %7.3fms (%5.1f GB/s)  "
               "zc vs Accel: %.1fx\n",
               N, "rvec_cmul",
               r.accel_ms, gbps(totalBytes, r.accel_ms),
               r.metal_ms, gbps(totalBytes, r.metal_ms),
               r.metal_zc_ms, gbps(totalBytes, r.metal_zc_ms),
               (r.metal_zc_ms > 0.0) ? r.accel_ms / r.metal_zc_ms : 0.0);
        std::free(W); std::free(X); std::free(C);
    }
    printf("\n");

    // =====================================================================
    // 8. Complex dot product: sum(conj(A) .* B)
    // =====================================================================
    printf("--- Complex dot product (cdot) ---\n");
    for (size_t N : sizes) {
        float* A = pageAlloc(2 * N);
        float* B = pageAlloc(2 * N);
        for (size_t i = 0; i < 2 * N; i++) {
            A[i] = (float)(i % 1000) * 0.001f - 0.5f;
            B[i] = (float)((i + 333) % 1000) * 0.001f - 0.5f;
        }

        BenchResult r;
        float reA, imA, reM, imM, reZ, imZ;
        r.accel_ms = median_ms([&]{
            pgComplex<float> result;
            pg_accel::try_accel_cdot(
                reinterpret_cast<const pgComplex<float>*>(A),
                reinterpret_cast<const pgComplex<float>*>(B),
                &result, N);
            reA = result.real(); imA = result.imag();
            doNotOptimize(result);
        });
        r.metal_ms = median_ms([&]{
            metal_cvec_cdot(mctx, A, B, &reM, &imM, N);
            doNotOptimize(reM);
        });
        r.metal_zc_ms = median_ms([&]{
            metal_cvec_cdot_zc(mctx, A, B, &reZ, &imZ, N);
            doNotOptimize(reZ);
        });

        size_t totalBytes = N * 2 * sizeof(float) * 2; // 2 input vectors
        printf("  N=%-8zu  %-20s  Accel: %7.3fms (%5.1f GB/s)  "
               "Metal: %7.3fms (%5.1f GB/s)  "
               "Metal-zc: %7.3fms (%5.1f GB/s)  "
               "zc vs Accel: %.1fx\n",
               N, "cdot",
               r.accel_ms, gbps(totalBytes, r.accel_ms),
               r.metal_ms, gbps(totalBytes, r.metal_ms),
               r.metal_zc_ms, gbps(totalBytes, r.metal_zc_ms),
               (r.metal_zc_ms > 0.0) ? r.accel_ms / r.metal_zc_ms : 0.0);
        std::free(A); std::free(B);
    }
    printf("\n");

    // =====================================================================
    // 9. Complex norm2sq: sum(|A|^2)
    // =====================================================================
    printf("--- Complex norm2sq (sum |A|^2) ---\n");
    for (size_t N : sizes) {
        float* A = pageAlloc(2 * N);
        for (size_t i = 0; i < 2 * N; i++)
            A[i] = (float)(i % 1000) * 0.001f - 0.5f;

        BenchResult r;
        r.accel_ms = median_ms([&]{
            float n2;
            pg_accel::try_accel_norm2sq(
                reinterpret_cast<const pgComplex<float>*>(A), &n2, N);
            doNotOptimize(n2);
        });
        r.metal_ms = median_ms([&]{
            float n2 = metal_cvec_norm2sq(mctx, A, N);
            doNotOptimize(n2);
        });
        r.metal_zc_ms = median_ms([&]{
            float n2 = metal_cvec_norm2sq_zc(mctx, A, N);
            doNotOptimize(n2);
        });

        size_t totalBytes = N * 2 * sizeof(float); // 1 input vector
        printf("  N=%-8zu  %-20s  Accel: %7.3fms (%5.1f GB/s)  "
               "Metal: %7.3fms (%5.1f GB/s)  "
               "Metal-zc: %7.3fms (%5.1f GB/s)  "
               "zc vs Accel: %.1fx\n",
               N, "norm2sq",
               r.accel_ms, gbps(totalBytes, r.accel_ms),
               r.metal_ms, gbps(totalBytes, r.metal_ms),
               r.metal_zc_ms, gbps(totalBytes, r.metal_zc_ms),
               (r.metal_zc_ms > 0.0) ? r.accel_ms / r.metal_zc_ms : 0.0);
        std::free(A);
    }
    printf("\n");

    // =====================================================================
    // 10. Real add: C = A + B
    // =====================================================================
    printf("--- Real element-wise add (C = A + B) ---\n");
    for (size_t N : sizes) {
        float* A = pageAlloc(N);
        float* B = pageAlloc(N);
        float* C = pageAlloc(N);
        for (size_t i = 0; i < N; i++) {
            A[i] = (float)(i % 1000) * 0.001f;
            B[i] = (float)((i + 500) % 1000) * 0.001f;
        }

        BenchResult r;
        r.accel_ms = median_ms([&]{
            vDSP_vadd(A, 1, B, 1, C, 1, (vDSP_Length)N);
            doNotOptimize(C[0]);
        });
        r.metal_ms = median_ms([&]{
            metal_vec_add(mctx, A, B, C, N);
            doNotOptimize(C[0]);
        });
        r.metal_zc_ms = median_ms([&]{
            metal_vec_add_zc(mctx, A, B, C, N);
            doNotOptimize(C[0]);
        });

        printRow(N, "real_add", sizeof(float), r);
        std::free(A); std::free(B); std::free(C);
    }
    printf("\n");

    // =====================================================================
    // 11. Real mul: C = A .* B
    // =====================================================================
    printf("--- Real element-wise mul (C = A .* B) ---\n");
    for (size_t N : sizes) {
        float* A = pageAlloc(N);
        float* B = pageAlloc(N);
        float* C = pageAlloc(N);
        for (size_t i = 0; i < N; i++) {
            A[i] = (float)(i % 1000) * 0.001f;
            B[i] = (float)((i + 500) % 1000) * 0.001f;
        }

        BenchResult r;
        r.accel_ms = median_ms([&]{
            vDSP_vmul(A, 1, B, 1, C, 1, (vDSP_Length)N);
            doNotOptimize(C[0]);
        });
        r.metal_ms = median_ms([&]{
            metal_vec_mul(mctx, A, B, C, N);
            doNotOptimize(C[0]);
        });
        r.metal_zc_ms = median_ms([&]{
            metal_vec_mul_zc(mctx, A, B, C, N);
            doNotOptimize(C[0]);
        });

        printRow(N, "real_mul", sizeof(float), r);
        std::free(A); std::free(B); std::free(C);
    }
    printf("\n");

    // =====================================================================
    // 12. Real sum: sum(A)
    // =====================================================================
    printf("--- Real sum (sum(A)) ---\n");
    for (size_t N : sizes) {
        float* A = pageAlloc(N);
        for (size_t i = 0; i < N; i++)
            A[i] = (float)(i % 1000) * 0.001f - 0.5f;

        BenchResult r;
        r.accel_ms = median_ms([&]{
            float s;
            vDSP_sve(A, 1, &s, (vDSP_Length)N);
            doNotOptimize(s);
        });
        r.metal_ms = median_ms([&]{
            float s = metal_vec_sum(mctx, A, N);
            doNotOptimize(s);
        });
        r.metal_zc_ms = median_ms([&]{
            float s = metal_vec_sum_zc(mctx, A, N);
            doNotOptimize(s);
        });

        size_t totalBytes = N * sizeof(float);
        printf("  N=%-8zu  %-20s  Accel: %7.3fms (%5.1f GB/s)  "
               "Metal: %7.3fms (%5.1f GB/s)  "
               "Metal-zc: %7.3fms (%5.1f GB/s)  "
               "zc vs Accel: %.1fx\n",
               N, "real_sum",
               r.accel_ms, gbps(totalBytes, r.accel_ms),
               r.metal_ms, gbps(totalBytes, r.metal_ms),
               r.metal_zc_ms, gbps(totalBytes, r.metal_zc_ms),
               (r.metal_zc_ms > 0.0) ? r.accel_ms / r.metal_zc_ms : 0.0);
        std::free(A);
    }
    printf("\n");

    printf("Done.\n");
    return 0;
}

#else // METAL_COMPUTE not defined

#include <cstdio>
int main() {
    puts("METAL_COMPUTE is not enabled. Rebuild with -DMETAL_COMPUTE=ON.");
    return 1;
}

#endif // METAL_COMPUTE
