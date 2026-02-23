
/*****************************************************************************

    File Name   [Tracer.hpp]

    Synopsis    [Tracer class for inserting nvtx range markers using RAII idiom.
                Look at NVIDIA dev blog for more details.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/10/2018]

*****************************************************************************/

/// @file Tracer.hpp
/// @brief RAII profiling markers for NVIDIA Nsight (NVTX) and Apple Instruments (os_signpost).

#ifndef PowerGrid_Tracer_hpp
#define PowerGrid_Tracer_hpp

// Token-paste helper (two-level indirection for __LINE__ expansion)
#define PG_CAT_(a, b) a ## b
#define PG_CAT(a, b) PG_CAT_(a, b)

#if defined(USE_NVTX)
// ---- NVIDIA NVTX backend ----
/// @brief RAII wrapper that pushes/pops an NVTX named range for GPU profiling.
class Tracer {
public:
    Tracer(const char* name = "") { nvtxRangePushA(name); }
    ~Tracer() { nvtxRangePop(); }
};
#define RANGE(...) Tracer PG_CAT(_nvtx_, __LINE__){__VA_ARGS__};

#elif defined(__APPLE__) && defined(USE_SIGNPOST)
// ---- Apple os_signpost backend (Instruments.app) ----
#include <os/log.h>
#include <os/signpost.h>

/// @brief RAII wrapper that emits os_signpost intervals for Apple Instruments.
///
/// Construct a `SignpostTracer` at the start of a code region to begin a
/// signpost interval; the destructor ends it. Intervals appear in the
/// "Points of Interest" track in Instruments.app.
class SignpostTracer {
public:
    SignpostTracer(const char* name = "(scope)") : name_(name) {
        log_ = os_log_create("com.powergrid.profiling", "Recon");
        spid_ = os_signpost_id_generate(log_);
        os_signpost_interval_begin(log_, spid_, "PowerGrid", "%s", name);
    }
    ~SignpostTracer() {
        os_signpost_interval_end(log_, spid_, "PowerGrid", "%s", name_);
    }
private:
    os_log_t log_;
    os_signpost_id_t spid_;
    const char* name_;
};
#define RANGE(...) SignpostTracer PG_CAT(_sp_, __LINE__){__VA_ARGS__};

#else
// ---- No profiling ----
#define RANGE(...)
#endif

#endif
