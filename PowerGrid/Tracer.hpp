
/*****************************************************************************

    File Name   [Tracer.hpp]

    Synopsis    [Tracer class for inserting nvtx range markers using RAII idiom. 
                Look at NVIDIA dev blog for more details.]

    Description []

    Revision    [0.1.0; Alex Cerjanic, BIOE UIUC]

    Date        [4/10/2018]

*****************************************************************************/

/// @file Tracer.hpp
/// @brief RAII NVTX range marker for GPU profiling with NVIDIA Nsight.

#ifndef PowerGrid_Tracer_hpp
#define PowerGrid_Tracer_hpp

#ifdef USE_NVTX
/// @brief RAII wrapper that pushes/pops an NVTX named range for GPU profiling.
///
/// Construct a `Tracer` at the start of a code region to open an NVTX range;
/// the destructor automatically closes it when the object goes out of scope.
/// Use the `RANGE(name)` macro for convenient single-line usage.
///
/// Only active when compiled with `USE_NVTX` defined (i.e., `OPENACC_GPU` builds).
/// In all other builds, `RANGE(name)` expands to nothing.
class Tracer {
public:
    /// @brief Open an NVTX range with the given name.
    ///
    /// @param name  Null-terminated string label shown in the profiler timeline.
    Tracer(const char* name)
    {
        nvtxRangePushA(name);
    }
    /// @brief Close the NVTX range.
    ~Tracer()
    {
        nvtxRangePop();
    }
};
#define RANGE(name) Tracer uniq_name_using_macros(name);
#else
#define RANGE(name)
#endif

#endif
