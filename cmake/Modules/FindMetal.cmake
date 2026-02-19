# FindMetal.cmake
# Locate the Metal compiler and required frameworks on macOS.
#
# Imported targets:
#   (none — frameworks linked via METAL_LIBRARY and ACCELERATE_LIBRARY variables)
#
# Result variables:
#   METAL_FOUND           - True if xcrun metal + Metal.framework + Accelerate.framework found
#   METAL_COMPILER        - Full path to the 'metal' compiler executable
#   METAL_LIBRARY         - Metal.framework
#   ACCELERATE_LIBRARY    - Accelerate.framework

execute_process(
    COMMAND xcrun --sdk macosx --find metal
    OUTPUT_VARIABLE METAL_COMPILER
    ERROR_QUIET
    OUTPUT_STRIP_TRAILING_WHITESPACE)

find_library(METAL_LIBRARY Metal)
find_library(ACCELERATE_LIBRARY Accelerate)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Metal
    DEFAULT_MSG
    METAL_COMPILER
    METAL_LIBRARY
    ACCELERATE_LIBRARY)

mark_as_advanced(METAL_COMPILER METAL_LIBRARY ACCELERATE_LIBRARY)
