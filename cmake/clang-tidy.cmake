# cmake/clang-tidy.cmake
# Optional clang-tidy integration. Default OFF to preserve NVHPC compatibility.

option(ENABLE_CLANG_TIDY "Enable clang-tidy static analysis (default OFF)" OFF)

if(ENABLE_CLANG_TIDY)
    find_program(CLANG_TIDY_EXE NAMES clang-tidy clang-tidy-17 clang-tidy-16)
    if(CLANG_TIDY_EXE)
        message(STATUS "clang-tidy found: ${CLANG_TIDY_EXE}")
        set(CMAKE_CXX_CLANG_TIDY "${CLANG_TIDY_EXE}")
    else()
        message(WARNING "ENABLE_CLANG_TIDY=ON but clang-tidy not found in PATH")
    endif()
endif()
