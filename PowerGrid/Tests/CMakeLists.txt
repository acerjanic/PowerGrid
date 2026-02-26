cmake_minimum_required(VERSION 3.14) #CMake 3.14 required for FetchContent

## Flags and settings that change the compiler (i.e. from default to C++11 or C++14)
## need to be set BEFORE calling project command (which does introspection and sets
## variables about the behavior of the compiler).
set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_EXTENSIONS ON)

project("PowerGridTests" VERSION 1.9.0 LANGUAGES CXX)

# Generate Version.h so PGLog.hpp can resolve POWERGRID_VERSION_STRING
# in the standalone test build (which doesn't use the root CMakeLists.txt).
configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/../Core/Version.h.in"
    "${CMAKE_CURRENT_BINARY_DIR}/Version.h"
    @ONLY
)
include_directories(${CMAKE_CURRENT_BINARY_DIR})

# Fetch spdlog and indicators for structured logging
include(FetchContent)
FetchContent_Declare(spdlog
    GIT_REPOSITORY https://github.com/gabime/spdlog.git
    GIT_TAG        v1.14.1)
FetchContent_Declare(indicators
    GIT_REPOSITORY https://github.com/p-ranav/indicators.git
    GIT_TAG        v2.3)
FetchContent_MakeAvailable(spdlog indicators)


# Prepare "Catch" library for other executables
set(CATCH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../Support/Catch)
add_library(Catch INTERFACE)
target_include_directories(Catch INTERFACE ${CATCH_INCLUDE_DIR})

# Find dependencies needed by utilityTests and operatorTests
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

find_path(FFTW_INCLUDE_DIR fftw3.h
    HINTS $ENV{FFTW_INC} /usr/include /usr/local/include)
find_library(FFTW_LIBRARIES  NAMES fftw3
    HINTS $ENV{FFTW_DIR} /usr/lib/x86_64-linux-gnu/ /usr/lib)
find_library(FFTWF_LIBRARIES NAMES fftw3f
    HINTS $ENV{FFTW_DIR} /usr/lib/x86_64-linux-gnu/ /usr/lib)
if(FFTW_INCLUDE_DIR)
    include_directories(${FFTW_INCLUDE_DIR})
endif()

# Support headers (ArmaExtensions, etc.)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../Support)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/../../Support/ArmaExtensions)

# Additional compiled sources for utilityTests
set(UTILITY_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/../Gridding/griddingSupport.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../FFT/fftCPU.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../Penalties/Robject.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../Penalties/QuadPenalty.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../Penalties/TVPenalty.cpp
)

# Operator sources needed for encoding operator tests
set(OPERATOR_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/../FFT/ftCpu.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../FFT/ftCpuWithGrads.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../Operators/Gdft.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../Gridding/gridding.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../Operators/Gnufft.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../Operators/GdftR2.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../Gridding/TimeSegmentation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/../Operators/SENSE.cpp
)

set(TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pgComplexTests.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pgColTests.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pgColCplxTests.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pgMatTests.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/pgPhase2Tests.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/penaltyTests.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/utilityTests.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/operatorTests.cpp
    ${UTILITY_SOURCES}
    ${OPERATOR_SOURCES}
)

# Only add PGI/NVHPC-specific flags when using nvc++ or pgc++
if(CMAKE_CXX_COMPILER_ID STREQUAL "PGI" OR CMAKE_CXX_COMPILER_ID STREQUAL "NVHPC")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -acc -gpu=cuda12.6 -fast")
else()
    # Suppress unknown pragma warnings from OpenACC directives when using g++/clang
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas")
endif()

# Make test executable
add_executable(tests ${TEST_SOURCES})
target_include_directories(tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_link_libraries(tests Catch ${ARMADILLO_LIBRARIES} ${FFTW_LIBRARIES} ${FFTWF_LIBRARIES} spdlog::spdlog indicators::indicators)
