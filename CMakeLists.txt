cmake_minimum_required(VERSION 3.17) #CMake version 3.2 required to properly support CMAKE_CXX_STANDARD

## Flags and settings that change the compiler (i.e. from default to C++11 or C++14)
## need to be set BEFORE calling project command (which does introspection and sets
## variables about the behavior of the compiler).
set(CMAKE_CXX_STANDARD 17)

project("PowerGrid" CXX)

option(OPENACC_GPU "Enable OpenACC support on NVIDIA GPU" ON)
option(OPENACC_MP "Enable OpenACC support on CPU (PGI Only) " OFF)
option(MPISupport "Enable MPI Support" ON)

if(OPENACC_GPU AND OPENACC_MP)
message(FATAL_ERROR "You cannot build with both OPENACC_GPU AND OPENACC_MP at the same time")
endif()

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# Find Boost
set(Boost_USE_STATIC_LIBS OFF)
# Boost version >= 1.70.0 have improved support for cmake built into the boost install
# Unfortunately, this support breaks when used with PGI compilers
# As a work around, ignore the generated cmake modules
#set(Boost_NO_BOOST_CMAKE ON)

if(MPISupport)
    find_package(MPI)
    find_package(Boost 1.43 REQUIRED COMPONENTS program_options serialization mpi)
else(MPISupport)
    find_package(Boost 1.43 REQUIRED COMPONENTS program_options serialization)
endif(MPISupport)

if(Boost_FOUND)
	set(LIBS ${LIBS} ${Boost_LIBRARIES})
	include_directories(${Boost_INCLUDE_DIRS})
endif()

#set(CUDAVER 11.5)
#set(PGI_YEAR $ENV{HPCSDK_RELEASE})
if(OPENACC_GPU) 
	find_package(CUDAToolkit)
	set(LIBS ${LIBS} CUDA::cufft)
	set(LIBS ${LIBS} CUDA::cudart)
	set(LIBS ${LIBS} CUDA::nvToolsExt)
	#find_path(cuFFT_INCLUDE_DIR cufft.h HINTS $ENV{NVHPC}/Linux_x86_64/$ENV{HPCSDK_VERSION}/math_libs/${CUDAVER}/include)
	#include_directories(${cuFFT_INCLUDE_DIR})
	#find_library(cuFFT_LIBRARIES cufft HINTS $ENV{NVHPC}/Linux_x86_64/$ENV{HPCSDK_VERSION}/math_libs/${CUDAVER}/lib64)
	#set(LIBS ${LIBS} ${cuFFT_LIBRARIES})
	#find_library(cudaRT_LIBRARIES cudart HINTS $ENV{NVHPC}/Linux_x86_64/$ENV{HPCSDK_VERSION}/cuda/${CUDAVER}/targets/x86_64-linux/lib /usr/local/cuda/lib64)
	#set(LIBS ${LIBS} ${cudaRT_LIBRARIES})
	#find_library(nvToolsExt_LIBRARIES nvToolsExt HINTS $ENV{NVHPC}/Linux_x86_64/$ENV{HPCSDK_VERSION}/cuda/${CUDAVER}/targets/x86_64-linux/lib /usr/local/cuda/lib64)
	#set(LIBS ${LIBS} ${nvToolsExt_LIBRARIES})
	find_path(FFTW_INCLUDE_DIR fftw3.h HINTS $ENV{FFTW_INC})
	find_library(FFTW_LIBRARIES NAMES fftw3 HINTS $ENV{FFTW_DIR} /usr/lib/x86_64-linux-gnu/)
	find_library(FFTWF_LIBRARIES NAMES fftw3f HINTS $ENV{FFTW_DIR} /usr/lib/x86_64-linux-gnu/)
	include_directories(${FFTW_INCLUDE_DIR})
	set(LIBS ${LIBS} ${FFTW_LIBRARIES})
	set(LIBS ${LIBS} ${FFTWF_LIBRARIES})
else()
	#find_package(FFTW REQUIRED)
	find_path(FFTW_INCLUDE_DIR fftw3.h HINTS $ENV{FFTW_INC})
    find_library(FFTW_LIBRARIES NAMES fftw3 HINTS $ENV{FFTW_DIR} /usr/lib/x86_64-linux-gnu/)
    find_library(FFTWF_LIBRARIES NAMES fftw3f HINTS $ENV{FFTW_DIR} /usr/lib/x86_64-linux-gnu/)
	include_directories(${FFTW_INCLUDE_DIR})
	set(LIBS ${LIBS} ${FFTW_LIBRARIES})
	set(LIBS ${LIBS} ${FFTWF_LIBRARIES})
endif()

# Find Armadillo
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
message(${ARMADILLO_INCLUDE_DIRS})
set(LIBS ${LIBS} ${ARMADILLO_LIBRARIES})


find_package(HDF5 COMPONENTS CXX)
include_directories(${HDF5_INCLUDE_DIRS})
set(LIBS ${LIBS} ${HDF5_C_LIBRARIES})
set(LIBS ${LIBS} ${HDF5_CXX_LIBRARIES})

include_directories(Support/ArmaExtensions)

find_library(ISMRMRD_LIBRARIES ismrmrd HINTS /usr/lib/)
set(LIBS ${LIBS} ${ISMRMRD_LIBRARIES})

if(MPISupport)
	include_directories(MPI)
endif(MPISupport)

message(STATUS "Libraries to be linked in: ${LIBS}")

if(APPLE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -fno-math-errno -fno-trapping-math -flto")

	# ---------------------------------------------------------------------------
	# Apple Metal GPU compute (Apple Silicon / MTLGPUFamilyApple6+)
	# ---------------------------------------------------------------------------
	option(METAL_COMPUTE "Enable Apple Metal GPU compute (Apple Silicon only)" ON)
	if(METAL_COMPUTE)
		include(FindMetal)
		if(METAL_FOUND)
			# if constexpr (used in Gnufft.cpp for T1==float guard) requires C++17
			set(CMAKE_CXX_STANDARD 17)
			add_definitions(-DMETAL_COMPUTE)
			find_library(FOUNDATION_LIBRARY Foundation REQUIRED)
			set(LIBS ${LIBS} ${METAL_LIBRARY} ${ACCELERATE_LIBRARY} ${FOUNDATION_LIBRARY})

			# Compile .metal source → .air → .metallib → C byte array (xxd embed)
			set(METAL_SRC  ${CMAKE_SOURCE_DIR}/PowerGrid/Metal/gridding_metal.metal)
			set(METAL_AIR  ${CMAKE_BINARY_DIR}/PowerGrid/Metal/gridding.air)
			set(METALLIB   ${CMAKE_BINARY_DIR}/PowerGrid/Metal/gridding.metallib)
			set(METALEMBED ${CMAKE_BINARY_DIR}/PowerGrid/Metal/gridding_metallib.h)

			# Pre-create the output directory at configure time so the
			# WORKING_DIRECTORY cd in the custom command below always succeeds.
			file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/PowerGrid/Metal)

			add_custom_command(
				OUTPUT ${METALEMBED}
				COMMAND ${METAL_COMPILER} -O3 -o ${METAL_AIR} -c ${METAL_SRC}
				COMMAND xcrun -sdk macosx metallib -o ${METALLIB} ${METAL_AIR}
				COMMAND ${CMAKE_COMMAND}
					-DINPUT=${METALLIB}
					-DOUTPUT=${METALEMBED}
					-DVAR_NAME=gridding_metallib
					-P ${CMAKE_SOURCE_DIR}/cmake/embed_binary.cmake
				WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/PowerGrid/Metal
				DEPENDS ${METAL_SRC}
				COMMENT "Compiling and embedding Metal shaders")

			add_custom_target(MetalShaders ALL DEPENDS ${METALEMBED})
			include_directories(${CMAKE_BINARY_DIR}/PowerGrid/Metal)

			# Vector ops Metal shaders (separate metallib)
			set(VECOPS_SRC    ${CMAKE_SOURCE_DIR}/PowerGrid/Metal/vectorops_metal.metal)
			set(VECOPS_AIR    ${CMAKE_BINARY_DIR}/PowerGrid/Metal/vectorops.air)
			set(VECOPS_LIB    ${CMAKE_BINARY_DIR}/PowerGrid/Metal/vectorops.metallib)
			set(VECOPS_EMBED  ${CMAKE_BINARY_DIR}/PowerGrid/Metal/vectorops_metallib.h)

			add_custom_command(
				OUTPUT ${VECOPS_EMBED}
				COMMAND ${METAL_COMPILER} -O3 -o ${VECOPS_AIR} -c ${VECOPS_SRC}
				COMMAND xcrun -sdk macosx metallib -o ${VECOPS_LIB} ${VECOPS_AIR}
				COMMAND ${CMAKE_COMMAND}
					-DINPUT=${VECOPS_LIB}
					-DOUTPUT=${VECOPS_EMBED}
					-DVAR_NAME=vectorops_metallib
					-P ${CMAKE_SOURCE_DIR}/cmake/embed_binary.cmake
				WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/PowerGrid/Metal
				DEPENDS ${VECOPS_SRC}
				COMMENT "Compiling and embedding Metal vector ops shaders")

			add_custom_target(MetalVecOpsShaders ALL DEPENDS ${VECOPS_EMBED})
		else()
			message(WARNING "METAL_COMPUTE=ON but Metal toolchain not found; disabling Metal compute.")
			set(METAL_COMPUTE OFF)
		endif()
	endif()
else()
	if(CMAKE_COMPILER_IS_GNUCXX)

	else()
	add_definitions(-DBOOST_NO_ALIGNMENT)
		if(OPENACC_GPU)
			add_definitions(-DOPENACC_GPU)
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -cuda -Minfo=accel -fast -acc -gpu=cuda12.6,cc60,cc70,cc75,cc86,cc89 ")
		elseif(OPENACC_MP)
			add_definitions(-DOPENACC_MP)
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Mnollvm -O1 -Minfo=accel -acc  -ta=multicore -pgf90libs")
		else()
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O1 -Minfo=opt -pgf90libs")
		endif()

    endif()
endif()

add_definitions(-DARMA_NO_DEBUG)
add_definitions(-DARMA_DONT_USE_HDF5)
#add_definitions(-DARMA_DONT_USE_OPENMP)
if(OPENACC_GPU)
	add_definitions(-DUSE_NVTX)
endif()

add_subdirectory(PowerGrid)

add_subdirectory(apps)

# ---------------------------------------------------------------------------
# Metal unit test executable (built only when METAL_COMPUTE is ON)
# ---------------------------------------------------------------------------
if(METAL_COMPUTE AND METAL_FOUND)
    set(CATCH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/Support/Catch)
    add_executable(metal_tests
        ${CMAKE_SOURCE_DIR}/PowerGrid/Tests/main.cpp
        ${CMAKE_SOURCE_DIR}/PowerGrid/Tests/MetalNufftTests.cpp
        ${CMAKE_SOURCE_DIR}/PowerGrid/Tests/MetalVectorOpsTests.cpp
        ${CMAKE_SOURCE_DIR}/PowerGrid/Tests/pgComplexTests.cpp
        ${CMAKE_SOURCE_DIR}/PowerGrid/Tests/pgColTests.cpp
        ${CMAKE_SOURCE_DIR}/PowerGrid/Tests/pgColCplxTests.cpp
        ${CMAKE_SOURCE_DIR}/PowerGrid/Tests/pgMatTests.cpp
        ${CMAKE_SOURCE_DIR}/PowerGrid/Tests/pgPhase2Tests.cpp
        ${CMAKE_SOURCE_DIR}/PowerGrid/Tests/SyntheticReconTests.cpp
        ${CMAKE_SOURCE_DIR}/PowerGrid/Tests/ISMRMRDWriteTests.cpp)
    target_include_directories(metal_tests PRIVATE
        ${CATCH_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/PowerGrid
        ${CMAKE_SOURCE_DIR}/Support/ArmaExtensions)
    target_link_libraries(metal_tests PGCommon ${LIBS})
    # Disable LTO for test binary: some linked libs (FFTW SIMD) use inline
    # assembly patterns that LTO's back-end doesn't handle on Apple clang.
    target_compile_options(metal_tests PRIVATE -fno-lto)
    target_link_options(metal_tests PRIVATE -fno-lto)
    add_dependencies(metal_tests MetalShaders)

    # Metal throughput benchmark (standalone main, no Catch)
    add_executable(metal_bench
        ${CMAKE_SOURCE_DIR}/PowerGrid/Tests/MetalBenchmarks.cpp)
    target_include_directories(metal_bench PRIVATE
        ${CMAKE_SOURCE_DIR}/PowerGrid
        ${CMAKE_SOURCE_DIR}/Support/ArmaExtensions)
    target_link_libraries(metal_bench PGCommon ${LIBS})
    target_compile_options(metal_bench PRIVATE -fno-lto)
    target_link_options(metal_bench PRIVATE -fno-lto)
    add_dependencies(metal_bench MetalShaders)
endif()

# ---------------------------------------------------------------------------
# ISMRMRD integration tests (non-Apple builds, e.g. Docker/Linux)
# ---------------------------------------------------------------------------
if(NOT APPLE)
    set(CATCH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/Support/Catch)
    add_executable(ismrmrd_tests
        ${CMAKE_SOURCE_DIR}/PowerGrid/Tests/main.cpp
        ${CMAKE_SOURCE_DIR}/PowerGrid/Tests/ISMRMRDWriteTests.cpp)
    target_include_directories(ismrmrd_tests PRIVATE
        ${CATCH_INCLUDE_DIR}
        ${CMAKE_SOURCE_DIR}/PowerGrid
        ${CMAKE_SOURCE_DIR}/Support/ArmaExtensions)
    target_link_libraries(ismrmrd_tests PGCommon ${LIBS})
endif()

# CPU-only synthetic recon tests (all platforms)
set(CATCH_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/Support/Catch)
add_executable(cpu_tests
    ${CMAKE_SOURCE_DIR}/PowerGrid/Tests/main.cpp
    ${CMAKE_SOURCE_DIR}/PowerGrid/Tests/SyntheticReconTests.cpp)
target_include_directories(cpu_tests PRIVATE
    ${CATCH_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/PowerGrid
    ${CMAKE_SOURCE_DIR}/Support/ArmaExtensions)
target_link_libraries(cpu_tests PGCommon ${LIBS})

