cmake_minimum_required(VERSION 2.8)
project("PowerGrid" CXX)


option(OPENACC "Enable OpenACC support" ON)
option(ISMRMD "Enable ISMRMRD Suppport" OFF)

#set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_STANDARD_REQUIRED 1)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(MATIO)
include_directories(${MATIO_INCLUDE_DIRS})
set(LIBS ${LIBS} ${MATIO_LIBRARIES})
if(OPENACC)
	find_path(cuFFT_INCLUDE_DIR cufft.h HINTS /usr/local/cuda/include)
	include_directories(${cuFFT_INCLUDE_DIR})
	find_library(cuFFT_LIBRARIES cufft HINTS /usr/local/cuda/lib64)
	set(LIBS ${LIBS} ${cuFFT_LIBRARIES})

else(OPENACC)
	#find_package(FFTW REQUIRED)
	find_path(FFTW_INCLUDE_DIR fftw3.h)
	find_library(FFTW_LIBRARIES NAMES fftw3)
	include_directories(${FFTW_INCLUDE_DIR})
	set(LIBS ${LIBS} ${FFT_LIBRARIES})
endif(OPENACC)
if(${USE_ACML})
	find_package(ACML)
	set(LIBS ${LIBS} ${ACML_LIBRARIES})
endif()

if (ISMRMRD)


endif (ISMRMRD)
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
MESSAGE(${ARMADILLO_INCLUDE_DIRS})
set(LIBS ${LIBS} ${ARMADILLO_LIBRARIES})
find_package(Xerces)
include_directories(${XERCESC_INCLUDE_DIRS})
set(LIBS ${LIBS} ${XERCESC_LIBRARY})
if(ISMRMRD)
    find_package(HDF5)
    include_directories(${HDF5_INCLUDE_DIRS})
endif(ISMRMRD)
set(LIBS ${LIBS} ${HDF5_C_LIBRARIES})
set(LIBS ${LIBS} ${HDF5_CXX_LIBRARIES})

set(POWERGRID_TEST_GDFT
    ./PowerGrid/PowerGridGdft.cpp
    ./PowerGrid/config.cxx
)

set(POWERGRID_TEST_GGRID
	./PowerGrid/PowerGridGgrid.cpp
	./PowerGrid/config.cxx
)
add_executable(PowerGridGdft ${POWERGRID_TEST_GDFT})
add_executable(PowerGridGgrid ${POWERGRID_TEST_GGRID})

target_link_libraries(PowerGridGdft ${LIBS})
target_link_libraries(PowerGridGgrid ${LIBS})

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/include")
if(APPLE)
	#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lgfortran")
else()
	if(CMAKE_COMPILER_IS_GNUCXX)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${FFT_LIBRARIES}")
	else()
		if(OPENACC)
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g  -Minfo -ta=tesla -fpermissive -pgf77libs")
		else()
			set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -pgf77libs")
		endif()
	endif()
endif()
#Needed for CMake 2.8 or so to not add the -rdynamic flag to compilers like PGI
#PGI chokes on link because it does not understand the -rdynamic flag
#Why CMake insists on passing this flag, I don't know...
#SET(CMAKE_SHARED_LIBRARY_LINK_C_FLAGS "")
#SET(CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS "${CMAKE_SHARED_LIBRARY_LINK_CXX_FLAGS} -lfftw3")

#include(CheckCXXCompilerFlag)
#CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
#CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
#if(COMPILER_SUPPORTS_CXX11)
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
#elseif(COMPILER_SUPPORTS_CXX0X)
#	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
#else()
#        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
#endif()
if(APPLE)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Ofast -flto")
endif(APPLE)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -lfftw3")
add_definitions(-D__STRICT_ANSI__)
add_definitions(-DARMA_DONT_USE_CXX11)
