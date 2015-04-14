cmake_minimum_required(VERSION 3.1)
project("PowerGrid")
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(MATIO)
include_directories(${MATIO_INCLUDE_DIRS})
set(LIBS ${LIBS} ${MATIO_LIBRARIES})
find_package(FFTW)
include_directories(${FFTW_INCLUDE_DIRS})
set(LIBS ${LIBS} ${FFTW_LIBRARIES})
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
MESSAGE(${ARMADILLO_INCLUDE_DIRS})
set(LIBS ${LIBS} ${ARMADILLO_LIBRARIES})


set(POWERGRID_ALL_SOURCES
    ./PowerGrid/main.cpp
    ./PowerGrid/test_pwls_pcg.cpp
)
add_executable(PowerGridMain ${POWERGRID_ALL_SOURCES})
target_link_libraries(PowerGridMain ${LIBS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/include")

if(CMAKE_COMPILER_IS_GNUCXX)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lgfortran")
endif()
 
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-std=c++11" COMPILER_SUPPORTS_CXX11)
CHECK_CXX_COMPILER_FLAG("-std=c++0x" COMPILER_SUPPORTS_CXX0X)
if(COMPILER_SUPPORTS_CXX11)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
elseif(COMPILER_SUPPORTS_CXX0X)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
else()
        message(STATUS "The compiler ${CMAKE_CXX_COMPILER} has no C++11 support. Please use a different C++ compiler.")
endif()
