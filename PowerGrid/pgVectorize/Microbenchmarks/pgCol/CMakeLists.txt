cmake_minimum_required(VERSION 3.2) #CMake version 3.2 required to properly support CMAKE_CXX_STANDARD

## Flags and settings that change the compiler (i.e. from default to C++11 or C++14)
## need to be set BEFORE calling project command (which does introspection and sets
## variables about the behavior of the compiler).
set(CMAKE_CXX_STANDARD 11)
#set(CMAKE_CXX_EXTENSIONS ON)

project("PowerGridpgColBench" CXX)

# Find Armadillo
find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})
# set(LIBS ${LIBS} ${ARMADILLO_LIBRARIES})

set(TEST_SOURCES
    ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp
)

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fast -fpermissive")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
add_definitions(-DARMA_NO_DEBUG)
# Make test executable
add_executable(pgCol ${TEST_SOURCES})
#target_link_libraries(tests Catch)