cmake_minimum_required(VERSION 3.1)
project("PowerGrid")
add_executable(PowerGrid ./PowerGrid/main.cpp)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")
find_package(MATIO)
find_package(Armadillo
	HINTS "${CMAKE_PREFIX_PATH}/software/armadillo-4.650.4/}")


