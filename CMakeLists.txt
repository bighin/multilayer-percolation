cmake_minimum_required(VERSION 3.17)
project(bilayer C)

set(CMAKE_C_STANDARD 11)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIR})

add_executable(bilayer
        bonds.c
        bonds.h
        clusters.c
        clusters.h
        common.h
        percolation.c
        sites.c
        sites.h jumps.c jumps.h)

target_link_libraries(bilayer ${GSL_LIBRARIES})
