# Changxi Zheng (cxz@cs.columbia.edu)
# Copyright @ Jan 2012
project(TimeDepHJB)

cmake_minimum_required(VERSION 2.6)

# if we compile the code into 64bit app
option(USE_64BIT_COMPILE "Compile 64bit executable" ON)
option(USE_OPENMP "Turn on the OpenMP feature when compiling" OFF)
option(USE_DEBUG "Turn on the debug mode" OFF)
option(USE_TIME_PROFILING "Turn on timing" OFF)

if ( USE_64BIT_COMPILE )
    add_definitions(-m64)
    set(LINK_FLAGS -m64)
endif ( USE_64BIT_COMPILE )

if ( USE_OPENMP )
    add_definitions(-DUSE_OPENMP)
endif ( USE_OPENMP )

if ( USE_TIME_PROFILING )
    add_definitions(-DENABLE_TIME_PROFILING)
endif ( USE_TIME_PROFILING )
#===================================================================
## Compiler
# set compiler flags for debug/release
if ( USE_DEBUG )
    add_definitions(-DDEBUG)
    set(CMAKE_BUILD_TYPE Debug)
else ( USE_DEBUG )
    add_definitions(-DNDEBUG -Wno-deprecated)
    set(CMAKE_BUILD_TYPE Release)
endif ( USE_DEBUG )

#===================================================================
## Libraries
# check the pthread
find_package(Threads 1.0 REQUIRED)

# check boost
set(Boost_USE_MULTITHREAD OFF)
find_package(Boost 1.40 REQUIRED COMPONENTS "system" "program_options")

add_definitions(-Wall)

add_subdirectory(src)

