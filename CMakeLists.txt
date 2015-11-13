cmake_minimum_required (VERSION 2.8)

project (CHMprj)

# The version number.
set (CHM_VERSION_MAJOR 0)
set (CHM_VERSION_MINOR 1)

# Options. Turn on with 'cmake -Dmyvarname=ON'.
option(USE_OMP "Enable OpenMP" ON)
option(PV_FILTER "Build Paraview time filter" ON)
option(test "Build all tests." ON) # Makes boolean 'test' available.
option(MATLAB "Enable Matlab linkage" OFF)

macro(remove_cxx_flag flag)
    string(REPLACE "${flag}" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endmacro()

#Setup output paths
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH  ${CMAKE_SOURCE_DIR}/lib)

#include addition findXXX.cmake files
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/")

#Setup CXX flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -D_FORTIFY_SOURCE=0")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 ")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 ")
endif()


find_package(CGAL REQUIRED)
if(CGAL_FOUND)
    message(STATUS "Found CGAL")
endif()
#as per http://cgal-discuss.949826.n4.nabble.com/CMake-and-flags-td949906.html
#but doesn't seem to work
set( CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE )
include(${CGAL_USE_FILE}) #as per https://www.cgal.org/releases.html release 4.2

#remove this flag as it destroys the build on -O1 and spams warnings
remove_cxx_flag("-Werror=format-security -Wp,-D_FORTIFY_SOURCE=2")

#C#####
#GAL requires strict rounding
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -fp-model strict -msse4 -finline ") #-frounding-math needed for cgal
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -frounding-math")
endif()
########

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++ -Wno-c++11-extensions")
  set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -stdlib=libc++")
endif()


find_package(OpenMP)

#if(OPENMP_FOUND)
#    if(USE_OMP)
#        if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") #Clang doesn't have mainline support for openmp
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#        endif()
#    endif()
#endif()


if(MATLAB)
#for matlab to run you need csh installed
find_package(Matlab REQUIRED)
endif()

find_package(Armadillo REQUIRED)  #at least v3

find_package(Boost
             1.54.0
             COMPONENTS
             system
             filesystem
             date_time
             log
             thread
             regex
             program_options
             REQUIRED)

find_package(MeteoIO REQUIRED)

find_package(sundials REQUIRED)

#need this for calling the logger from multiple modules
ADD_DEFINITIONS(-DBOOST_LOG_DYN_LINK)

find_package (Threads REQUIRED)

find_package(TBB REQUIRED)
#SET(ENV{LD_PRELOAD} "libtbbmalloc_proxy.so.2")

#find_package(EIGEN3 REQUIRED)

find_package(GSL REQUIRED)

#can't use both, so if we are building the plugin we need to compile against paraview (which brings in VTK)
#otherwise use just VTK
if(PV_FILTER)
    FIND_PACKAGE(ParaView REQUIRED)
    INCLUDE(${PARAVIEW_USE_FILE})
else()
    find_package(VTK REQUIRED)
    include(${VTK_USE_FILE})
endif()

if (test)
    message(STATUS "Tests enabled. Run with make check")
    enable_testing()
    ADD_SUBDIRECTORY (tests)
endif()

# add a target to generate API documentation with Doxygen
# generate this via 'make doc'
find_package(Doxygen)
if(DOXYGEN_FOUND)
	configure_file(${CMAKE_SOURCE_DIR}/Doxyfile.in ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)

	add_custom_target(doc #ALL
		${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile
		WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
		COMMENT "Generating API documentation with Doxygen" VERBATIM 
	) 

endif(DOXYGEN_FOUND)

find_program(GVPR gvpr)
find_program(DOT dot)

if(NOT GVPR OR NOT DOT)
    message("Could not find gvpr and/or dot, module output to pdf will not function.")
else()
    message("gvpr and dot present, module output to pdf will function.")
endif()

add_subdirectory(src)

if(PV_FILTER)
    add_subdirectory(pv_filter)
endif()