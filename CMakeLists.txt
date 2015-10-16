cmake_minimum_required (VERSION 2.8)


project (CHMprj)
# The version number.
set (CHM_VERSION_MAJOR 0)
set (CHM_VERSION_MINOR 1)

set(USE_OMP FALSE)

macro(remove_cxx_flag flag)
    string(REPLACE "${flag}" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endmacro()

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
set(CMAKE_BUILD_TYPE Debug)

#message(FATAL_ERROR "${CMAKE_CXX_COMPILER_ID}")

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -std=c++11 -Wall -D_FORTIFY_SOURCE=0")

# Options. Turn on with 'cmake -Dmyvarname=ON'.
option(test "Build all tests." ON) # Makes boolean 'test' available.


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



find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

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

if(OPENMP_FOUND)
    if(OMP)
        if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") #Clang doesn't have mainline support for openmp
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        endif()
    endif()
endif()

set(MATLAB false)
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



if (test)
    message(STATUS "Tests enabled. Run with make check")
    enable_testing()
    ADD_SUBDIRECTORY (tests)
endif()

add_subdirectory(src)

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

