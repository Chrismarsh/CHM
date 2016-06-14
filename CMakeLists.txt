cmake_minimum_required (VERSION 2.8)

project (CHMprj)

#http://stackoverflow.com/questions/31561309/cmake-warnings-under-os-x-macosx-rpath-is-not-specified-for-the-following-targe
#https://cmake.org/cmake/help/v3.0/policy/CMP0042.html
set(CMAKE_MACOSX_RPATH 1)

# The version number.
set (CHM_VERSION_MAJOR 0)
set (CHM_VERSION_MINOR 1)


# Options. Turn on with 'cmake -Dmyvarname=ON'.
option(USE_OMP "Enable OpenMP"  ON )
option(PV_FILTER "Build Paraview time filter"  OFF )
option(test "Build all tests."  ON ) # Makes boolean 'test' available.
option(MATLAB "Enable Matlab linkage"  OFF )

#meteoio plugins
SET(PLUGIN_A3DIO ON CACHE BOOL "Compilation A3DIO ON or OFF")
SET(PLUGIN_ALPUG ON CACHE BOOL "Compilation ALPUG ON or OFF")
SET(PLUGIN_ARCIO ON CACHE BOOL "Compilation ARCIO ON or OFF")
SET(PLUGIN_ARPSIO ON CACHE BOOL "Compilation ARPSIO ON or OFF")
SET(PLUGIN_BORMAIO OFF CACHE BOOL "Compilation BORMAIO ON or OFF")
SET(PLUGIN_CNRMIO OFF CACHE BOOL "Compilation CNRMIO ON or OFF")
SET(PLUGIN_COSMOXMLIO OFF CACHE BOOL "Compilation COSMOXMLIO ON or OFF")
SET(PLUGIN_GEOTOPIO ON CACHE BOOL "Compilation GEOTOPIO ON or OFF")
SET(PLUGIN_GRASSIO ON CACHE BOOL "Compilation GRASSIO ON or OFF")
SET(PLUGIN_GRIBIO OFF CACHE BOOL "Compilation GRIBIO ON or OFF")
SET(PLUGIN_GSNIO OFF CACHE BOOL "Compilation GSNIO ON or OFF")
SET(PLUGIN_IMISIO OFF CACHE BOOL "Compilation IMISIO ON or OFF")
SET(PLUGIN_NETCDFIO OFF CACHE BOOL "Compilation NETCDFIO ON or OFF")
SET(PLUGIN_PGMIO ON CACHE BOOL "Compilation PGMIO ON or OFF")
SET(PLUGIN_PNGIO OFF CACHE BOOL "Compilation PNGIO ON or OFF")
SET(PLUGIN_PSQLIO OFF CACHE BOOL "Compilation PSQLIO ON or OFF")
SET(PLUGIN_SMETIO ON CACHE BOOL "Compilation SMETIO ON or OFF")
SET(PLUGIN_SNIO ON CACHE BOOL "Compilation SNIO ON or OFF")
SET(PLUGIN_SASEIO OFF CACHE BOOL "Compilation SASEIO ON or OFF")
SET(PROJ4 OFF CACHE BOOL "Use PROJ4 for the class MapProj ON or OFF")
SET(DATA_QA OFF CACHE BOOL "Data Quality Assurance outputs ON or OFF")

#the mio and snowpack includes behave odly if you don't do the Cmake_source_dir like below, so treat these differently than the rest.
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake/")

#some machines need to be told where odd system libraries are, so put them here
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake_osx/")

#for meton
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake_meton/")

## Any custom boost path. Few of the hpc systems have out of date boost libs
set(Boost_NO_SYSTEM_PATHS FALSE)
if (Boost_NO_SYSTEM_PATHS)
    set(BOOST_ROOT "/opt/boost/boost_1_61_0")
endif()
#set(Boost_NO_SYSTEM_PATHS FALSE)
#if (Boost_NO_SYSTEM_PATHS)
#    set(BOOST_ROOT "/home/cbm038/boost_1_58_0")
##    set(BOOST_LIBRARYDIR /home/cbm038/boost158/lib)
#
#endif()


#########
# Helper macro that removes cxx flag flags
macro(remove_cxx_flag flag)
    string(REPLACE "${flag}" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endmacro()

#Setup output paths
set(EXECUTABLE_OUTPUT_PATH ${CMAKE_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH  ${CMAKE_SOURCE_DIR}/lib)


#Setup CXX flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall ")

#ignore these two under Clion as CGAL will complain
if(CMAKE_BUILD_TYPE MATCHES RelWithDebInfo OR
        CMAKE_BUILD_TYPE MATCHES MinSizeRel OR
        NOT CMAKE_BUILD_TYPE)

    set(CMAKE_BUILD_TYPE "Release")
endif()



#as per http://cgal-discuss.949826.n4.nabble.com/CMake-and-flags-td949906.html
#don't override internal settings
set( CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE )
find_package(CGAL REQUIRED)
if(CGAL_FOUND)
    message(STATUS "Found CGAL ")
endif()
include(${CGAL_USE_FILE}) #as per https://www.cgal.org/releases.html release 4.2

#remove this flag as it destroys the build on -O1 and spams warnings, #thankscgal
remove_cxx_flag("-Werror=format-security -Wp,-D_FORTIFY_SOURCE=2")

remove_cxx_flag("-O2")
remove_cxx_flag("-O3")

#reset these back
if (CMAKE_BUILD_TYPE MATCHES Debug)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 ")
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 ")
endif()

#Setup CXX flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -D_FORTIFY_SOURCE=0 -Wno-unused-variable")


#GAL requires strict rounding
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -fp-model strict -msse4 -finline ") #-frounding-math needed for cgal
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") # OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -frounding-math")
endif()
########

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++ -Wno-c++11-extensions")
  set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -stdlib=libc++")
endif()


find_package(OpenMP)

if(OPENMP_FOUND AND CMAKE_BUILD_TYPE MATCHES Release)
    if(USE_OMP)
        if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") #osx Clang doesn't have mainline support for openmp
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        endif()
    endif()
endif()


if(MATLAB)
    #for matlab to run you need csh installed
    find_package(Matlab REQUIRED)
endif()

find_package(Armadillo REQUIRED)  #at least v3

find_package(Boost
             1.57.0
             COMPONENTS
             system
             filesystem
             date_time
             log
             thread
             regex
             program_options
             REQUIRED)
message(STATUS "Found boost incl: ${Boost_INCLUDE_DIR}")
message(STATUS "Found boost lib: ${Boost_LIBRARY_DIRS}")
#need this for calling the logger from multiple modules
ADD_DEFINITIONS(-DBOOST_LOG_DYN_LINK)


#find_package(MeteoIO REQUIRED)
#find_package(Libsnowpack REQUIRED)

find_package (Threads REQUIRED)

find_package(TBB REQUIRED)

find_package(GSL REQUIRED)

find_package(Curses REQUIRED)

#if we are building the plugin we need to compile against paraview (which brings in VTK)
#otherwise use just VTK. You cannot link against both at the same time
if(PV_FILTER)
    FIND_PACKAGE(ParaView REQUIRED)
    INCLUDE(${PARAVIEW_USE_FILE})
else()
    find_package(VTK REQUIRED)
    include(${VTK_USE_FILE})
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
    message(STATUS "Could not find gvpr and/or dot, module output to pdf will not function.")
else()
    message(STATUS "gvpr and dot present, module output to pdf will function.")
endif()


#setup src dirs
add_subdirectory(src)

if (test)
    message(STATUS "Tests enabled. Run with make check")
    enable_testing()
#    ADD_SUBDIRECTORY (src/tests)
endif()

if(PV_FILTER)
    add_subdirectory(pv_filter)
endif()

add_subdirectory(triangle)

