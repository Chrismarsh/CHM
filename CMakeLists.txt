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
option(USE_OCL "Enable OpenCL with ViennaCL"  OFF )
option(PV_FILTER "Build Paraview time filter"  OFF )
option(test "Build all tests."  OFF ) # Makes boolean 'test' available.
option(MATLAB "Enable Matlab linkage"  OFF )
option(STATIC_ANLAYSIS "Enable PVS static anlaysis" OFF)

option(USE_TCMALLOC "Use tcmalloc from gperftools " ON)

set(SUPERBUILD FALSE CACHE BOOL "Do superbuild")

set(ENABLE_SAFE_CHECKS FALSE CACHE BOOL "Enable variable map checking. Runtime perf cost. Enable to debug")

#enable with cmake -DSUPERBUILD=TRUE .
if(SUPERBUILD)
    set(BUILD_ARMADILLO TRUE CACHE BOOL "Build Armadillo")
    set(BUILD_CGAL TRUE CACHE BOOL "Build CGAL")
    set(BUILD_GDAL TRUE  CACHE BOOL "Build GDAL")
    set(BUILD_GSL TRUE  CACHE BOOL "Build GSL")
    set(BUILD_VIENNACL TRUE  CACHE BOOL "Build ViennaCL")
    set(BUILD_BOOST TRUE  CACHE BOOL "Build Boost")
    set(BUILD_PROJ4 TRUE  CACHE BOOL "Build proj4")
    set(BUILD_NETCDF TRUE CACHE BOOL "Build netCDF")
    set(BUILD_VTK TRUE CACHE BOOL "Build VTK")
endif()

if(BUILD_ARMADILLO)
    set(Armadillo_DIR "${CMAKE_BINARY_DIR}/lib/armadillo/share/Armadillo/CMake")
endif()

if(BUILD_CGAL)
    set(CGAL_DIR "${CMAKE_BINARY_DIR}/lib/CGAL/lib64/CGAL")

        #this is defined in the BUILD_CGAL step
  #  LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_BINARY_DIR}/lib/CGAL/lib64/CGAL")
endif()

if(BUILD_GDAL)
    set(ENV{GDAL_ROOT} "${CMAKE_BINARY_DIR}/lib/gdal")  #this one really does have to be ENV
    set(BUILD_PROJ4 TRUE) #build our own proj4 if we are building gdal
endif()

if(BUILD_GSL)
    set(GSL_ROOT "${CMAKE_BINARY_DIR}/lib/gsl/bin")
endif()

if(BUILD_VIENNACL)
    set(VCL_ROOT "${CMAKE_BINARY_DIR}/lib/ViennaCL")
endif()

#To avoid issues #104, newer VTK version needs to be used. Until VTK7+ becomes available in the LTS/&c releases, allow for easily building a new VTK
if(BUILD_VTK)
    set(BUILD_PROJ4 TRUE)
    set(VTK_DIR "${CMAKE_BINARY_DIR}/lib/VTK/lib/cmake/vtk-8.0")
endif()

if(BUILD_BOOST)
    set(Boost_NO_SYSTEM_PATHS TRUE)
    set(BOOST_ROOT "${CMAKE_BINARY_DIR}/lib/boost")
endif()

if(BUILD_PROJ4)
    set(PROJ4_ROOT "${CMAKE_BINARY_DIR}/lib/proj4")
endif()

if(BUILD_NETCDF)
    set(NETCDF_DIR "${CMAKE_BINARY_DIR}/lib/netcdf-c")
    set(NETCDF_INCLUDE_DIR "${NETCDF_DIR}/include")
    set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${NETCDF_DIR}")
    set(ENV{HDF5_ROOT} "${CMAKE_BINARY_DIR}/lib/hdf5")

    if(NOT BUILD_GDAL)
        message(WARNING "You should do a GDAL build to avoid linking against system netCDF and HDF5.")
    endif()

endif()
#set(VCL_ROOT "${CMAKE_BINARY_DIR}/lib/ViennaCL")

## Any custom boost path.
#set(Boost_NO_SYSTEM_PATHS TRUE)
#set(BOOST_ROOT "/opt/boost_1_62_0")


#    Super build doesn't do TBB at the moment. If a custom location is needed, specify here
#set(TBB_INSTALL_DIR /cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/tbb/2017.2.132/tbb)

LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake/")

#########
# Helper macro that removes cxx flag flags
macro(remove_cxx_flag flag)
    string(REPLACE "${flag}" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endmacro()

configure_file(${CMAKE_SOURCE_DIR}/gdal-2.2.1.patch ${CMAKE_BINARY_DIR}/gdal-2.2.1.patch COPYONLY)
configure_file(${CMAKE_SOURCE_DIR}/vtk-8.0.0.patch ${CMAKE_BINARY_DIR}/vtk-8.0.0.patch COPYONLY)

#things covered by the superbuild
# super build things we are missing
configure_file(CMakeLists_external.txt.in
    ${CMAKE_BINARY_DIR}/lib/CMakeLists.txt)


#Trigger superbuild
#If superbuild is not needed, we still need to build sparsehash, gperftools, and eigen3
execute_process(COMMAND ${CMAKE_COMMAND}
        -DBUILD_VTK=${BUILD_VTK}
        -DBUILD_ARMADILLO=${BUILD_ARMADILLO}
        -DBUILD_CGAL=${BUILD_CGAL}
        -DBUILD_GDAL=${BUILD_GDAL}
        -DBUILD_GSL=${BUILD_GSL}
        -DBUILD_VIENNACL=${BUILD_VIENNACL}
        -DBUILD_BOOST=${BUILD_BOOST}
        -DBUILD_PROJ4=${BUILD_PROJ4}
        -DBUILD_NETCDF=${BUILD_NETCDF}

        -DUSE_TCMALLOC=${USE_TCMALLOC} -G "${CMAKE_GENERATOR}" .
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
execute_process(COMMAND ${CMAKE_COMMAND}  --build .
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/lib )


# lovely CMake script to integrate git hashes
# http://xit0.org/2013/04/cmake-use-git-branch-and-commit-details-in-project/
# Get the current working branch
# Generate gitrevision.hh if Git is available
# and the .git directory is present
# this is the case when the software is checked out from a Git repo
find_program(GIT_SCM git DOC "Git version control")
mark_as_advanced(GIT_SCM)
find_file(GITDIR NAMES .git PATHS ${CMAKE_SOURCE_DIR} NO_DEFAULT_PATH)
execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${GITDIR}
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
)

# Get the latest abbreviated commit hash of the working branch
execute_process(
        COMMAND git log -1 --format=%h
        WORKING_DIRECTORY ${GITDIR}
        OUTPUT_VARIABLE GIT_COMMIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
)


configure_file(
        src/version.h.in
        src/version.h
)


#ignore these two under Clion as CGAL will complain
if(CMAKE_BUILD_TYPE MATCHES RelWithDebInfo OR
        CMAKE_BUILD_TYPE MATCHES MinSizeRel OR
        NOT CMAKE_BUILD_TYPE)

    set(CMAKE_BUILD_TYPE "Release")
endif()


#Setup output paths

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin/${CMAKE_BUILD_TYPE})

#as per http://cgal-discuss.949826.n4.nabble.com/CMake-and-flags-td949906.html
#don't override internal settings
set( CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE )
if(BUILD_CGAL)
    find_package(CGAL REQUIRED HINTS ${CGAL_DIR})
else()
    find_package(CGAL REQUIRED)
endif()

if(CGAL_FOUND)
    message(STATUS "Found CGAL ")
    include(${CGAL_USE_FILE}) #as per https://www.cgal.org/releases.html release 4.2
endif()

#remove this flag as it destroys the build on -O1 and spams warnings, #thankscgal
remove_cxx_flag("-Werror=format-security -Wp,-D_FORTIFY_SOURCE=2")

remove_cxx_flag("-O2")
remove_cxx_flag("-O3")


#reset these back
if (CMAKE_BUILD_TYPE MATCHES Debug)

    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "4.7")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3 -O0 -fno-inline") #new to 4.8 https://gcc.gnu.org/gcc-4.8/changes.html
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3 -O0 -fno-inline -ggdb")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSAFE_CHECKS")
else()

    if(ENABLE_SAFE_CHECKS)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSAFE_CHECKS")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 ")
endif()

#Setup CXX flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -D_FORTIFY_SOURCE=0 -Wno-unused-variable -Wno-unknown-pragmas")


#GAL requires strict rounding
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qoverride-limits -fp-model strict -msse4 -finline ") #-frounding-math needed for cgal
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") # OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -frounding-math")
endif()
########

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++ -Wno-c++11-extensions")
  set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -stdlib=libc++")
endif()



if(USE_OMP)
    find_package(OpenMP)
    if(OPENMP_FOUND AND CMAKE_BUILD_TYPE MATCHES Release AND NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") #osx Clang doesn't have mainline support for openmp
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        add_definitions(-DVIENNACL_WITH_OPENMP)
    endif()
endif()


if(MATLAB)
    #for matlab to run you need csh installed
    find_package(Matlab REQUIRED)
endif()


find_package(OpenCL)
if(USE_OCL AND OPENCL_FOUND)
    message(STATUS "Enabling ViennaCL with OpenCL, vcl_scalar_type=float")

    add_definitions(-Dvcl_scalar_type=float) #This sets the precision of the floating point precision used for ViennaCL matrixes.
    add_definitions(-DVIENNACL_WITH_OPENCL)
    if(APPLE)
       set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -framework OpenCL")
    endif()
else()
    add_definitions(-Dvcl_scalar_type=float) #might as well default to double precision if no OpenCL device found.
    unset(OPENCL_LIBRARIES CACHE) ##zero these out so we can use them later without trouble
    unset(OPENCL_INCLUDE_DIRS CACHE)
endif()


if(BUILD_ARMADILLO)
    find_package(Armadillo REQUIRED  HINTS ${Armadillo_DIR} )
else()
    find_package(Armadillo REQUIRED)
endif()


set(EIGEN3_INCLUDE_DIR "${CMAKE_BINARY_DIR}/lib/")
message("${CMAKE_BINARY_DIR}/lib/eigen3/include/")

set( SPARSEHASH_ROOT "${CMAKE_BINARY_DIR}/lib/sparsehash")
find_package(Sparsehash)

#sparsehash doesn't compile with gcc 4.6.x - 4.8.x because of this
# https://bugs.launchpad.net/ubuntu/+source/gcc-4.6/+bug/1269803
# so we need to optionally disable it
if( SPARSEHASH_FOUND )
    add_definitions(-DUSE_SPARSEHASH)
else()
    set(SPARSEHASH_INCLUDE_DIRS "")
endif()

set(Tcmalloc_INCLUDE_DIR "${CMAKE_BINARY_DIR}/lib/gperftools/include")
find_package(tcmalloc)
if(Tcmalloc_FOUND)  #gperftools may not compile on machines w/o nanosleep so we need to optionall disable if it fails to compile
    message(STATUS "Found Tcmalloc, disabling builtin malloc, free")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free")
else()
    set(Tcmalloc_INCLUDE_DIR "")
    set(Tcmalloc_LIBRARY "")
endif()

if(BUILD_VIENNACL)
    find_package(ViennaCL 1.7.1 REQUIRED HINTS ${VCL_ROOT})
else()
    find_package(ViennaCL 1.7.1 REQUIRED)
endif()

message(STATUS "Found ViennaCL at ${VIENNACL_INCLUDE_DIR}")

find_package(Boost
             1.60.0
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
message(STATUS "Boost libs to link: ${Boost_LIBRARIES}")
message(STATUS "Found boost lib: ${Boost_LIBRARY_DIRS}")
#need this for calling the logger from multiple modules
ADD_DEFINITIONS(-DBOOST_LOG_DYN_LINK)

find_package(PROJ4 REQUIRED) #gdal dependency, only need to link against it
find_package(GDAL 2.1 REQUIRED )

find_package (Threads REQUIRED)

find_package(TBB REQUIRED)

find_package(GSL REQUIRED)

find_package(Curses REQUIRED)


set(NETCDF_CXX "YES")
find_package(NetCDF REQUIRED)

message(STATUS "netcdf lib ${NETCDF_LIBRARIES}")

#if we are building the plugin we need to compile against paraview (which brings in VTK)
#otherwise use just VTK. You cannot link against both at the same time
if(PV_FILTER)
    FIND_PACKAGE(ParaView REQUIRED)
    INCLUDE(${PARAVIEW_USE_FILE})
else()

    #N.B: The find VTK doesn't allow for compatibility between major versions. That is, if we ask for v7, v8 won't be allowed.
    #Therefore we need to ensure this version is up-to-date with what we are building. There are likely issues with VTK6 as per issue #104
    # so for now, we will accept finding *any* VTK version, but will flag it and error if we find <=6
    if(BUILD_VTK)
        find_package(VTK REQUIRED HINTS ${VTK_DIR})
    else()
        find_package(VTK REQUIRED)
    endif()
    include(${VTK_USE_FILE})

    if(VTK_MAJOR_VERSION LESS 7)
        #https://github.com/Chrismarsh/CHM/issues/104
        message( FATAL_ERROR "VTK > 6 is required, please see issue #104. You can build a compatible VTK using -DBUILD_VTK=ON." )
    endif()
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
endif()

if(PV_FILTER)
    add_subdirectory(pv_filter)
endif()


message(STATUS "${CMAKE_CXX_FLAGS}")