cmake_minimum_required (VERSION 2.8)

project (CHMprj)

#http://stackoverflow.com/questions/31561309/cmake-warnings-under-os-x-macosx-rpath-is-not-specified-for-the-following-targe
#https://cmake.org/cmake/help/v3.0/policy/CMP0042.html
set(CMAKE_MACOSX_RPATH 1)

# The version number.
set (CHM_VERSION_MAJOR 0)
set (CHM_VERSION_MINOR 1)

#on gcc 4.6 or 4.8 sparsehash and tcmalloc do not build so can turn off with
#cmake -DCMAKE_DISABLE_FIND_PACKAGE_tcmalloc=TRUE -DCMAKE_DISABLE_FIND_PACKAGE_Sparsehash=TRUE .

# Options. Turn on with 'cmake -Dmyvarname=ON'.
option(USE_OMP "Enable OpenMP"  ON )
option(PV_FILTER "Build Paraview time filter"  OFF )
option(test "Build all tests."  OFF ) # Makes boolean 'test' available.
option(MATLAB "Enable Matlab linkage"  OFF )
option(STATIC_ANLAYSIS "Enable PVS static anlaysis" OFF)

set(Armadillo_DIR "${CMAKE_BINARY_DIR}/lib/armadillo/share/Armadillo/CMake")
set(CGAL_DIR "${CMAKE_BINARY_DIR}/lib/CGAL/lib64/CGAL")
set(ENV{GDAL_DIR} "${CMAKE_BINARY_DIR}/lib/gdal")
set(GSL_ROOT "${CMAKE_BINARY_DIR}/lib/gsl/bin")
set(VTK_DIR "${CMAKE_BINARY_DIR}/lib/VTK/lib/cmake/vtk-8.0")
set(VCL_ROOT "${CMAKE_BINARY_DIR}/lib/ViennaCL")
set(TBB_INSTALL_DIR /cvmfs/soft.computecanada.ca/easybuild/software/2017/Core/tbb/2017.2.132/tbb)
# TBB_INSTALL_DIR or $ENV{TBB21_INSTALL_DIR} or $ENV{TBB_INSTALL_DIR}
# TBB_INCLUDE_DIRS, where to find task_scheduler_init.h, etc.
#the mio and snowpack includes behave odly if you don't do the Cmake_source_dir like below, so treat these differently than the rest.
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCEDIR}/CMake/")

LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_BINARY_DIR}/lib/CGAL/lib64/CGAL")
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake/")
#some machines need to be told where odd system libraries are, so put them here
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake_osx/")
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake_fedora/")

#for meton
#set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/CMake_meton/")

## Any custom boost path. Few of the hpc systems have out of date boost libs
set(Boost_NO_SYSTEM_PATHS TRUE)
if (Boost_NO_SYSTEM_PATHS)
#    set(BOOST_ROOT "/opt/boost/boost_1_61_0")
#    set(BOOST_ROOT "/home/cbm038/boost_1_63_0")
#    set(BOOST_ROOT "/usr/local/Cellar/boost/1.62.0")
set(BOOST_ROOT "${CMAKE_BINARY_DIR}/lib/boost")
endif()


#########
# Helper macro that removes cxx flag flags
macro(remove_cxx_flag flag)
    string(REPLACE "${flag}" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endmacro()


#things covered by the superbuild
# super build things we are missing
configure_file(CMakeLists_external.txt.in
        ${CMAKE_BINARY_DIR}/lib/CMakeLists.txt)

# lovely CMake script to integrate git hashes
# http://xit0.org/2013/04/cmake-use-git-branch-and-commit-details-in-project/
# Get the current working branch

execute_process(COMMAND ${CMAKE_COMMAND} -G "${CMAKE_GENERATOR}" .
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/lib )
execute_process(COMMAND ${CMAKE_COMMAND} --build .
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/lib )


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
        mesher/version.h.in
        mesher/version.h
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
find_package(CGAL REQUIRED)
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
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -Og ") #new to 4.8 https://gcc.gnu.org/gcc-4.8/changes.html
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 ")
    endif()
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DDEBUG")
else()

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

if(OPENMP_FOUND)
    message(STATUS "Enabling ViennaCL with OpenCL")
    add_definitions(-DVIENNACL_WITH_OPENMP)
endif()

find_package(OpenCL)
if(OpenCL_FOUND)
    message(STATUS "Enabling ViennaCL with OpenCL")
    add_definitions(-Dvcl_scalar_type=float) #This sets the precision of the floating point precision used for ViennaCL matrixes.
    add_definitions(-DVIENNACL_WITH_OPENCL)
        if(APPLE)
             set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -framework OpenCL")
        endif()
else()
    add_definitions(-Dvcl_scalar_type=double) #might as well default to double precision if no OpenCL device found.
    set(OPENCL_LIBRARIES "")
endif()

find_package(Armadillo REQUIRED  PATH ${CMAKE_BINARY_DIR}/lib/armadillo/share/Armadillo/CMake NO_DEFAULT_PATH )

set(EIGEN3_INCLUDE_DIR "${CMAKE_BINARY_DIR}/lib/")
message("${CMAKE_BINARY_DIR}/lib/eigen3/include/")
#find_package(Eigen3 REQUIRED)
#if(NOT EIGEN3_INCLUDE_DIR)
#    set(Eigen3_DIR "${CMAKE_BINARY_DIR}/lib/eigen3")
#    find_package(Eigen3 REQUIRED)
#endif()

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

find_package(tcmalloc)
if(Tcmalloc_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free")
    set(Tcmalloc_INCLUDE_DIR "${CMAKE_BINARY_DIR}/lib/gperftools/include")
else()

    set(Tcmalloc_INCLUDE_DIR "")
    set(Tcmalloc_LIBRARY "")
endif()
find_package(ViennaCL REQUIRED)

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

find_package(GDAL 2.1 REQUIRED)

find_package (Threads REQUIRED)

find_package(TBB REQUIRED)

find_package(GSL REQUIRED)

find_package(Curses REQUIRED)

#find_package(Sundials REQUIRED)

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

add_subdirectory(mesher)

