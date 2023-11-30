cmake_minimum_required (VERSION 3.16)
project (CHM-project C CXX Fortran)

# In CMake 3.12 and above the find_package(<PackageName>) command now searches prefixes specified by the <PackageName>_ROOT
# CMake variable and the <PackageName>_ROOT environment variable.
#https://cmake.org/cmake/help/latest/policy/CMP0074.html
cmake_policy(SET CMP0074 NEW)

#https://cmake.org/cmake/help/latest/policy/CMP0042.html#policy:CMP0042
# set RPATH ON
cmake_policy(SET CMP0042 NEW)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
    set(CMAKE_FIND_FRAMEWORK LAST) # framework gdal, etc can end up superceeding other search paths
endif()


# Options. Turn on with 'cmake -Dmyvarname=ON'.
option(USE_MPI "Enable MPI"  OFF )
option(USE_OMP "Enable OpenMP"  ON )
option(OMP_SAFE_EXCEPTION "Enables safe exception handling from within OMP regions." ON)
option(ENABLE_SAFE_CHECKS "Enable variable map checking. Runtime perf cost. Allows for ensuring a variable is indeed available to be lookedup." ON)
option(BUILD_TESTS "Build all tests."  OFF ) # Makes boolean 'test' available.
option(STATIC_ANLAYSIS "Enable PVS static anlaysis" OFF)
option(USE_TCMALLOC "Use tcmalloc from gperftools " OFF)
option(USE_JEMALLOC "Use jemalloc" ON)
option(BUILD_DOCS "Builds documentation" OFF)


message(STATUS "This is an MPI build? ${USE_MPI}")



if(BUILD_WITH_CONAN AND USE_JEMALLOC)
    message(WARNING "With a Conan build, USE_JEMALLOC will use a system jemalloc" )
endif()


option (FORCE_COLORED_OUTPUT "Always produce ANSI-colored output (GNU/Clang only)." TRUE)
if (${FORCE_COLORED_OUTPUT})
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        add_compile_options (-fdiagnostics-color=always)
    elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
        add_compile_options (-fcolor-diagnostics)
    endif ()
endif ()

if(USE_JEMALLOC AND USE_TCMALLOC)
    message(FATAL_ERROR "Cannot have both USE_TCMALLOC=TRUE and USE_JEMALLOC=TRUE")
endif()

LIST(APPEND CMAKE_PREFIX_PATH "${CMAKE_BINARY_DIR}")
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_BINARY_DIR}")
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake/")
LIST(APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/CMake/")


#########
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

macro(create_target package)
    if(NOT TARGET ${package}::${package})
        add_library(${package}::${package} INTERFACE IMPORTED)

        set_target_properties(${package}::${package} PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                    "${${package}_INCLUDE_DIRS};${${package}_INCLUDE_DIR}")


        set_property(TARGET ${package}::${package} PROPERTY INTERFACE_LINK_LIBRARIES
                "${${package}_LIBRARIES};${${package}_LINKER_FLAGS_LIST}")
#        set_property(TARGET ${package}::${package} PROPERTY INTERFACE_COMPILE_DEFINITIONS
#                ${${package}_COMPILE_DEFINITIONS})
#        set_property(TARGET ${package}::${package} PROPERTY INTERFACE_COMPILE_OPTIONS
#                "${${package}_COMPILE_OPTIONS_LIST}")
    endif()
endmacro()

# ---------------------------------------------------------------
# ------ Start of find_package section --------------------------
# ---------------------------------------------------------------

# Need to check for MPI before building external libraries
find_package(MPI REQUIRED COMPONENTS C CXX Fortran) # MPI_C is needed for the mpi that HDF can bring in
if(MPI_FOUND AND USE_MPI)
    message(STATUS "Found MPI at ${MPI_CXX_INCLUDE_DIRS}")
    add_definitions(-DUSE_MPI)
else()
  message(WARNING "MPI not enabled. Building for single process.")
endif()

if(USE_MPI AND NOT MPI_FOUND)
    message(FATAL_ERROR "MPI was asked for, but not found")
endif()

set(Boost_USE_MULTITHREADED ON) #https://stackoverflow.com/a/58085634/410074
if(USE_MPI AND MPI_FOUND)
    find_package(Boost
            1.71.0
            COMPONENTS
            system
            filesystem
            date_time
            log
            thread
            regex
            iostreams
            program_options
            mpi
            serialization
            REQUIRED)
else()
    find_package(Boost
            1.71.0
            COMPONENTS
            system
            filesystem
            date_time
            log
            thread
            regex
            iostreams
            program_options
            REQUIRED
            )
endif()

#need this for calling the logger from multiple modules
add_definitions(-DBOOST_LOG_DYN_LINK)

# Note the boost::bind deprecation warnings for the global namespace:
#
#   ‘#pragma message: The practice of declaring the Bind placeholders (_1, _2, ...) in the global namespace is deprecated. Please use <boost/bind/bind.hpp> + using namespace boost::placeholders, or define BOOST_BIND_GLOBAL_PLACEHOLDERS to retain the current behavior.’
add_definitions(-DBOOST_BIND_GLOBAL_PLACEHOLDERS)


find_package(CGAL 5.0 REQUIRED)

# There is a bug in clang <11 (so far), that has a compiler parse error w/ boost 1.71
# https://bugs.llvm.org/show_bug.cgi?id=43266
# https://github.com/Oslandia/SFCGAL/issues/188
# https://stackoverflow.com/questions/57857572/why-does-boost-log-break-boost-gmp-multiprecision
# But CGAL uses boost/mp which appears to cause this problem. So if we see clang, disable boost MP usage.
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" )
    add_definitions(-DCGAL_DO_NOT_USE_BOOST_MP)

    message(WARNING "Disabling CGAL's use of boost mp to avoid clang compiler parse error.")
endif()


#these are the flags we will use to build CHM. Used with target properties instead of modifying CXX_FLAGS
set(CHM_BUILD_FLAGS "-Wall -Wno-unused-variable -Wno-unknown-pragmas")
#reset these back
if (CMAKE_BUILD_TYPE MATCHES Debug)
    message(WARNING "Debug mode")
    add_definitions(-DSAFE_CHECKS)

    set(DEBUG_FLAGS "")

    if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
        # A new general optimization level, -Og, has been introduced. It addresses the need for fast compilation and
        # a superior debugging experience while providing a reasonable level of run-time performance.
        # Overall experience for development should be better than the default optimization level -O0.
        # new to 4.8 https://gcc.gnu.org/gcc-4.8/changes.html
        set(DEBUG_FLAGS "${DEBUG_FLAGS} -g3 -Og -fno-inline -ggdb")
    else("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" )
        set(DEBUG_FLAGS "${DEBUG_FLAGS} -g3 -O0 -ggdb")
    endif()

    set(CGAL_DO_NOT_WARN_ABOUT_CMAKE_BUILD_TYPE TRUE)

    set(CHM_BUILD_FLAGS "${CHM_BUILD_FLAGS} ${DEBUG_FLAGS}")
else()
    message(STATUS "Release mode")
    set(RELEASE_FLAGS "-g -O3 -ggdb")

    if(ENABLE_SAFE_CHECKS)
        add_definitions(-DSAFE_CHECKS)
    endif()

    set(CHM_BUILD_FLAGS "${CHM_BUILD_FLAGS} ${RELEASE_FLAGS}")
endif()


#CGAL requires strict rounding
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CHM_BUILD_FLAGS "${CHM_BUILD_FLAGS} -qoverride-limits -fp-model strict -msse4 -finline ")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR
        "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CHM_BUILD_FLAGS "${CHM_BUILD_FLAGS}  -frounding-math")
endif()
########

if(USE_OMP)
    find_package(OpenMP REQUIRED) # This will produce the target for us

    if(OpenMP_FOUND)
        set(CHM_BUILD_FLAGS "${CHM_BUILD_FLAGS} ${OpenMP_CXX_FLAGS}")
        if(OMP_SAFE_EXCEPTION)
            message(STATUS "Enabling safe OMP exception handling")
            add_definitions(-DOMP_SAFE_EXCEPTION)
        endif()
     #should correctly work with cmake 3.12+ https://iscinumpy.gitlab.io/post/omp-on-high-sierra/
    endif()
endif()

#already makes a target
find_package(MeteoIO REQUIRED)

# needs to be title case
find_package(Armadillo REQUIRED)

# needs to be upper-case
# https://cmake.org/cmake/help/latest/module/FindArmadillo.html
create_target(ARMADILLO)

find_package(Eigen3 REQUIRED NO_MODULE)

find_package(BLAS REQUIRED)
create_target(BLAS) # cmake < 3.18


find_package(Sparsehash)
if( Sparsehash_FOUND )
    add_definitions(-DUSE_SPARSEHASH)
endif()

if( ${USE_TCMALLOC})
    find_package(Gperftools)

    #gperftools may not compile on machines w/o nano sleep so we need to optionally disable if it fails to compile
    if(Gperftools_FOUND)
        message(STATUS "Found Tcmalloc, disabling builtin malloc, free")
        set(CHM_BUILD_FLAGS "${CHM_BUILD_FLAGS} -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free")
    endif()
endif()

if( USE_JEMALLOC )
    find_package(Jemalloc 5.3 REQUIRED)
    message(STATUS "Found JeMalloc, disabling builtin malloc, free")
    set(CHM_BUILD_FLAGS "${CHM_BUILD_FLAGS} -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free")
endif()


find_package(Func REQUIRED)

find_package(GDAL 3.0 REQUIRED )
create_target(GDAL)

find_package(PROJ REQUIRED)
create_target(PROJ)

find_package (Threads REQUIRED)

find_package(TBB REQUIRED)

find_package(GSL REQUIRED)
create_target(GSL)

find_package(HDF5 COMPONENTS CXX REQUIRED)
create_target(HDF5)

find_package(NetCDF COMPONENTS C CXX REQUIRED)


#To avoid issues #104, newer VTK version needs to be used.
#8.1+ needs to be used so that proj4 name-collisions are resolved
# http://vtk.1045678.n5.nabble.com/Don-t-build-vtkproj4-td5746073.html
#N.B: The find VTK doesn't allow for compatibility between major versions. That is, if we ask for v7, v8 won't be allowed.
#Therefore we need to ensure this version is up-to-date with what we are building. There are likely issues with VTK6 as per issue #104
# so for now, we will accept finding *any* VTK version, but will flag it and error if we find <=6
find_package(VTK
        9.0
            COMPONENTS
            CommonCore
            CommonDataModel
            CommonTransforms
            FiltersCore
            FiltersGeneral
            FiltersGeometry
            FiltersSources
            IOXML

            REQUIRED)
message(STATUS "Found VTK: ${VTK_DIR}")

#create_target(${VTK})

if(VTK_MAJOR_VERSION EQUAL 8 AND VTK_MINOR_VERSION LESS 1)
    message( FATAL_ERROR "VTK >= 8.1 is required to avoid libproj4 name collision." )
endif()

find_program(GVPR gvpr)
find_program(DOT dot)

if(NOT GVPR OR NOT DOT)
    message(STATUS "Could not find gvpr and/or dot, module output to pdf will not function.")
else()
    message(STATUS "gvpr and dot present, module output to pdf will function.")
endif()

find_package(Trilinos 14.0 REQUIRED)
create_target(Trilinos)

#setup src dirs
include(third_party/CMakeLists.txt)
add_subdirectory(src)

if(BUILD_DOCS)
    add_subdirectory(docs)
endif()