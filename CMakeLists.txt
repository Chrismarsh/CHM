cmake_minimum_required (VERSION 3.16)
project (CHM-project CXX Fortran)

# In CMake 3.12 and above the find_package(<PackageName>) command now searches prefixes specified by the <PackageName>_ROOT
# CMake variable and the <PackageName>_ROOT environment variable.
#https://cmake.org/cmake/help/latest/policy/CMP0074.html
cmake_policy(SET CMP0074 NEW)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

#http://stackoverflow.com/questions/31561309/cmake-warnings-under-os-x-macosx-rpath-is-not-specified-for-the-following-targe
#https://cmake.org/cmake/help/v3.0/policy/CMP0042.html
set(CMAKE_MACOSX_RPATH 1)
if (${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    set(MACOSX TRUE)
    set(CMAKE_FIND_FRAMEWORK NEVER) #Can cause frameworks to take priority over superbuild libraries
endif()


# Options. Turn on with 'cmake -Dmyvarname=ON'.
option(USE_MPI "Enable MPI"  OFF )
option(USE_OMP "Enable OpenMP"  ON )
option(BUILD_TESTS "Build all tests."  OFF ) # Makes boolean 'test' available.
option(MATLAB "Enable Matlab linkage"  OFF )
option(STATIC_ANLAYSIS "Enable PVS static anlaysis" OFF)
option(USE_TCMALLOC "Use tcmalloc from gperftools " ON)
option(USE_JEMALLOC "Use jemalloc" OFF)
option(BUILD_DOCS "Builds documentation" OFF)
option(BUILD_WITH_CONAN "Uses Conan for dependencies" ON)

set(ENABLE_SAFE_CHECKS FALSE CACHE BOOL "Enable variable map checking. Runtime perf cost. Enable to debug")

if(BUILD_WITH_CONAN AND USE_JEMALLOC)
    message(WARNING "With a Conan build USE_JEMALLOC will use a system jemalloc" )
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

#LIST(APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/CMake/")
#LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake/")

if(NOT BUILD_WITH_CONAN)
    LIST(APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/CMake-without-conan/")
    LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake-without-conan/")
endif()

#########
# Helper macro that removes cxx flag flags
macro(remove_cxx_flag flag)
    string(REPLACE "${flag}" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endmacro()

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

if(NOT BUILD_WITH_CONAN)
    set(CMAKE_BUILD_TYPE NoConan)
endif()

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
include(${CMAKE_CURRENT_SOURCE_DIR}/cmake_find_names.cmake)

# Need to check for MPI before building external libraries
find_package(MPI)
if(MPI_FOUND AND USE_MPI)
    message(STATUS "Found MPI at ${MPI_CXX_INCLUDE_DIRS}")
    add_definitions(-DUSE_MPI)

    if(BUILD_WITH_CONAN)
        # ensure boost MPI support was built
        if(MPI_FOUND AND USE_MPI)
            file(READ ${CMAKE_BINARY_DIR}/FindBoost.cmake txt)
            string(FIND "${txt}" "boost_mpi" RETCODE)

            if(${RETCODE} EQUAL -1) # -1 is not found retcode from find
                message(FATAL_ERROR "MPI was asked for, however boost was not built with mpi. Please build boost (--build boost) locally with a working MPI. See wiki for more information.")
            endif()
        endif()
    endif()
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
ADD_DEFINITIONS(-DBOOST_LOG_DYN_LINK)

find_package(${CGAL} 5.0 REQUIRED)

# There is a bug in clang <11 (so far), that has a compiler parse error w/ boost 1.71
# https://bugs.llvm.org/show_bug.cgi?id=43266
# https://github.com/Oslandia/SFCGAL/issues/188
# https://stackoverflow.com/questions/57857572/why-does-boost-log-break-boost-gmp-multiprecision
# But CGAL uses boost/mp which appears to cause this problem. So if we see clang, disable boost MP usage.
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" )
    ADD_DEFINITIONS(-DCGAL_DO_NOT_USE_BOOST_MP)

    message(WARNING "Disabling CGAL's use of boost mp to avoid clang compiler parse error.")
endif()

#reset these back
if (CMAKE_BUILD_TYPE MATCHES Debug)
    message(WARNING "Debug mode")
    if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND CMAKE_CXX_COMPILER_VERSION VERSION_GREATER "4.7")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3 -O1 -fno-inline") #new to 4.8 https://gcc.gnu.org/gcc-4.8/changes.html
    else()
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g3 -O1 -fno-inline -ggdb")
    endif()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSAFE_CHECKS")
else()

    if(ENABLE_SAFE_CHECKS)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DSAFE_CHECKS")
    endif()
    set(CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -g -O3 ")
endif()

#Setup CXX flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-variable -Wno-unknown-pragmas")


#GAL requires strict rounding
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qoverride-limits -fp-model strict -msse4 -finline ") #-frounding-math needed for cgal
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU") # OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -frounding-math")
endif()
########

if(USE_OMP)
    find_package(OpenMP ) # This will produce the target for us

    if(OpenMP_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
     #should correctly work with cmake 3.12+ https://iscinumpy.gitlab.io/post/omp-on-high-sierra/
    endif()
endif()


if(MATLAB)
    #for matlab to run you need csh installed
    find_package(Matlab REQUIRED)
endif()


#already makes a target
find_package(${METEOIO} REQUIRED)

find_package(${ARMADILLO} REQUIRED)
create_target(ARMADILLO) #needs to be upper-case

if(BUILD_WITH_CONAN)
    find_package(${EIGEN3} REQUIRED)
else()
    find_package(${EIGEN3} REQUIRED NO_MODULE)
endif()

# The trilinos library may have been built against a custom location blas/openblas
# If so, conaninfo.txt has information on the custom blas location if given with the -o trilinos:blas_root=<root> option
# If conaninfo.txt has a custom blas location, extract it and pass it as a hint to find_package(blas) which will detect any blas.
# Then use the find_package link target in CHM to fully avoid any problems of linking against the wrong blas or not being able to find it.
if(BUILD_WITH_CONAN)
    if(EXISTS "${CMAKE_BINARY_DIR}/conaninfo.txt")
        file(READ ${CMAKE_BINARY_DIR}/conaninfo.txt conaninfotxt)
        string(REGEX MATCH "trilinos:blas_root=([^ ]+)lib" blasloc "${conaninfotxt}") # because cmake regex is only greedy
        if(NOT "${blasloc}" STREQUAL "")
            string(REPLACE "\n" "" blasloc ${CMAKE_MATCH_1}) # strip out a new line we've picked up
            message(WARNING "Using the conaninfo.txt trilinos:blas_root=${blasloc}")
            set(BLAS_ROOT ${blasloc})
        endif()
    elseif()
        set(blasloc "")
    endif()
endif()

if(MACOSX AND "${blasloc}" STREQUAL "")

    execute_process(
            COMMAND brew --prefix openblas
            RESULT_VARIABLE BREW_BLAS
            OUTPUT_VARIABLE BREW_BLAS_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
    )

    if(BREW_BLAS EQUAL 0)
        message(WARNING "Found homebrew OpenBLAS ${BREW_BLAS_PREFIX}" )
        set(BLAS_ROOT ${BREW_BLAS_PREFIX})
        set(BLA_VENDOR Generic)
    endif()

endif()

find_package(${BLAS}
            REQUIRED
            )

# cmake < 3.18 doesn't produce a BLAS target
create_target(${BLAS})


find_package(${SPARSEHASH})
if( ${SPARSEHASH}_FOUND )
    add_definitions(-DUSE_SPARSEHASH)
endif()

if( ${USE_TCMALLOC})
    find_package(${GPERFTOOLS})

    #gperftools may not compile on machines w/o nano sleep so we need to optionally disable if it fails to compile
    if(${GPERFTOOLS}_FOUND)
        message(STATUS "Found Tcmalloc, disabling builtin malloc, free")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free")
    endif()
endif()

if( USE_JEMALLOC )
    find_package(Jemalloc REQUIRED)
    message(STATUS "Found JeMalloc, disabling builtin malloc, free")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free")
endif()


find_package(${FUNC} REQUIRED)

find_package(${GDAL} 3.0 REQUIRED )
create_target(${GDAL})

find_package(${PROJ} REQUIRED)
create_target(${PROJ})

find_package (Threads REQUIRED)

find_package(${TBB} REQUIRED)

find_package(${GSL} REQUIRED)
create_target(${GSL})

find_package(${HDF5} COMPONENTS CXX REQUIRED)
create_target(${HDF5})

find_package(${NETCDFCXX} COMPONENTS C CXX REQUIRED)


#To avoid issues #104, newer VTK version needs to be used.
#8.1+ needs to be used so that proj4 name-collisions are resolved
# http://vtk.1045678.n5.nabble.com/Don-t-build-vtkproj4-td5746073.html
#N.B: The find VTK doesn't allow for compatibility between major versions. That is, if we ask for v7, v8 won't be allowed.
#Therefore we need to ensure this version is up-to-date with what we are building. There are likely issues with VTK6 as per issue #104
# so for now, we will accept finding *any* VTK version, but will flag it and error if we find <=6
find_package(${VTK}
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
if(NOT ${BUILD_WITH_CONAN})
    message(STATUS "Found VTK: ${VTK_INCLUDE_DIRS}")
endif()

#create_target(${VTK})

if(VTK_MAJOR_VERSION EQUAL 8 AND VTK_MINOR_VERSION LESS 1)
    message( FATAL_ERROR "VTK >= 8.1 is required to avoid libproj4 name collision. You can build a compatible VTK using -DBUILD_VTK=ON." )
endif()

find_program(GVPR gvpr)
find_program(DOT dot)

if(NOT GVPR OR NOT DOT)
    message(STATUS "Could not find gvpr and/or dot, module output to pdf will not function.")
else()
    message(STATUS "gvpr and dot present, module output to pdf will function.")
endif()

find_package(${TRILINOS} REQUIRED)
create_target(${TRILINOS})

#setup src dirs
include(third_party/CMakeLists.txt)
add_subdirectory(src)

if(BUILD_DOCS)
    add_subdirectory(docs)
endif()