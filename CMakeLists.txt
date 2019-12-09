cmake_minimum_required (VERSION 3.12)

# In CMake 3.12 and above the find_package(<PackageName>) command now searches prefixes specified by the <PackageName>_ROOT CMake variable and the <PackageName>_ROOT environment variable.
#https://cmake.org/cmake/help/latest/policy/CMP0074.html
cmake_policy(SET CMP0074 NEW)
set(CMAKE_FIND_PACKAGE_PREFER_CONFIG TRUE)


project (CHM-project)

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


# Some machines require machine-specific work arounds when building the external libraries. Special casing to a per-machine basis should not be done often
# and really should only be done if 100% required.
# End of December 2018 caused the Westgrid machine Graham to require a special CMAKE_PREFIX_PATH from a system ENV var. This special cases that dramatic change until it's
# fixed in a more seamless manner

#cmake_host_system_information(RESULT hostname QUERY HOSTNAME)
#if(hostname MATCHES "gra-login.*")
#    message(WARNING "Detected machine hostname = Compute Canada Graham")
#    set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};$ENV{CMAKE_SYSTEM_PREFIX_PATH}")
#endif()

# Options. Turn on with 'cmake -Dmyvarname=ON'.
option(USE_MPI "Enable MPI"  OFF )
option(USE_OMP "Enable OpenMP"  ON )
option(BUILD_TESTS "Build all tests."  OFF ) # Makes boolean 'test' available.
option(MATLAB "Enable Matlab linkage"  OFF )
option(STATIC_ANLAYSIS "Enable PVS static anlaysis" OFF)

option(USE_TCMALLOC "Use tcmalloc from gperftools " ON)


set(ENABLE_SAFE_CHECKS FALSE CACHE BOOL "Enable variable map checking. Runtime perf cost. Enable to debug")

# Need to check for MPI before building external libraries
find_package(MPI)
if(MPI_FOUND AND USE_MPI)
  message(STATUS "Found MPI at ${MPI_CXX_INCLUDE_DIRS}")
  add_definitions(-DUSE_MPI)

  file(READ ${CMAKE_BINARY_DIR}/conaninfo.txt txt)
  string(FIND "${txt}" "boost:without_mpi=True" RETCODE)

  # -1 is not found retcode from find
  if(${RETCODE} GREATER "-1")
        message(FATAL_ERROR "MPI was asked for, however boost was not built with mpi. Please set boost:without_mpi=True in conanfile.txt.")
    endif()
else()
  message(WARNING "MPI not enabled. Building for single process.")
endif()


LIST(APPEND CMAKE_PREFIX_PATH "${CMAKE_BINARY_DIR}")
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_BINARY_DIR}")

LIST(APPEND CMAKE_PREFIX_PATH "${CMAKE_SOURCE_DIR}/CMake/")
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake/")

#########
# Helper macro that removes cxx flag flags
macro(remove_cxx_flag flag)
    string(REPLACE "${flag}" "" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endmacro()

if(MACOSX)
    find_program(BREW brew)
    if(BREW)
        execute_process(COMMAND brew --prefix icu4c  OUTPUT_VARIABLE ICU)
        string(REGEX REPLACE "\n$" "" ICU "${ICU}")
        set(ICU "--with-icu=${ICU}" CACHE STRING "icu path for boost build" FORCE)
        message(WARNING "icu = ${ICU}")

        if(ICU STREQUAL "")
            message(WARNING "icu4c not found via brew, likely going to result in boost linking issues" )
        endif()
    endif()
endif()


#on macos, if qt5 is installed with brew, it won't be picked up

#The FindQt5.cmake is the “old way”, and it is still searched for. Qt 5
#distributes CMake config files, such as Qt5Config.cmake that contain paths,
#macros, functions, etc. They are unfortunately a little deeper in the tree, so
#you would have to set Qt5_DIR to “/usr/lib64/cmake/Qt5” in my case. I use this
#approach on Windows, macOS, and Linux. The reduction in Qt5*_DIR variables, and
#the assurance that all modules come from the same Qt5 is worth it.
#https://blog.kitware.com/cmake-finding-qt5-the-right-way/#comment-7191


if(MACOSX)
    find_program(BREW brew)
    if(BREW)
        execute_process(COMMAND brew --prefix qt5  OUTPUT_VARIABLE qt5path)
        message(WARNING "Ammending CMAKE_PREFIX_PATH to include brew qt5 -- ${qt5path}" )
        set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};${qt5path}/lib/cmake/Qt5")
    endif()
endif()

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

find_package(cgal REQUIRED)

if(CGAL_FOUND)
    message(STATUS "Found CGAL ")

    # There is a bug in clang <11 (so far), that has a compiler parse error w/ boost 1.71
    # https://bugs.llvm.org/show_bug.cgi?id=43266
    # https://github.com/Oslandia/SFCGAL/issues/188
    # https://stackoverflow.com/questions/57857572/why-does-boost-log-break-boost-gmp-multiprecision

    # But CGAL uses boost/mp which appears to cause this problem. So if we see clang, disable boost MP usage.
    if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang" )
        ADD_DEFINITIONS(-DCGAL_DO_NOT_USE_BOOST_MP)

        message(WARNING "Disabling CGAL's use of boost mp to avoid clang compiler parse error.")
    endif()
endif()


#remove this flag as it destroys the build on -O1 and spams warnings, #thankscgal
remove_cxx_flag("-Werror=format-security -Wp,-D_FORTIFY_SOURCE=2")

remove_cxx_flag("-O2")
remove_cxx_flag("-O3")


#reset these back

if (CMAKE_BUILD_TYPE MATCHES Debug)
    message(WARNING "Debug mode")
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
    find_package(OpenMP )
    if(OpenMP_FOUND)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
     #should correctly work with cmake 3.12+ https://iscinumpy.gitlab.io/post/omp-on-high-sierra/
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
    add_definitions(-Dvcl_scalar_type=double) #might as well default to double precision if no OpenCL device found.
    unset(OPENCL_LIBRARIES CACHE) ##zero these out so we can use them later without trouble
    unset(OPENCL_INCLUDE_DIRS CACHE)
endif()

find_package(meteoio REQUIRED)

find_package(armadillo REQUIRED)

find_package(eigen3 REQUIRED)

find_package(sparsehash)

#sparsehash doesn't compile with gcc 4.6.x - 4.8.x because of this
# https://bugs.launchpad.net/ubuntu/+source/gcc-4.6/+bug/1269803
# so we need to optionally disable it
if( SPARSEHASH_FOUND )
    add_definitions(-DUSE_SPARSEHASH)
else()
    set(SPARSEHASH_INCLUDE_DIRS "")
endif()

find_package(gperftools)
if(gperftools_FOUND)  #gperftools may not compile on machines w/o nanosleep so we need to optionall disable if it fails to compile
    message(STATUS "Found Tcmalloc, disabling builtin malloc, free")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-builtin-malloc -fno-builtin-calloc -fno-builtin-realloc -fno-builtin-free")
else()
    set(Tcmalloc_INCLUDE_DIR "")
    set(Tcmalloc_LIBRARY "")
endif()


find_package(viennacl  REQUIRED)


message(STATUS "Found ViennaCL at ${VIENNACL_INCLUDE_DIR}")

if(BUILD_FUNC)
    find_package(func REQUIRED)
endif()
message(STATUS "Found func at ${FUNC_INCLUDE_DIR}")
message(STATUS "Found func libraries at ${FUNC_LIBRARY}")
message(STATUS "Found func libraries at ${FUNC_IMPLS_LIBRARY}")


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
    REQUIRED)
endif()
message(STATUS "Found boost incl: ${Boost_INCLUDE_DIR}")
message(STATUS "Boost libs to link: ${Boost_LIBRARIES}")
#need this for calling the logger from multiple modules
ADD_DEFINITIONS(-DBOOST_LOG_DYN_LINK)

find_package(Iconv REQUIRED) #needed for gdal
find_package(gdal 2.4 REQUIRED )

find_package (Threads REQUIRED)

find_package(tbb REQUIRED tbb)

find_package(gsl REQUIRED)

set(NETCDF_CXX "YES")
find_package(netcdf-cxx  REQUIRED)
message(STATUS "netcdf lib ${NETCDF_LIBRARIES}")


#To avoid issues #104, newer VTK version needs to be used.
#8.1+ needs to be used so that proj4 name-collisions are resolved
# http://vtk.1045678.n5.nabble.com/Don-t-build-vtkproj4-td5746073.html
#N.B: The find VTK doesn't allow for compatibility between major versions. That is, if we ask for v7, v8 won't be allowed.
#Therefore we need to ensure this version is up-to-date with what we are building. There are likely issues with VTK6 as per issue #104
# so for now, we will accept finding *any* VTK version, but will flag it and error if we find <=6
find_package(vtk REQUIRED)

if(VTK_MAJOR_VERSION LESS 7)
    #https://github.com/Chrismarsh/CHM/issues/104
    message( FATAL_ERROR "VTK > 6 is required, please see issue #104. You can build a compatible VTK using -DBUILD_VTK=ON." )
endif()

if(VTK_MAJOR_VERSION EQUAL 8 AND VTK_MINOR_VERSION LESS 1)
    message( FATAL_ERROR "VTK >= 8.1 is required to avoid libproj4 name collision. You can build a compatible VTK using -DBUILD_VTK=ON." )
endif()

#endif()

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

if (BUILD_TESTS)
    message(STATUS "Tests enabled. Run with make check")
    enable_testing()
endif()


message(STATUS "${CMAKE_CXX_FLAGS}")
