cmake_minimum_required(VERSION 3.4)
project(mesher)

#as per http://cgal-discuss.949826.n4.nabble.com/CMake-and-flags-td949906.html
#don't override internal settings
set( CGAL_DONT_OVERRIDE_CMAKE_FLAGS TRUE )
find_package(CGAL REQUIRED)
if(CGAL_FOUND)
    message(STATUS "Found CGAL ")
endif()
include(${CGAL_USE_FILE}) #as per https://www.cgal.org/releases.html release 4.2

#ignore these two under Clion as CGAL will complain
if(CMAKE_BUILD_TYPE MATCHES RelWithDebInfo OR
        CMAKE_BUILD_TYPE MATCHES MinSizeRel OR
        NOT CMAKE_BUILD_TYPE)

    set(CMAKE_BUILD_TYPE "Release")
endif()

#reset these back
#if (CMAKE_BUILD_TYPE MATCHES Debug)
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 ")
#else()
#    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 ")
#endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -g -O3")

find_package(GDAL REQUIRED)
find_package(Boost
        COMPONENTS
        program_options
        filesystem
        REQUIRED)

include_directories(${GDAL_INCLUDE_DIR}
        ${Boost_INCLUDE_DIR}
        )


set(SOURCE_FILES mesher.cpp triangle.cpp triangle.h raster.cpp raster.h)
add_executable(mesher ${SOURCE_FILES})
target_link_libraries(
        mesher
        m
        ${GDAL_LIBRARY}
        ${Boost_LIBRARIES}
)
