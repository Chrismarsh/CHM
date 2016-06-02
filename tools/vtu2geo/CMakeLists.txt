cmake_minimum_required(VERSION 2.8)
 
PROJECT(ReadUnstructuredGrid)
 
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

find_package(GDAL REQUIRED)


find_package(Boost
        1.54.0
        REQUIRED)


include_directories(${GDAL_INCLUDE_DIR} ${Boost_INCLUDE_DIR} )
add_executable(ReadUnstructuredGrid main.cpp )
target_link_libraries(ReadUnstructuredGrid ${VTK_LIBRARIES} ${GDAL_LIBRARY} ${Boost_LIBRARIES})
