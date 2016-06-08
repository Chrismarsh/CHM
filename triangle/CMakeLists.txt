cmake_minimum_required (VERSION 2.8)

project (triangle)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 ")

include_directories(${GDAL_INCLUDE_DIR})
find_package(gdal)
add_executable(
	triangle
	triangle.c
 )

target_link_libraries(
	triangle
	m
		${GDAL_LIBRARY}
)