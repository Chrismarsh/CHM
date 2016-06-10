cmake_minimum_required (VERSION 2.8)

project (triangle)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -g ")

find_package(gdal 2.0)

include_directories(${GDAL_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(
	triangle
		triangle.c
		triunsuitable.c)

target_link_libraries(
	triangle
	m
		${GDAL_LIBRARY}
)