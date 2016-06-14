cmake_minimum_required (VERSION 2.8)

project (triangle)

set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O0 -g -fPIC") #fPIC is needded for building under CLion on linux...no real idea why.

find_package(GDAL REQUIRED)

include_directories(${GDAL_INCLUDE_DIR} ${CMAKE_CURRENT_SOURCE_DIR})

add_executable(
	triangle
		triangle.c
		triunsuitable.c
		tri.c
		masked_array.c error_metrics.c
)

target_link_libraries(
	triangle
	m
		${GDAL_LIBRARY}
)