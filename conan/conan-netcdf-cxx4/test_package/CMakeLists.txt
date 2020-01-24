cmake_minimum_required(VERSION 2.8.12)
project(PackageTest C CXX)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup(KEEP_RPATHS)

add_executable(example simple_xy_wr.cpp)
target_link_libraries(example ${CONAN_LIBS} )
