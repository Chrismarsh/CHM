cmake_minimum_required(VERSION 3.5)
project(PackageTest)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

set(examples
    exact.cpp
    array_convex_hull_2.cpp
    vector_convex_hull_2.cpp
    example.cpp)

add_executable(example ${examples})

target_link_libraries(example ${CONAN_LIBS})

add_library(example-static STATIC ${examples})
target_link_libraries(example-static ${CONAN_LIBS})

add_library(example-shared SHARED ${examples})
target_link_libraries(example-shared ${CONAN_LIBS})