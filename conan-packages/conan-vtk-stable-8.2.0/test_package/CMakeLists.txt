cmake_minimum_required(VERSION 2.8.12)
project(test_package)

set(CMAKE_CXX_STANDARD 11)

include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
conan_basic_setup()

if(APPLE)
    # Workaround for rpath error
    set(CMAKE_SKIP_RPATH 0)
endif(APPLE)

add_executable(${PROJECT_NAME} test_package.cpp )

target_link_libraries(${PROJECT_NAME} ${CONAN_LIBS} ${CMAKE_DL_LIBS})

if(UNIX AND NOT APPLE)
    TARGET_LINK_LIBRARIES(${PROJECT_NAME} m)
endif()
