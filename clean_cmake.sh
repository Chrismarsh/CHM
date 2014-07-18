#!/bin/bash
 rm -rf \
CMakeCache.txt CMakeFiles/ Makefile cmake_install.cmake \
src/CMakeFiles/ src/Makefile src/cmake_install.cmake\

cmake -D CMAKE_CXX_COMPILER=icpc CMAKE_C_COMPILER=icc .
