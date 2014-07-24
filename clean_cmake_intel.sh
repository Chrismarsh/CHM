#!/bin/bash
 rm -rf \
CMakeCache.txt CMakeFiles/ Makefile cmake_install.cmake \
src/CMakeFiles/ src/Makefile src/cmake_install.cmake\

source /opt/intel/bin/compilervars.sh intel64
cmake -D CMAKE_CXX_COMPILER=icpc -D CMAKE_C_COMPILER=icc .


