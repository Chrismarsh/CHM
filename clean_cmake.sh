#!/bin/bash

CMAKECACHE=`find . -type f  -name "CMakeCache.txt"`
rm -rf $CMAKECACHE

CMAKEFILES=`find . -type d -name "CMakeFiles"`
rm -rf $CMAKEFILES

MAKEFILES=`find . -type f -name "Makefile"`
rm -rf $MAKEFILES

CMAKEINSTALL=`find . -type f -name "cmake_install.cmake"`
rm -rf $CMAKEINSTALL 

CMAKETEST=`find .  -type f -name "CTestTestfile.cmake"`
rm -rf $CMAKETEST

source /opt/intel/bin/compilervars.sh intel64
cmake -D CMAKE_CXX_COMPILER=icpc -D CMAKE_C_COMPILER=icc .
#cmake .

