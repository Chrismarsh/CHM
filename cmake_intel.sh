#!/bin/bash

source /opt/intel/bin/compilervars.sh intel64
cmake -D CMAKE_CXX_COMPILER=icpc -D CMAKE_C_COMPILER=icc .


