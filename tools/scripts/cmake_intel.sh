#!/bin/bash

source /opt/intel/bin/compilervars.sh intel64
cmake -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc .


