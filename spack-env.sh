#!/bin/bash
export SPACK_PYTHON=/usr/bin/python
. /home/cmarsh/science/spack/spack/share/spack/setup-env.sh
spack env activate ~/science/CHM

spack -E load gcc@14

MPI=/home/cmarsh/science/spack/spack-opt/linux-skylake/openmpi-5.0.7-5yrtznppid7ha5ghb7aodiifiibfj2ny

export CXX=$MPI/bin/mpicxx
export CC=$MPI/bin/mpicc
export FC=$MPI/bin/mpifort