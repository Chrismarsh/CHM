#!/usr/bin/env bash

# script to export dep versions from all (submodule) dependencies

CHM_TYPE=CHM/stable

conan export conan/conan-armadillo 10.2.0@$CHM_TYPE
conan export conan/conan-boost 1.75.0@$CHM_TYPE
conan export conan/conan-cgal 5.2@$CHM_TYPE
conan export conan/conan-eigen3 3.3.9@$CHM_TYPE
conan export conan/conan-func 0.1@$CHM_TYPE
conan export conan/conan-gdal 3.2.2@$CHM_TYPE
conan export conan/conan-gmp 6.2.1@$CHM_TYPE
conan export conan/conan-gperftools 2.8.1@$CHM_TYPE
conan export conan/conan-gsl 2.6@$CHM_TYPE
conan export conan/conan-hdf5 1.12.0@$CHM_TYPE
conan export conan/conan-meteoio 2.8.0@$CHM_TYPE
conan export conan/conan-mpfr 4.1.0@$CHM_TYPE
conan export conan/conan-netcdf-c 4.7.4@$CHM_TYPE
conan export conan/conan-netcdf-cxx4 4.3.1@$CHM_TYPE
conan export conan/conan-proj 7.2.1@$CHM_TYPE
conan export conan/conan-sparsehash 2.0.3@$CHM_TYPE
conan export conan/conan-tbb 2020.3@$CHM_TYPE
conan export conan/conan-trilinos 12.18.1@$CHM_TYPE
conan export conan/conan-vtk 9.0.1@$CHM_TYPE
