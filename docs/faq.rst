Running CHM
===========

Westgrid’s Graham
-----------------

If the google memory allocator is used, adding it to the LD_LIBRARY_PATH
is needed, otherwise an incorrect allocator is use

``LD_LIBRARY_PATH=~/build/lib/VTK/lib:~/build/lib/gperftools/lib:~/build/lib/gsl/lib:~/build/lib/boost/lib:~/build/lib/netcdf-c/lib:~/build/lib/ViennaCL/lib:~/build/lib/sparsehash/lib:~/build/lib/gdal/lib:~/build/lib/hdf5/lib:~/build/lib/proj4/lib ./CHM``
