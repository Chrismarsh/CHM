Tools
============


vtu2geo
--------
The conversion of the vtu format to arbitrary GIS formats is provided by
``vtu2geo`` located in ``tools/vtu2geo/main.py``. Requires vtk and gdal
+ python bindings for each.

This tool produces an internal shp file that corresponds to the
triangulation and uses GDAL to rasterize this to an output geotiff.

The vtu files contain multiple variables. Therefore, each output geotiff
is a 1-band file corresponding to the selected output. The variables of
interest as set in the ``variables`` list.

For vtu variables that are parameters (and therefore constant with
time), only 1 output file is needed. These are defined in the
``parameters`` list. Only 1 geotiff will be produced from these.

``input_path`` should be either a single .vtu file or a the .pvd file.
If a pvd file is given it will produce a tiff for each vtu.

``output_path`` is the output path. If output into the current folder is
wanted, use ``'.'``.

``pixel_size`` is the size of the raster cells in m^2.
``var_resample_method`` and ``param_resample_method`` determine what
resampling method to use when calculating on a clipped raster.

.. code:: python

   # Configuration file for vtu2geo tool


   # Input path to where output vtu files are located
   input_path = 'meshes/SC1506124800.vtu'
   output_path = 'meshes/'


   # Output variables
   variables = ['snowdepthavg']  #set to None to dump all variables
   var_resample_method    = {'snowdepthavg':'average'} # Methods to use when calculating clipped raster

   # Output parameters
   parameters = ['Elevation'] # paramters are one offs we want to extract from the vtu files
   param_resample_method    = {'Elevation':'average'} # Methods to use when calculating clipped raster


   # Output pixel size that the mesh is interpolated to (?)
   pixel_size = 30 # (m)
