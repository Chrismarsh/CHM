mesh gen
==========
Mesh generation is handled by the
`Mesher <https://github.com/Chrismarsh/mesher>`__ program. Mesher
depends heavily upon GDAL to handle the geospatial data and the GDAL
python bindings. Mesher’s input rasters can be in any 1-band raster that
GDAL can open. The triangulation is performed using the Delaunay
triangulation code in CGAL.

Configuration parameters are set in a second .py file and passed as an
argument to ``main.py`` on the command line. For example:

.. code:: bash

   python main.py example_config.py

Experimental support is now enabled for lat/long input files. Due to the
diversity of input data, all input parameters and DEM files to mesher
are set to unified datum, defined by the EPSG number of
``dem_filename``. Further, all files’ nodata value is set to -9999.

The coordinate system of ``dem_filename`` is used for all other files.
However, the user may specify ``EPSG`` in the configuration file, which
will override the coordinate system of ``dem_filename``.

*Although lat long support works, it has implications for normal vector
calculations that are not currently resolved. Therefore all input files
are reprojected to a North American Albers Conic Conformal.*

[STRIKEOUT:If a lat/long (i.e., geographic) dataset is found, a few
things happen:] - [STRIKEOUT:Shortcuts in the meshing step cannot be
taken, so meshing will likely take a bit longer] - [STRIKEOUT:CHM will
scale all lat/long values by 100000 for the vtu output as Paraview seems
to struggle rendering points so close to together.] - [STRIKEOUT:When
triangles are checked during meshing to get the parameter and initial
conditions, the triangle is temporarily projected to an equal-area conic
so-as to determine the area in m^2.]

The extent of ``dem_filename`` is used to define the simulation extent.
Input parameters are constrained to this extent. However, parameters
need not cover the entire extent. Therefore modules *must* check that
paramters are not NaN.

``max_area`` Is a constraint on the maximum size (m^2) of a triangle.

``max_tolerance`` The maximum difference (vertical distance) between the
triangle and the underlying raster

``min_area`` A minimum area (m^2) past which mesher should not refine a
triangle further. A good setting is the square area of a DEM cell. This
does not mean a triangle won’t be smaller than this; rather, if a
triangle is below this threshold it will automatically be accepted as
valid. This will override the tolerance setting. For example, if the
threshold is 3m^, and a 2m^ triangle is checked for validity, it will
automatically be accepted, without checking the tolerance. A triangle
may end up smaller than this threshold due to other splitting that
occurs in order to guarantee triangle quality.

``errormetric`` Assigned an integer value that determines the error
metric to use. ‘mean_tol’ = Mean elevation difference ‘max_tol’ = Max
elevation difference ‘rmse’ = RMSE tolerance

``parameter_files`` is a dictionary that lists additional parameters.
Because a triangle may cover more than one raster cell, the ``method``
variable specifies either ‘mode’ or ‘mean’. This controls how the >1
cells are merged and assigned to a triangle. ‘Mode’ sets the triangle to
be the value that is most common out of all cells.

``initial_conditions`` is a dictionary that lists additional parameters.
Because a triangle may cover more than one raster cell, the ``method``
variable specifies either ‘mode’ or ‘mean’. This controls how the >1
cells are merged and assigned to a triangle. ‘Mode’ sets the triangle to
be the value that is most common out of all cells.

For both initial conditions and parameter files, an optional ‘tolerance’
can be set. If ‘method’ is ‘mode’, then this is a fractional percent of
the dominate landcover to cover the triangle area to not split the
triangle. Otherwise, it is RMSE in the units of the raster’s value.

If the input rasters are of different resolutions, then ``min_area``
should be set to approximately 1/4 the size of the coarsest raster. For
example, if the dem is LiDAR at 1m x 1m, and landcover is present at 30m
x 30m, then ``min_area = 30^2 / 4`` is a good starting point. This is to
prevent the creation of small triangles along the coarse mesh edge.

.. code:: python


   def my_classifier(value):
       if value < .98:
           value = 0
       else:
           value = 1
       return value

       parameter_files = {
           'landcover': { 'file' : 'eosd.tif',
                          'method':'mode',
                           'tolerance':0.8
   },  # mode, mean
           'svf':{'file':'wolf_svf1.tif',
                  'method':'mean',
                  'classifier':my_classifier
                  }
       }
   initial_conditions={
        'sm' : {'file': 'granger_sm_2000.tif', 'method': 'mean'},
       'swe' : {'file': 'granger_swe_2001.tif', 'method':'mean'}
   }

The classifier allows for reclassifying data within mesher so that
data-loss from classifying on a raster is minimized.

Complex basin shapes might result in the creation of many triangles
along the complex edges. Thus ``simplify=True`` and ``simplify_tol`` can
be used to simplify the basin outline. ``simplify_tol`` is the
simplification tolerance specified in meters. Be careful as too high a
tolerance will cause many lines to be outside of the bounds of the
raster.

.. code:: python

   # Configuration file for Mesher
   dem_filename = 'bow_srtm1.tif'
   max_area=1000000
   max_tolerance=50
   min_area = 30**2
   parameter_files={ }
   initial_conditions={ } 
   errormetric = 1 
   simplify     =   False
   simplify_tol =   5   

Mesher creates a directory with the same name as the input dem. This
directory has the reprojected files (``*_projected``), Triangle’s
intermediary files (.node, .elem, .neigh), and the triangulation shape
file (``*_USM.shp``). A ``*.vtu`` file is also created for visualizing
in 3D in Paraview.
