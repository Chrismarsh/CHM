Point Mode
===========

Point mode is a special mode in CHM that evaluates the set of modules only at a specific set of points, i.e., 1D mode. Specifically,
only the triangles associated with the output points are run. Four criteria must be met to take advantage of this mode

1. ``point_mode`` is defined in the ``option`` configuration section
2. The module ``point_mode`` is used
3. Only ``parallel::data`` modules are used
4. Have only timeseries outputs at named points

The last constraint limits the module selection to only those modules that operate in a column mode and do not have
any requirement on the surrounding mesh elements. This prohibits the use of modules such as blowing snow and shadowing.



Step 1. Set options flag
--------------------------

In the ``option`` section of the configuration, set ``point_mode``

.. code:: json

   option:
   {
          "point_mode": {
               //empty
            },
   }


Point mode works with either netcdf or ascii file input forcing. If used with ASCII forcing data, a station name may be
optionally given to force that specific station to be used. If this is done then only 1 output point / triangle is
possible to ensure a 1:1 mapping of the two.

.. code:: json

   option:
   {
          "point_mode": {
               "forcing": "my_station"
            },
   }

Please see :ref:`configuration:option` for more details.

Step 2. Add ``point_mode`` to modules
---------------------------------------

The ``point_mode`` module acts as a shim and allows for passing met from the forcing file through to the other modules un touched.
This is useful if you have an insitu met station collection data that you wish to use directly in, e.g., the snow module. Please see the
``point_mode`` module for more information.

Step 3. Only use data parallel
---------------------------------
Only column mode (1D) models are available in ``point_mode``. Currently the only way to confirm this is by checking the source code of a module,
although CHM will inform you at runtime if a ``parallel::domain`` module is chosend

Step 4. Timeseries output
--------------------------
Only timeseries output at named points can be used. E.g.,

.. code:: json

   "output":
   {
      "UpperClearing":
      {
         "longitude": -115.175362,
         "latitude": 50.956547,
         "file": "uc.txt",
         "type": "timeseries"
      },
      "VistaView":
      {
         "longitude": -115.172169,
         "latitude": 50.970944,
         "file": "vv.txt",
         "type": "timeseries"
      }
   },


