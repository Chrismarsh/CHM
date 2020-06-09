Configuration
===============


.. note::

   Support for full geographic meshes exists, however there is currently an issue with computing surface normals. Thus projected coordinate systems should only be used.


.. note::
   Regardless of input coordinate system **all** input points are specified in latitude and longitude in WSG84.


Schema
------------

The config file is a JSON file. However, it does support C-style comments: ``//`` and ``/** **/`` are both valid. 

There are a few required sections: ``modules``, ``meshes``, ``forcing``.

The general layout of a CHM config JSON file is


.. code:: json 

   {
      "option":
      {
         "option_a": true,
         "option_b": 1234,
         ...
      },
      "modules":
      [
         "module1",
         "module2",
         ...
      ],
      "config":
      {
         "module1":
         {
            ...
         },

         ...
      }
      "meshes":
      {
         ...
      },
      "forcing":
      {
        ...
      },
      "output":
      {
         ...
      }
   }

For every section, if a top-level key:value pair is found and the value contains ".json", that file is loaded and inserted into this option. The key-value is not used, and may be anything. Key names are enclosed in quotes (" "). Although it tends to make more sense to arrange the keys in the shown order, the order does not matter (anywhere) and will be read correctly.


.. warning::
   The ``modules`` key is an array and requires the use of [ ]


.. warning::

   Do not prefix a number with zero (0). This is the
   octal prefix and it will cause the JSON parser to choke on lines that
   otherwise look fine.

.. note::
   A user can specify a number as ``"5"`` or ``5``. Internally to CHM it will be converted to a numeric type. Thus, both are fine, however a non-string should be preferred. This is similar for ``"true"`` and ``true``. 

.. note::
   Boolean types are case sensitive.

Sections
---------

option
********

This section contains options for CHM and the simulation in general. 

.. code:: json

   {
      "option":
      {
           "station_N_nearest": 1,
           "interpolant": "nearest",
           "per_triangle_timeseries": false,
           "ui": false,
           "debug_level": "debug",
           "prj_name": "SnowCast",
           "enddate": "20180501T050000"

      }
   }

.. confval:: station_search_radius

   :type: double
   :default: None

   
   The search radius (meters) surrounding any given triangle within which to search for a station. This is used to ensure only "close" stations are used. Cannot be used when ``station_N_nearest`` is set. Based off the center of the triangle. 


.. confval:: station_N_nearest
   
   :type: int
   :default: 5

   Use the nearest N stations to include for the interpolation at a triangle. Based off the center of the triangle. 

   Both ``station_search_radius`` and ``station_N_nearest`` cannot be
   simultaneously specified. If neither is specific, then ``station_N_nearest:5`` is used as default. If the :confval:`interpolant` mode is ``nearest``, then this is automatically set to 1.


.. confval:: interpolant

   :type: string
   :default: "spline"

   Chooses either thin plate spline with tension (spline) or inverse
   distance weighting (idw). Nearest selects the closest
   station and only uses that with no interpolation. 

   .. code:: json 

      "interpolant" : "idw"
      "interpolant" : "spline"
      "interpolant" : "nearest"

.. confval::  point_mode
   
   :type: ``{ }``
   :required: No

   Point mode selects that the model should be run in point mode, versus
   distributed mode. 

   There are two keys that need to be specified:

   - ``output`` (string)
   - ``forcing`` (string)

   ``output`` needs to correspond to a specific output point of type timeseries as defined in the output section.

   ``forcing`` needs to correspond to a specific input point as defined in the forcing section

   If other points are specified, these points will be ignored.

   Usage of this key also requires adding ``point_mode`` to the module list. Lastly, no
   modules which are defined ``parallel:domain`` may be used when point_mode is enabled.

.. code:: json 

       "point_mode":
       {
         "output":"UpperClearing",
         "forcing":"UpperClearing"
       },

.. confval:: notification_script

   :type: string
   :default: None

   Path to a script to call upon model execution. This is useful
   for sending a notification to a computer or phone upon the completion of
   a long model run.

.. code:: json 

       "notification_script":"./finished.sh"

And example of what ``finished.sh`` might do is below, which triggers a
notification to Pushbullet thus showing up on all computers and phones
that the account is active on:

.. code::  

   #!/bin/bash

   curl -s -u <token here>: https://api.pushbullet.com/v2/pushes -d type=note -d title="Finished model run" >/dev/null



.. confval:: debug_level

   :type: string
   :default: "Debug"

   This controls the verbosity of the output. Options are: 

   - verbose [ all messages ] 
   - debug [ most messages useful for debugging ] 
   - warning [only warnings] 
   - error [ only errors which terminate model execution ]

   Currently most useful internal messages are debug level.

.. code:: json 

       "debug_level":"debug"


.. confval:: startdate
   
   :type: string
   :default: None


Allows for a different start start time that that specified by the input timeseries.
In the same ISO format as the forcing data: ``YYYYMDTHMS``.

.. code:: json 

   "startdate":"20010501T000000"

.. confval:: enddate
   
   :type: string
   :default: None

Allows for a different start start time that that specified by the input timeseries.
In the same ISO format as the forcing data: ``YYYYMDTHMS``.

.. code:: json 

   "enddate":"20010502T000000"

modules
********

Modules to run. These are a comma separated list of keys. This is a required section.

A few notes:

- order as defined in this list has no bearing on the order modules execute
- may be commented out to remove them from execution
- names are case sensitive
- ``point_mode`` module is required to enable point mode, in addition to being enabled in ``option.point_mode``.

.. note::
   Modules are in a list (``[ ]``) 

.. code:: json 

     "modules":
     [
        "Liston_wind",
        "Burridge_iswr",
        "slope_iswr",
        "Liston_monthly_llra_ta",
        "kunkel_rh",
        "Thornton_p",
        "Walcek_cloud",
        "Sicart_ilwr",
        "Harder_precip_phase",
        "snobal",
        "Gray_inf",
        "Richard_albedo"

     ]

remove_depency
***************

   Under some edge cases, a cyclic dependency is created when a module B
   depends on module A’s output, and module A depends on module B’s output. There is no way to
   automatically resolve this. It requires the modeller to manually break
   the cycle and force one module to run ahead of another (essentially
   time-lagged).

   An example of this occurring is that the albedo model requires knowledge
   of SWE, provided by the snowmodel. However, the snowmodel requires
   albedo to run. Therefore, the modeller may define that the albedo
   routine is run first, then the snowpack model.

   Specifically: if module A depends on B (A->B), then to remove the decency
   of B from A, specify it as ``"A":"B"``

   This can be thought of as ``A`` needs to come before ``B``. If the specified modules are not added to the modules list, they are ignored.

   .. code:: json 

        "remove_depency":
        {
          "Richard_albedo":"snobal"
        }

   

config
*******

Each module, upon creation is provided a configuration instance. These configuration data are set by creating a
key that exactly matches the module name. If a section is added, but that module isn't specified, the section is ignored.

.. confval:: module_name

   :type: ``{ }``



For example:

.. code:: json 

   "config":
   {

      "slope_iswr":
          {
            "no_slope":true
          }
   },


If the configuration is sufficiently large or cumbersome, it may be best
to have it in a separate file. This can be specified as

.. code:: json 

   //consider this in CHM.json
   "config":
   {
       "simple_canopy":"canopy.json"   
   }

   ​
And ``canopy.json`` is 

.. code:: json

   "canopy": 
   {
     "LAI":3 
   }
   


Note that the sub-keys for a module's configuration are entirely dependent upon the module. Please see the module's help for specific options.

meshes
*******

This section defines the mesh and optional the parameter files to use. It is a require section.
This section has two keys:

.. confval:: mesh

   :type: string


   File path  to the ``.mesh`` file produced by mesher.

.. confval:: parameters

   :type: ``{ }``

   Optionally, A set of key:value pairs to other ``.param`` files that contain extra parameters to be used.
   These are in the format ``{ "file":"<path>"" }``


.. code:: json

   "meshes":
   {
    "mesh":"meshes/granger30.mesh",
    "parameters":
    {
      "file":"meshes/granger30.param",
      "file":"meshes/granger30_surface.param"
    }
   }





parameter_mapping
******************

The parameters may be classified values for use in a look-up table. For example, the landcover may be a numeric class value and values such as LAI need to be obtained from a lookup table. These parameters may be either specified directly in the file or located in another file:

.. code:: json 

     "parameter_mapping":
     {
       "soil":"parameters/wolf_soil_param.json"
     }

or as a key:value pair. In all cases, the parameter name is how it will
be referenced in the module that is looking for it. Please see the module's documentation for what the expected format is.

.. code:: json

      {
         "landcover":
         {
            "20":
            {
              "desc":"lake",
              "is_waterbody":true
            },
            "31":
            {
              "desc":"snow ice"
            }
         }
      }

output
*********

Output may be either to an ascii-timeseries for a specific triangle on the mesh
or it may be the entirety of the mesh. The two output types are set by:

   - a key named ``"mesh":{ ... }`` will enable the entire mesh output
   - all other keys (``"some_name":{...}```) are assumed to be the names of output timeseries

Both mesh and timeseries can be used together.


.. confval:: output_dir

   :type: string
   :default: "output"

   The output directory name.



timeseries output
~~~~~~~~~~~~~~~~~~

The name of the ``timeseries`` key is used to uniquely identify this output: ``"output_name"{ ... }``. 

If using ``point_mode``, this name corresponds to the ``output`` key. If a lot of stations are to be
output, consider keeping them in a separate file and inserting using the top-level ".json" behaviour.

.. confval:: longitude

   :type: float

   WGS84 longitude of output point. The triangle that contains this point is then selected for output. An error is raised if no triangle contains the point.

.. confval:: latitude

   :type: float

   WGS84 latitude of output point. The triangle that contains this point is then selected for output. An error is raised if no triangle contains the point.

.. confval:: file

   The output file name. The output is in csv format and each column is a variable.


.. code:: json 

     "output":
     {
        "more_stations":"mystations.json",
        "UpperClearing": 
        {
            "longitude": "-115.175362",
            "latitude": "50.956547",
            "file": "uc.txt"
        }
    }

where ``mystations.json`` would look like

.. code:: json

   {
        "some_station": 
        {
            "longitude": "-115.175362",
            "latitude": "50.956547",
            "file": "somestation.txt"
        }
   }

mesh
~~~~

The entire mesh may be written to Paraview’s vtu format for
visualization in Paraview and for analysis. This is denoted by a ``"mesh":{ ... }`` key.

For more details, please see the `output`__ section.

.. confval:: base_name
   
   :type: string
   
   The base file name to be used. 

.. confval:: variables

   :type: ``[ "variable", ... ]``

   The default behaviour to is write every variable at each timestep. This may produce an undesirable amount of output. This takes a list of variables to output.

.. code:: json

   "variables": [
                "t",
                "U_2m_above_srf",
                "swe",
                "iswr"
            ],

.. confval frequency::

   :type: int
   :default: 1

   Frequency can be set to write ever *N* timesteps. 

.. confval write_parameters::

   :type: boolean
   :default: true

   Disables/enables writing parameters to the output.

Example:

.. code:: json

   "output":
   {
    "mesh": {
            "base_name": "SC",
            "variables": [
                "t",
                "U_2m_above_srf",
                "swe",
                "iswr"
            ],
            "frequency": "24",
            "write_parameters": false
        }
   }







Forcing
--------

Input forcing can be either a ASCII timeseries or a NetCDF.

Input forcing stations do not need to be located within the simulation
domain. Therefore they can act as ‘virtual stations’ so-as to use
reanalysis data, or met stations located outside of the basin.

An example of this is shown below, where each black point is a virtual station, representing the center for a NetCDF grid cell from a NWP product.

.. image:: images/netcdf.png


.. confval:: UTC_offset

   :type: int
   :default: 0

 If the input timeseries it not UTC, then this is the correction to account for UTC offset (all solar radiation calculations are in UTC).
 This is Positive west!

.. confval:: use_netcdf

   :type: boolean
   :default: false

   Specify if a NetCDF (.nc) file will be used. Cannot be used along with ASCII inputs!



.. note::

   ASCII and NetCDF inputs cannot be mixed. It is one or the other.

ASCII timeseries
*****************

Forcing data are defined as a delineated format with ISO datetime. This is given as ``"station_name":{ ... }``. If using ``point_mode``, then the value ``station_name`` must exactly match the ``input`` used for ``option.point_mode``.

.. code:: json

 Date                 Rh   Tair  Precip
   20080220T000000    50      -12      2
   20080221T000015    40      -10      0

where the input variable names correspond to the variable names the selected modules expect. Please refer to those modules' documentation.

Some restrictions:
        - No more than 2147483647 steps. At 1s intervals, this equates to roughly 68 years.
        - Consistent units. You mustn't have mm on one line, then meters on the next, for the same observation
        - Has to be on a constant time step. The first interval is taken as the interval for the rest of the file
        - Missing values are not currently allowed - that is, each row must be complete with n entries where n is number of variables.
        - However, a missing value value (i.e., -9999) can be used to represent missing data
        - Whitespace, tab or comma delimited. Allows for mixed usage. ex 1234, 4543 890 is legal
        - Values must be numeric

Integer styles:
      -   +1234
      -   -1234
      -   1234567890

Floating point:
      -   12.34
      -   12.
      -   .34
      -   12.345
      -   1234.45
      -   +12.34
      -   -12.34
      -   +1234.567e-89
      -   -1234.567e89

Time:
        - Must be in one column in the following ISO 8601 date time form:
        - YYYYMMDDThhmmss   e.g., 20080131T235959


.. confval:: file

   A relative or absolute path to an input forcing file


.. confval:: latitude

   :type: double

   Latitude of the input station, WGS84. Positive North. Not "N" or "S" suffix

.. confval:: longitude

   :type: double

   Longitude of the input station, WGS84. Positive East. Not "N" or "S" suffix

.. confval:: elevation

   :type: double

   Elevation is given in meters. It does *not* need to be equal to the elevation of the triangle upon which it lies if the station is located in the simulation domain.
   This value is used in the lapse rate equations to interpolate the data.

If required, forcing station definitions can be located in an external
file. For the external file, the name of the key doesn’t matter. The
external file should contain the stations in the format as per above. It
does *not* require an addition ``"forcing":`` section definition.

.. code:: json 

   "forcing":
     {
       "some_station":
        {
          //definition
        },
       "reanalysis_extract_1": "external_file_1.json",
      "reanalysis_extract_2": "external_file_2.json",
   }

where ``external_file_*.json`` looks like

.. code:: json 

   {
    "station1":
      {
       //details here
      },
    "station2":
      {
       //details here
      }
   }


Filters
~~~~~~~~

Filters perform an operation on the data prior to being passed to a module. They allow for things such as wind-undercatch corrections to be done on the fly. 

If a filter is defined, it must be defined on the forcing file and operate upon a variable that exists in the forcing data. They are given as:

``"filter_name": { ... }```. The configuration values are filter-specific; please see the filter documentation for what is required. Multiple filters may be specified.

.. code-block:: json

  "buckbrush": 
  {
    "file": "bb_m_2000-2008",
    "latitude": 60.52163,
    "longitude": -135.197151,
    "elevation": 1305,
    "filter": 
    {
      "scale_wind_speed": {
        "Z_F": 4.6,
        "variable": "u"
      },
      "goodison_undercatch": {
        "variable": "p"
      }
    }
  }


Example
~~~~~~~~
.. code:: json

   "forcing": 
       {

         "UTC_offset": 8,

         "buckbrush": 
           {
             "file": "bb_m_2000-2008",
             "latitude": 60.52163,
             "longitude": -135.197151,
             "elevation": 1305,
             "filter":  
               {
               "scale_wind_speed": 
                   {
                   "Z_F": 4.6,
                   "variable": "u"
               },
               "goodison_undercatch":
               {
                   "variable":"p"
               }
            }
         },
         "alpine": 
           {
             "file": "alp_m_2000-2008",
             "latitude": 60.567267,
             "longitude": -135.184652,
             "elevation": 1559,
             "filter": {
            "scale_wind_speed": {
                "Z_F": 2.5,
                "variable": "u"
            },
            "goodison_undercatch":
            {
                "variable":"p"
            }

             }
         }
      }


NetCDF
********

The use NetCDF as input creates virtual stations at the cell-centres. The NetCDF file is lazy loaded as required for each triangle, so only the values required are loaded.
The variable names, like for ASCII inputs, needs to correspond to the values expected by the filters.


.. warning::
   
   NetCDF and ``point_mode`` are not supported.

Grid
~~~~~~

- the nc file is a regular grid of x,y values
- WGS84 lat/long
- consistent grid between model timesteps
- The underlying grid is specified in coordinates: ``ygrid_0`` and ``xgrid_0``
- The lat/long is specified in variables ``gridlat_0`` and ``gridlon_0``
- The elevation of the observation is given in ``HGT_P0_L1_GST`` (m)

Timesteps
~~~~~~~~~~

- at least two timesteps
- named ``datetime``
- time is in UTC+0
- the difference between these is used to determine model dt
- timesteps are offsets from an epoch (format ``YYYY-mm-dd HH:MM:SS`` or ``YYYY-mm-ddTHH:MM:SS``)
- units are hours, minutes, seconds
- This is specified as the units: ``datetime:units = "hours since 2017-09-01 06:00:00" ;``
- offset are given as ``int64``

Schema
~~~~~~~

In detail the following is the schema for the required NetCDF files:

.. code:: 

   dimensions:
      datetime = UNLIMITED ; 
      ygrid_0 = int ;
      xgrid_0 = int ;

   variables:
      double VAR_NAME(datetime, ygrid_0, xgrid_0) ;
         VAR_NAME:_FillValue = NaN ;
         VAR_NAME:coordinates = "gridlat_0 gridlon_0" ;

      double HGT_P0_L1_GST(datetime, ygrid_0, xgrid_0) ;
         HGT_P0_L1_GST:_FillValue = NaN ;
         HGT_P0_L1_GST:coordinates = "gridlat_0 gridlon_0" ;


      int64 datetime(datetime) ;
         datetime:standard_name = "time" ;
         datetime:long_name = "Validity time" ;
         datetime:axis = "T" ;
         datetime:units = "hours since 2017-09-01 06:00:00" ;
         datetime:calendar = "proleptic_gregorian" ;

      double gridlat_0(ygrid_0, xgrid_0) ;
         gridlat_0:_FillValue = NaN ;
         gridlat_0:long_name = "latitude" ;
         gridlat_0:standard_name = "latitude" ;
         gridlat_0:units = "degrees_north" ;

      double gridlon_0(ygrid_0, xgrid_0) ;
         gridlon_0:_FillValue = NaN ;


      double xgrid_0(xgrid_0) ;
         xgrid_0:_FillValue = NaN ;
         xgrid_0:long_name = "longitude in rotated pole grid" ;
         xgrid_0:standard_name = "grid_longitude" ;
         xgrid_0:units = "degrees" ;
         xgrid_0:axis = "X" ;

      double ygrid_0(ygrid_0) ;
         ygrid_0:_FillValue = NaN ;
         ygrid_0:long_name = "latitude in rotated pole grid" ;
         ygrid_0:standard_name = "grid_latitude" ;
         ygrid_0:units = "degrees" ;
         ygrid_0:axis = "Y" ;


Filters
~~~~~~~~

Filters are the same as for ASCII with one important distinction: every specified filter is run for every virtual station (i.e., grid cell centre).

.. code:: json

    "filter": {
               "scale_wind_speed": {
                   "Z_F": "40",
                   "variable": "u"
               }
           }


Example
~~~~~~~

.. code:: json

   "forcing": {
           "UTC_offset": "0",
           "use_netcdf": true,
           "file": "GEM-CHM_2p5_west_2017100106_2018080105.nc",
           "filter": {
               "scale_wind_speed": {
                   "Z_F": "40",
                   "variable": "u"
               }
           }
       }








