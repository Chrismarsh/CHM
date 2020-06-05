Getting Started
================

This page serves as a broad overview of 

For the general installation procedure, refer to the `build
page <https://github.com/Chrismarsh/CHM/wiki/Building-CHM>`__.

How it works
============

-  Configuration parser 
-  Terrain representation (raster data) 
-  Parameterization (triangulation of raster data)
-  Input filter 
-  Modular process 
-  Output/Data visualization

Configuration parser
--------------------

Every aspect of the model structure including initial conditions,
modules (and their parameters), forcing data, and output data can be
modified through either the command line or through a JSON configuration
file. Specifying parameters through the command line allows quick
on-the-fly testing without compromising some base configuration files.

A sample JSON configuration file is listed in Appendix.

Command line
~~~~~~~~~~~~

Specifying options through the `command
line <https://github.com/Chrismarsh/CHM/wiki/Command-line>`__ can be
done with either long- or short-form flags, i.e., –help, –version(-v),
–config-file(-f), –config(-c), –remove(-r), –remove-module(-d),
–add-module(-m).

“–config(-c)” overrides configurations in configuration files; however,
this parameter does not support list values, for example:

::

   -c config.Harder_precip_phase.const.b:1.5 -c config.debug.debug_level:"error"

“–remove(-r)” removes configurations, and it overrides configurations
specified by “–config”:

::

   -c nproc:2 -r nproc

“–remove-module(-d)” removes a module, and it overrides configurations
specified by “–config”:

::

   -d Marsh_shading_iswr

“–add-module(-m)” adds a module to the list:

::

   -m snobal -m Marsh_shading_iswr

Configuration files
~~~~~~~~~~~~~~~~~~~

`Configuration <https://github.com/Chrismarsh/CHM/wiki/Configuration>`__
for CHM is via a structured JSON file. The config file is structured
into “key:value” pairs separated by commas. Key names are enclosed in
quotes (“ ”), for example:

::

     {
       "key":
       {
         "value1":123
       },
       "key2":
       {
         "value2:"True"
       }
     }

Important notes:

1. Some configuration options are not compatible with others;

2. Do not prefix a number with zero (0). This is the octal prefix and it
   will cause the JSON parser to choke on lines that otherwise look
   fine;

3. Currently, ``mesher`` produces meshes that projected in Albers Conic
   equal area. Support for geographic exists in CHM, but there is an
   issue with computing surface normals;

4. Regardless of input coordinate system all input points are specified
   in latitude and longitude in WSG84.

Configurations include:

1. `option <https://github.com/Chrismarsh/CHM/wiki/Configuration#option>`__
   specifies additional options in the simulation;

2. `modules <https://github.com/Chrismarsh/CHM/wiki/Configuration#modules>`__
   specifies module components to load. Modules order as defined in this
   list has no bearing on the order they are run. Note modules are in a
   list (``[]``). Modules may be commented out to remove them from
   execution. Module names are case sensitive. The ``point_mode`` module
   is required to enable point mode, in addition to being enabled in
   ``option.point_mode``;

3. `remove_depency <https://github.com/Chrismarsh/CHM/wiki/Configuration#remove_depency>`__
   is to resolve circular dependencies among modules;

4. `config <https://github.com/Chrismarsh/CHM/wiki/Configuration#config>`__
   specifies module configurations;

5. `parameter_mapping <https://github.com/Chrismarsh/CHM/wiki/Configuration#parameter_mapping>`__
   specifies meta-data for on-mesh parameters;

6. `output <https://github.com/Chrismarsh/CHM/wiki/Configuration#output>`__
   specifies various output formats;

7. `global <https://github.com/Chrismarsh/CHM/wiki/Configuration#global>`__
   specifies a set of globally applicable parameters.

Mesh and ``mesher``
-------------------

General information on mesh representation refers to
`Mesh <https://github.com/Chrismarsh/CHM/wiki/Mesh>`__, `Mesh
generation <https://github.com/Chrismarsh/CHM/wiki/Mesh-generation>`__,
and an external tool for generating the meshes:
```Mesher`` <https://github.com/Chrismarsh/mesher>`__.

The mesh structure (``.mesh`` file) is formed as follows:

1. The triangle vertices are stored under the “vertex” key, e.g.,

   ::

           "vertex": [
           [
           488489.5,
           6713518.0,
           1525.4852294921875
           ]

2. Then a triangle is defined by indexing into that list of vertexes:

   ::

                  "elem": [
                  [
                  8033,
                  8160,
                  8043
                  ]

   So, the three edges of a triangle is from vertex 8033
   :math:`\rightarrow` vertex 8160; vertex 8160 :math:`\rightarrow`
   vertex 8043; vertex 8043 :math:`\rightarrow` vertex 8033.

3. Then “neigh” holds the neighbour topology:

   ::

                  "neigh": [
                  [
                  17687,
                  16277,
                  15812
                  ],

   | So triangle 0 has triangles 17687, 16277, and 15812 as neighbours.
   | If a triangle is an edge triangle, it’ll be missing a neighbour,
     denoted by -1:

   ::

                  [
                  -1,
                  19214,
                  11591
                  ],

Filters
-------

`Filters <https://github.com/Chrismarsh/CHM/wiki/Filters>`__ are a
mechanism whereby the input forcing data can be modified in some way
prior to the model run. For example, this could be use to apply a gauge
undercatch to precipitation. Filters modify the data of a station in
situ.

Note! Filters run in the order defined in the configuration file.

Input Timeseries
----------------

Time series data are input in a tab delimited format. Refer to
`Timeseries <https://github.com/Chrismarsh/CHM/wiki/Timeseries>`__ for
accepted format.

Modules and parallelization
---------------------------

`Modules <https://github.com/Chrismarsh/CHM/wiki/Modules>`__ are the
short-hand for a process representation. A principal design goal of a
module is that it may depend upon either some set of variables produced
by other modules or on input forcing data. Modules define a set of
variables which it provides globally to other modules. A module may not
overwrite a variable that another module declares. It should also not
overwrite the variables of another module. Implementation details on
modules can be found
`here <https://github.com/Chrismarsh/CHM/wiki/Modules#implementation-details>`__.

1. All ``module``\ s have pre-/post-conditions;

   Pre condition

   -  input forcing data or post-conditions from other ``module``\ s;

   Post condition

   -  see pre condition;

   Variables

   -  provide global variables to other ``module``\ s, but these
      variables do not change in other ``module``\ s.

2. There are two types of ``module``\ s:

   Forcing data interpolant

   -  depends upon point-scale input forcing data variables and
      interpolate these data onto every domain element;

   Standard process module

   -  depends only upon the output of interpolation ``module``\ s or
      other ``module``\ s’ output.

3. Parallelizations are offered in two ways, each module belongs to one
   of them:

   Data parallel

   -  point-scale models that are applied to every triangle;

   Domain parallel

   -  requires knowledge of surrounding mesh points.

   Parallelization process group ``module``\ s with same parallel type
   (data/domain) together and execute them simultaneously.

The class hierarchy of ``module`` looks like Figure [module_hier]:

Output Handling and Data visualization
--------------------------------------

Visualization is via `Paraview <https://www.paraview.org/>`__ if mesh
output is enabled in the configuration file. If ``PV_FILTER`` is enabled
in ``CMakeLists.txt``, a Paraview
`plugin <https://github.com/Chrismarsh/CHM/wiki/Visualization#datetime-plugin>`__
to show the date and time is built.

To convert the Paraview output (vtu files) to arbitrary GIS format,
refer to
`this <https://github.com/Chrismarsh/CHM/wiki/VTU-conversion>`__ page.

Known Issues
------------

1. CHM uses ``tcmalloc`` for memory allocation. There is a known
   deadlock
   `issue <https://github.com/gperftools/gperftools/issues/1037>`__ for
   ``tcmalloc`` used with ``gperftools`` (or any similar profiler, like
   Intel VTune). To resolve this, build CHM without ``tcmalloc``:
   ``-DUSE_TCMALLOC=OFF``

Resources
=========

TBA

Converting this document to PDF
===============================

To convert this document to PDF, follow the instructions below:

1. Install texlive, texlive-fonts-recommended, texlive-fonts-extra

2. Install pandoc (version > 2)

3. Copy ``eisvogel.latex`` to ``~/.pandoc/templates``

4. Execute
   ``pandoc CHM-tutorial.md -o CHM-tutorial.pdf --from markdown --template eisvogel --listings -V options``.
   Options can be found
   `here <https://pandoc.org/MANUAL.html#variables-for-latex>`__ and
   `here <https://github.com/Wandmalfarbe/pandoc-latex-template#custom-template-variables>`__.
   For example:

   ::

      pandoc CHM-tutorial.md -o CHM-tutorial.pdf --from markdown --template eisvogel --listings -V toc -V titlepage=true -V toc-own-page -V book -V title="CHM Tutorial"

A “Hello World” example
=======================

Source code (JSON file) of the example is listed as follows:

.. code:: json

         {
       
         "option":
         {
       
         // For point model to work, there must be an input and output station of the appropriate names. All other points will be ignored.
         // "point_mode":
         // {
         //   "output":"UpperClearing",
         //   "forcing":"UpperClearing"
         // },
       
         //      "notification_script":"./finished.sh",
         "per_triangle_timeseries":"false",
         "ui":false,
         "debug_level":"debug",
       
         "prj_name":"Marmot",
       
         "startdate":"20081001T140000",
         "enddate": "20081001T150000"
         //      "enddate":"20081001T000000"
         },
         "modules": //imporant these are [ ]
         [
         "solar",
         "iswr",
         "iswr_from_obs",
         // "point_mode",
         "Marsh_shading_iswr" // this is a domain parallel module
         //"fast_shadow" // this is a data parallel module
         // "scale_wind_vert",
       
         // "Harder_precip_phase",
       
         // "Sicart_ilwr",
         // "Walcek_cloud",
       
         //processes
         //    "snobal",
         //    "snowpack",
         // "Richard_albedo"
       
         ],
       
         // In case of a cycle depencency, remove dependencies between two modules.
         // If module A depends on B (A->B), then to remove the depency specify it as
         // "A":"B"
         // will remove the dependency on B from A.
         "remove_depency":
         {
         "scale_wind_vert":"snowpack",
         "scale_wind_vert":"snobal"
         },
         "config":
         {
         "Richard_albedo":
         {
         "min_swe_refresh":10,
         "init_albedo_snow":0.8
         },
         "point_mode":
         {
         "provide":
         {
         "iswr_diffuse":false,
         "iswr_direct":false,
         "iswr":false,
         "ilwr":false,
         "U_R":false,
         "vw_dir":false,
         "T_g":true
         }
       
         },
         "snobal":
         {
         "z_0":0.01
         },
         "snowpack":
         {
         "Snowpack":
         {
         "ROUGHNESS_LENGTH":0.01,
         "HEIGHT_OF_WIND_VALUE":2.96,
         "HEIGHT_OF_METEO_VALUES":2.6,
         "ATMOSPHERIC_STABILITY":"MO_MICHLMAYR"
         },
         "SnowpackAdvanced":
         {
         "ADJUST_HEIGHT_OF_WIND_VALUE":true,
         "ADJUST_HEIGHT_OF_METEO_VALUES":true,
         "HN_DENSITY":"MEASURED"
       
         }
         }
       
         },
         "meshes":
         {
         "mesh":"mesh/marmot1m.mesh"
         },
         "output":
         {
         "mesh":
         {
         "base_name":"shadow",
         "frequency":1
         }
       
         },
         "forcing":
         {
         "UTC_offset":6,
       
         "UpperClearing":
         {
         "file":"met/uc_2005_2018.txt",
         "longitude": -115.175362,
         "latitude":  50.956547,
         "elevation": 1844.6
         // "filter": {
         //    "scale_wind_speed": {
         //      "Z_F": 2,
         //      "variable": "u"
         //    }
         //  }
       
         }


         }
       
         }
