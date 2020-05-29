Timeseries
==========

ASCII file forcing data
=======================

Time series data are input in a tab delimitated format. The only
required column name is ``datetime``. Datetime format is ISO standard as
[year][month][day]T[hour][minute][second]. Column order does not matter.
The model may start at any time, however there is an assumption of (1)
constant time stepping and (2) all forcing files have the same start and
end dates. The difference in time_0 and time_1 is used to determine the
internal model timesteps. The variable names in the columns must be the
same names as what the interpolation modules expect. Missing values are
permitted and should be set to -9999. Timezone should match the
“UTC_offset” parameter in the configuration file.

::

   datetime       Qsi      Qli    g    t       rh      u       vw_dir  p
   20001001T000000 -0.237  276 -2.436  -10.98  95.7    2.599   308.2   0.0
   20001001T003000 -0.233  278 -2.42   -11.19  95.7    3.133   307.1   0.0

The forcing units are:

::

   datetime       Qsi      Qli    g    t       rh      u       vw_dir  p
   20001001T003000 (Wm/2)  (W/m2)  (Wm/2)  (C) (%) (m/s)   (degrees (clockwise from true north))   (mm/(time step))

This format is easily parseable with Pandas in Python

.. code:: python

   obs = pd.read_csv("uc_2005_2014.txt",sep="\t",parse_dates=[0])
   obs.set_index('datetime',inplace=True)
   obs.index = pd.to_datetime(obs.index)

Various conversion scripts for other models’ input/output are located in
``tools``

These are specific in the forcing section of the configuration file

::

       "station_name": {
           "elevation": 1580.0,
           "file": "05BA801_out.txt",
           "latitude": 51.416667,
           "longitude": -116.183333
       }

NetCDF
======

NetCDF files are specific on a lat/long grid. These are added as
follows. Filters can be run on them, but the same filter is run across
the entire netcdf grid

``"forcing":   {         "UTC_offset":0,         "use_netcdf":true,         "file":"GEM-CHM_2p5_snowcast_2017090106_2018083005.nc",         "filter": {             "scale_wind_speed": {                 "Z_F": 40,                 "variable": "u"        }     }   }``
