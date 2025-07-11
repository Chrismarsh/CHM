Output
=======

There are two main outputs from CHM: timeseries and mesh outputs.

The mesh outputs are either the Paraview vtu format or the netcdf ugrid format.

The vtu output has 1 file per MPI rank, per timestep.
For large domains, large MPI rank counts, and long time periods, this can produce a large number of files. HPC
file systems often do not perform well with millions of files. A further consideration is that to produce
raster tiff files or similar in post processing requires that the vtu files be converted to ugrid using pyCHM.
The vtu format does have some benifits as it allows for viewing a single ranks' subset, as well as visualisation of the
ghost regions. Lastly, because of limitations in the vtu format, a duplicate mesh used for outputting must be held in
memory making it more memory heavy than the ugrid output.


mesh (.vtu)
************

The ``.vtu`` format is a Paraview (vtk) unstructured mesh. The naming scheme of these files is 

   ``base_name`` + ``posix datetime`` + ``_MPIrank`` + ``.vtu``

For example: ``SC1506837600_0.vtu``

When running in MPI mode, only the part of the mesh being processed by the sub-MPI process is written out. This process' rank is thus suffixed to the file. If MPI is not used, a _0 will always be added for consistency. 

In addition to the vtu files, a ``base_name.pvd`` is written. This is an XML file that holds a reference to all
the vtu files:

.. code:: xml
   
   <?xml version="1.0" encoding="utf-8"?>
   <VTKFile type="Collection" version="0.1">
    <Collection>
        <DataSet timestep="1506837600" group="" part="0" file="SC1506837600_0.vtu"/> 
         ...

Although the ``vtu`` files may be loaded directly into Paraview, it is preferred to load the ``pvd`` file. Due to the ``timestep`` field, the :ref:`visualization:Datetime plugin` can then show an overlay with the date-time for easier analysis. 

.. note::

   If MPI is enabled, the ``pvd`` file is the only reasonable way of loading all the parts of the mesh into one view.

ugrid
*******

The ugrid is a netcdf standard https://ugrid-conventions.github.io/

For analysis, xarray can load this file, but uxarray is likely a better option

timeseries
***********

The timeseries output is a text, comma-separated-value (``.csv``) file. The first column is always ``datetime``, and the subsequent column names are the variables output from all the modules. This order is not guaranteed.

The datetime is in the format ``YYYYMMDDThhmmss``.

These files are output to the ``output_folder/points/`` subdirectory. The files are named ``station_name.txt``.

::

   datetime,ilwr,l,acc_snow
   20170901T060000,429.61,1.81749

