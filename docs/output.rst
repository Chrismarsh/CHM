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

vtu
************

The ``.vtu`` format is a Paraview (vtk) unstructured mesh. It is one file, per MPI rank, per timestep.

The naming scheme of these files is

   ``base_name`` + ``posix datetime`` + ``_MPIrank`` + ``.vtu``

For example: ``SC1506837600_0.vtu``

Each MPI rank writes its own subdomain into a single file. This process' rank is thus suffixed to the file.
If only 1 rank is used, a _0 will always be added for consistency.

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

   If more than 1 rank is used, the ``pvd`` file is the only reasonable way of loading all the parts of the mesh into one view.

ugrid
*******

The ugrid file is a netcdf `URGID <https://ugrid-conventions.github.io/>`_  schema. It defines the mesh in the same
 way as what is used internally to CHM where the triangles are composed indexes of 3-nodes into the vertex list.
 The values in CHM are almost exclusively face centered.

.. code::

    $ ncdump -h output_casr3_debug.nc
    netcdf output_casr3_debug {
    dimensions:
            nMesh2_node = 3761 ;
            nMesh2_face = 6828 ;
            Two = 2 ;
            Three = 3 ;
            time = 24 ;
    variables:
            double time(time) ;
                    time:standard_name = "time" ;
                    time:long_name = "Time" ;
                    time:units = "minutes since 1970-01-01 00:00:00" ;
            int Mesh2 ;
                    Mesh2:cf_role = "mesh_topology" ;
                    Mesh2:long_name = "Topology data of 2D unstructured mesh" ;
                    Mesh2:topology_dimension = 2 ;
                    Mesh2:node_coordinates = "Mesh2_node_x Mesh2_node_y" ;
                    Mesh2:face_node_connectivity = "Mesh2_face_nodes" ;
                    Mesh2:face_dimension = "nMesh2_face" ;
                    Mesh2:face_coordinates = "Mesh2_face_x Mesh2_face_y" ;
            uint Mesh2_face_nodes(nMesh2_face, Three) ;
                    Mesh2_face_nodes:cf_role = "face_node_connectivity" ;
                    Mesh2_face_nodes:long_name = "Maps every triangular face to its three corner nodes." ;
            double Mesh2_node_x(nMesh2_node) ;
                    Mesh2_node_x:standard_name = "longitude" ;
                    Mesh2_node_x:long_name = "Longitude of 2D mesh nodes." ;
                    Mesh2_node_x:units = "degrees_east" ;
                    Mesh2_node_x:_FillValue = NaN ;
            double Mesh2_node_y(nMesh2_node) ;
                    Mesh2_node_y:standard_name = "latitude" ;
                    Mesh2_node_y:long_name = "Latitude of 2D mesh nodes." ;
                    Mesh2_node_y:units = "degrees_north" ;
                    Mesh2_node_y:_FillValue = NaN ;
            double Mesh2_node_z(nMesh2_node) ;
                    Mesh2_node_z:standard_name = "altitude" ;
                    Mesh2_node_z:long_name = "Z coordinate of 2D mesh nodes." ;
                    Mesh2_node_z:units = "m" ;
                    Mesh2_node_z:mesh = "Mesh2" ;
                    Mesh2_node_z:location = "node" ;
                    Mesh2_node_z:_FillValue = NaN ;
            double Mesh2_node_z_paraview(nMesh2_node) ;
                    Mesh2_node_z_paraview:standard_name = "scaled_altitude" ;
                    Mesh2_node_z_paraview:long_name = "Scaled Z coordinate of 2D mesh nodes." ;
                    Mesh2_node_z_paraview:units = "m" ;
                    Mesh2_node_z_paraview:mesh = "Mesh2" ;
                    Mesh2_node_z_paraview:location = "node" ;
                    Mesh2_node_z_paraview:_FillValue = NaN ;
            uint64 global_id(nMesh2_face) ;
                    global_id:mesh = "Mesh2" ;
                    global_id:location = "face" ;
                    global_id:coordinates = "Mesh2_face_x Mesh2_face_y" ;
            double Mesh2_face_x(nMesh2_face) ;
                    Mesh2_face_x:standard_name = "latitude" ;
                    Mesh2_face_x:long_name = "Characteristics latitude of 2D mesh triangle (e.g. circumcenter coordinate)." ;
                    Mesh2_face_x:units = "degrees_north" ;
            double Mesh2_face_y(nMesh2_face) ;
                    Mesh2_face_y:standard_name = "latitude" ;
                    Mesh2_face_y:long_name = "Characteristics latitude of 2D mesh triangle (e.g. circumcenter coordinate)." ;
                    Mesh2_face_y:units = "degrees_north" ;
            double Mesh2_face_z(nMesh2_face) ;
                    Mesh2_face_z:standard_name = "altitude" ;
                    Mesh2_face_z:long_name = "Characteristics latitude of 2D mesh triangle (e.g. circumcenter coordinate)." ;
                    Mesh2_face_z:units = "m" ;
            double rh(time, nMesh2_face) ;
                    rh:mesh = "Mesh2" ;
                    rh:location = "face" ;
                    rh:_FillValue = NaN ;
            double Elevation(nMesh2_face) ;
                    Elevation:mesh = "Mesh2" ;
                    Elevation:location = "face" ;
                    Elevation:_FillValue = NaN ;


Variables are defined on ``dimensions = (face, time)``, and parameters are only on ``dimension = (face)``.


This ugrid file can be loaded into Paraview (chose the netcdf ugrid loader option in Paraview). By default the z coordinate is not shown.
From the Filters menu option, search for "Warp by Scalar" and choose the ``Mesh2_node_z_paraview`` field. Paraview struggles with lat/long x,y coords with
high elevation values thus, this is a scaled-for-visualization value. If you want to see the cell elevation, choose the "Elevation" from the variable colouring menu.

For analysis, xarray can load ugrid, but uxarray is likely a better option

timeseries
***********

The timeseries output is a text, comma-separated-value (``.csv``) file. The first column is always ``datetime``, and the subsequent column names are the variables output from all the modules. This order is not guaranteed.

The datetime is in the format ``YYYYMMDDThhmmss``.

These files are output to the ``output_folder/points/`` subdirectory. The files are named ``station_name.txt``.

::

   datetime,ilwr,l,acc_snow
   20170901T060000,429.61,1.81749

