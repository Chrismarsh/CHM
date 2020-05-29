Visualization
==============

Visualization is via `Paraview <http://www.paraview.org/>`__ if mesh
output is enabled in the configuration file.

Datetime plugin
===============

A `custom paraview
plugin <https://github.com/Chrismarsh/vtk-paraview-datetimefilter>`__ is
available to show the date and time. If ``BUILD_ParaView`` is enabled
paraview and the plugin are built.

To load the plugin in Paraivew:
``Tools -> Manage Plugins -> Load new -> Navigate to build directory``
After restarting Paraview, you will have to reload the plugin via
``Tools -> Manage Plugins -> Load Selected``

Optionally, copy the compiled filter into the plugins directory of
paraview and add ``<Plugin name="TimeAnnotate" auto_load="1"/>`` to the
``.plugins`` file. This will load the plugin automatically.

If mesh output is selected a pvd file, as well as multiple vtu files,
are generated. The pvd file is an XML file that links each output to the
julian date. This is required for showing the time.

To add the datetime filter to the view, load the pvd file and ensure
this is selected the left pane. Then, ``Filters->Search`` and search for
``datetime``.

Stations
========

For each forcing station/grid cell, a ``stations.vtp`` file is written
to the root of the output folder. The vtp files are a point dataset of
the x,y,z value of the forcing stations, as well as the station name as
a label. (See below for output point location file). To view the points
in Paraivew:

-  Load the vtp file
-  With the vtp file selected in the Pipeline Browser, choose Point
   Gaussian as the representation. Change the radius so the point is
   visible, or decrease it if it is too large.

To view the point labels: - Select vtp file in the pipeline browser -
Create a new spreadhseet layout - Select the points you wish to have
labels displayed for - View->Selection Display Inspector - Choose Point
Labels drop down and select ‘Station name’. - Use the cog next to
‘selection colour’ to change the display font (size, colour, etc)

|image0|

Output points
=============

If single triangle point-output is selected, these points are written to
a seperate vtp file (``output_points.vtp``) in the output/points
directory. To view, follow the above directions.

Remote visualization on Graham
==============================

A pvserver can be hosted on Compute Canada’s Graham cluster, allowing
for remote visualization without the need to copy files locally. CC
notes this is still experimental and crashes do occur. However, it seems
to generally work well. Instructions can be found
`here <https://docs.computecanada.ca/wiki/ParaView>`__. Contrary to
documentation, some users have found better performance not using
``--mesa-swr-avx2`` and some stability improvements were found too. It
is worth trying with and without this option.

For larger domains and complex multi-views, using the GPU nodes appears
to be more stable.

Note that the local Paraview client version must *exactly* match the pvserver version on Graham.
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. |image0| image:: images/viz_points.png
