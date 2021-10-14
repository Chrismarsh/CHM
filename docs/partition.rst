Partition tool
===============

The partition tool performs two functions:

1. Convert json (``.mesh``) format meshes to binary HDF5 (``.h5``)
2. Partition an HDF5 mesh format for *n* MPI ranks.

To use CHM in MPI mode *requires* using an HDF5 mesh file. When using CHM with MPI, the mesh is partitioned at run-time
for the appropriate number of MPI ranks. In this configuration CHM will  load the entire base mesh
(despite only operating on a small portion of it) into memory and load the
corresponding subset from the parameter file. As a result, the memory usage for each MPI rank is high.

The partition tool allows for pre-partitioning the mesh and parameter file into *n* chunks, one for each MPI rank.
Therefore only the mesh elements for this rank are loaded, dramatically reducing memory overhead.

.. warning::

   The mesh should be have been permuted prior to hdf5 conversion for maximum MPI performance.
   Failure to do so will result in significantly degraded runtime performance!

   Please see  :ref:`Mesh permutation <target mesh-permutation>`  for more details.



Usage
++++++

Options:

   - ``--help``
   - ``--mesh-file``, ``-m``,
   - ``--param-file``, ``-p``,
   - ``--max-ghost-distance``, ``-g``
   - ``--standalone``, ``-s``
   - ``--mpi-ranks``


``--help``
***********
Outputs the help message

``--mesh-file``
*****************
The mesh file to operate on. Either json or h5 format. However, formats cannot be mixed between ``mesh-file`` and ``param-file``.

.. code::

   --mesh-file <filename>
   --mesh-file my_mesh.mesh
   --mesh-file my_mesh.h5

``--param-file``
******************
A list of the parameter files. Either json or h5 format. However, formats cannot be mixed between ``mesh-file`` and ``param-file``.

.. code::

   --param-file <filename>
   --param-file my_mesh.param
   --param-file my_mesh_param.h5


``--max-ghost-distance``
**************************

The maximum distance at which to include ghost triangles. Certain modules need triangles at a distance, for example
to compute shadowing or fetch. ``--max-ghost-distance`` defaults to 100 m however this should be changed to match the
configuration option for various modules.

.. code::

   --max-ghost-distance <distance in meters>
   --max-ghost-distance 100

``--standalone``
******************
Produces a specific ranks' mesh in a way that allows it to be loaded by itself. This allows debugging on a
an individual rank. In order to be able to be loaded, the mesh is written **without** ghost regions. Zero-indexed.

.. code::

   --standalone <rank to output>
   --standalone 5

``--mpi-ranks``
*****************
Number of MPI ranks to split the mesh for. Ranks must be >1.

.. code::

   --mpi-ranks <n ranks>
   --mpi-ranks 32

Output
++++++++

.. note::

   If ``--mpi-ranks`` is specified and the input mesh is in json format (*.mesh*),
   the tool will, in two steps, convert the mesh to h5 and then partition the h5 mesh.

json to hdf5
*************
Output will be two ``.h5`` files called ``basename_mesh.h5`` and ``basename_param.h5``.
For example if the input mesh is ``granger1m.mesh`` then basename is ``granger1m`` resulting in the files

   - granger1m_mesh.h5
   - granger1m_param.h5

partition
**********
A json file is written ``basename_mesh.np<MPI_RANKS>.partition``. For example, the above file with 2 mpi ranks is
``granger1m_mesh.np2.partition``. The contents of this file describe the partition with various metadata and allow CHM
the ensure the file, when loaded, matches the runtime MPI configuration.

.. code:: json

   {
    "ranks": "2",
    "max_ghost_distance": "100",
    "num_global_faces": "37645",
    "meshes": [
        "granger1m_mesh.np2.partition.meshes\/granger1m_mesh.partition.0_mesh.h5",
        "granger1m_mesh.np2.partition.meshes\/granger1m_mesh.partition.1_mesh.h5"
    ],
    "parameters": [
        "granger1m_mesh.np2.partition.meshes\/granger1m_mesh.partition.0_param.h5",
        "granger1m_mesh.np2.partition.meshes\/granger1m_mesh.partition.1_param.h5"
    ]
   }

Secondly a new folder is created ``basename_mesh.np<MPI_RANKS>.partition.meshes`` that holds the underlying h5 mesh and param files.


Example Usage
++++++++++++++

.. code::

   ./partition --mesh-file granger1m.mesh --param-file granger1m.param # json to hdf5
   ./partition --mesh-file granger1m_mesh.h5 --param-file granger1m_param.h5 --mpi-ranks 2  # 2 MPI ranks
   ./partition --mesh-file granger1m_mesh.h5 --param-file granger1m_param.h5 --mpi-ranks 2 --standalone 1 # Make the 2nd rank standalone