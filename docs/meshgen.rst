Mesh generation
=================

CHM uses a variable resolution, unstructured mesh. This mesh is generated using the `Mesher <https://mesher-hydro.readthedocs.io/en/latest/>`__ software. It needs to be installed separately from CHM. 


An example of a mesh is shown below

.. image:: images/usm.png 


The output ``.mesh``, ``.param``, and ``.ic`` from mesher, as `documented here <https://mesher-hydro.readthedocs.io/en/latest/output.html>`__ are used as input, without modification, to CHM.


MPI compatibility
-------------------
In MPI mode, CHM needs to be able to only read part of the mesh on a per-MPI rank basis. Therefore a HDF5 mesh is used.
Mesher does not produce this format. However, CHM can convert from the Mesher json format to HDF5.

Conversion
++++++++++
1. Create mesh as normal with mesher

2. Create the CHM configuration mesh section using the .mesh and .param files as normal

.. code:: json

   "meshes":
   {
    "mesh":"meshes/granger.mesh",
    "parameters":
    {
      "file":"meshes/granger.param",
    }
   }

3. Run CHM with the flag `--convert-hdf5`. E.g.,

.. code::

   ./CHM -f config.json --convert-hdf5

CHM will run and then terminate after the conversion. This flag will produce .h5 files named `<original_name>.mesh_mesh.h5` and `<original_name>.mesh_param.h5`.

4. Edit configuration to use these files

.. code:: json

   "meshes":
   {
    "mesh":"meshes/granger.mesh_mesh.h5",
    "parameters":
    {
      "file":"meshes/granger.mesh_param.h5,
    }
   }

.. warning::

   HDF5 and json meshes/param files **cannot** be mixed.