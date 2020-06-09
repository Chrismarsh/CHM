Mesh generation
=================

CHM uses a variable resolution, unstructured mesh. This mesh is generated using the `Mesher <https://mesher-hydro.readthedocs.io/en/latest/>`__ software. It needs to be installed separately from CHM. 


An example of a mesh is shown below

.. image:: images/usm.png 


The output ``.mesh``, ``.param``, and ``.ic`` from mesher, as `documented here <https://mesher-hydro.readthedocs.io/en/latest/output.html>`__ are used as input, without modification, to CHM.