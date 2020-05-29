mesh
=====
The internal topographic representation is via an unstructured
triangular mesh (herein ‘the mesh’). An example of what this looks like
is below.

|image0|

Internally, the mesh is held in a CGAL structure that provides various
ease-of-use structures. The relevant files are in ``mesh/``.

Iteration
=========

Because the triangle iterators provided by CGAL have a non-deterministic
order, as well as being incompatible with OpenMP, the way to access the
i-th triangle is via

.. code:: cpp

   #pragma omp parallel for
       for (size_t i = 0; i < domain->size_faces(); i++)
       {
           auto elem = domain->face(i);
        ...
       }

Internal data structure
=======================

The mesh is represented as a vector of (x,y,z) points, connected via
edges which are indexes into the vertex array. Thus, two vectors are
held in memory as (pseudo code)

.. code:: python

   vertexes = [ 
     [x1,y1,z1],
     [x2,y2,z2],
     [x3,y3,z3]
   ]

   triangles = [
     [1,2,3]
   ]

Above, this represents a single triangle.

Input
=====

The input mesh is produced using the
`mesher <https://github.com/Chrismarsh/mesher>`__ tool. Please see the
mesher documentation for more details.

.. |image0| image:: images/mesh.png
