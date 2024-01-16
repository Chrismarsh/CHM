Compilation
============

.. warning::
    Conan is no longer used to build CHM. Spack is now used



CHM uses `spack <https://spack.readthedocs.io/>`__ to manage and build all
dependencies. Because of the various requirements on build
configuration, versions, and inter-dependencies, using system libraries
it not supported.


Build requirements
*******************

Linux (x86_64) and Macos (arm64) are the only supported environments.

The following have been tested

=======  =====  ========  =====
   Linux          Macos
--------------  ---------------
Ubuntu   18.04
  -      20.04
  -      22.04  Ventura  13
=======  =====  ========  =====        

Only gcc is currently support:
gcc (libc 2.27+): 9.3.0+


.. warning::
   Below, references to building with the Intel compiler are made. Unfortunately the Intel compiler doesn't currently work with applications that also
   link against GSL. This is being investigated. For now, please do no build CHM with Intel Compilers.

If using conan to build the dependencies, the only requirements are:

   - conan (via Python pip)
   - cmake >=3.21  (via apt-get/brew)
   - C++14 compiler (gcc 9.3.0+) (via apt-get/brew)
   - Fortran 90+ compiler (gfortran) (via apt-get/brew)






On MacOS, `homebrew <https://brew.sh/>`__ should be used to install
cmake. Macport based installs likely work, but have not been
tested.

Intel compiler
---------------

.. warning::
   As noted above, Intel compiler builds are currently not supported and do not work. Please use gcc at this time.

If the Intel compiler is being used (this is optional), ensure the Intel compilervars is sourced, e.g.,

::

   source /opt/intel/bin/compilervars.sh intel64

prior to running the conan. Use the gcc settings for conan.

Build dependencies
*********************

OpenMP
~~~~~~~~

On MacOS, the openmp library should be installed via homebrew:

::

   brew install libomp


OpenMPI
~~~~~~~~~
.. note::

   No MPI support is now deprecated. MPI is now required.


Other MPI versions probably work, but only OpenMPI has been tested.
      Ensure MPI is installed on linux:
      ``libopenmpi-dev`` and ``openmpi-bin``




Install spack
-------------

Install `spack <https://spack-tutorial.readthedocs.io/en/latest/tutorial_basics.html>`__.

Clone the CHM spac-repo

::

   git clone https://github.com/Chrismarsh/spack-repo.git /some/path/here/


then create `repos.yml` in `~/.spack` and add the path to the above cloned `spack-repo`.
It will look like this

::

    $ cat ~/.spack/repos.yaml
    	repos:
    	  - /some/path/here//spack-repo
    	  - $spack/var/spack/repos/builtin


Setup CHM source folders
------------------------

An out of source build should be used. That is, build in a separate folder removed from the CHM source. This makes it easier to clean up
and start from scratch. An example is given below:

::

   cd ~/
   git clone --recurse-submodules https://github.com/Chrismarsh/CHM

   mkdir ~/build-CHM

Build dependencies
---------------------

This step will build and install the dependencies via spack.

::

    spack env create chm path/to/CHM/spack.yml
    spack concretize
    spack install -j 8 # number of parallel builds, adjust accordingly


Intel MKL
~~~~~~~~~

.. warning::
   Using MKL with Trilinos is not supported as the final CHM link will conflict with the internal BLAS in GSL.


Build CHM
***********

.. note::
   The follow instructions assume that they are invoked from within ``~/build-CHM`` (or your equivalent).


Run cmake
---------

You can set the install prefix to be anywhere, such as shown in the
example below

::

   cmake ~/CHM -DCMAKE_INSTALL_PREFIX=/opt/chm-install

Both ``ninja`` and ``make``
(this is the default) are supported. To use ``ninja``, add

::

   cmake ~/CHM -DCMAKE_INSTALL_PREFIX=/opt/chm-install -G "Ninja"

Ninja speeds up compilation of CHM by ~6%.

The default build option creates an optimizted “release” build. To build
a debug build, use ``-DCMAKE_BUILD_TYPE=Debug``.




Intel compiler
~~~~~~~~~~~~~~

If the Intel compiler is used, add the following cmake flags:

::

   -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc -DCMAKE_FORTRAN_COMPILER=ifort

High performance allocators
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default jemalloc is used.

``-DUSE_TCMALLOC=FALSE -DUSE_JECMALLOC=TRUE``.

Building
--------

Using make

::

   make -jN CHM

where N is the number of parallel jobs (e.g., total core count).

Using Ninja

::

   ninja -C . 

Run tests
---------

Tests can be enabled with ``-DBUILD_TESTS=TRUE`` and run with
``make check``/ ``ninja check``

Install
-------

``make install``/``ninja install``

Build docs
***********
To build the documentation requires `Doxygen <https://www.doxygen.nl/download.html>`__ and Sphinx+Breathe+Exhale.

.. code::

   pip install sphinx
   pip install sphinx-rtd-theme
   pip install breathe<4.13.0
   pip install exhale

The Breathe version requirement is for Read the Docs compatibility. See `issue#89 <https://github.com/svenevs/exhale/issues/89>`__.

The documentation can be built with:

::

   cd CHM/docs
   READTHEDOCS="True" make html


The env var is required to ensure the correct directories are searched for in-source builds. 


Troubleshooting
***************

Disable allocators
--------------------

The high performance allocators may need to be disabled and can be done via
``-DUSE_TCMALLOC=FALSE -DUSE_JEMALLOC=FALSE``



Building on Compute Canada (WestGrid)
******************************************

To build on Compute Canada stack machines, such as Graham, all dependencies must be built
from source to ensure the correct optimizations are used. This should be done with the Compute Canada easybuild system.

Only the ``gcc/9.3.0`` environment is supported. This can be enabled with

::

   module load gcc/9.3.0


easybuild
-----------

Build all dependencies that are not available from compute canada stack

::

   git clone https://github.com/Chrismarsh/easy_build.git
   cd easy_build
   chmod +x install-all.sh
   ./install-all.sh

Building CHM
------------

Ensure the environment is correctly setup

::

   module load armadillo/10.4.1
   module load cgal/5.2.1
   module load hdf5/1.10.6
   module load meteoio
   module load func
   module load netcdf/4.7.4
   module load gdal/3.2.3
   module load boost-mpi/1.72.0
   module load openblas/0.3.9
   module load gsl/2.6
   module load eigen/3.3.7
   module load sparsehash
   module load tbb/2020.2
   module load trilinos/13.3.0
   module load netcdf-c++4/4.3.1
   module load vtk/9.0.1
   module load proj/9.0.1
   module load jemalloc/5.3.0
   module load cmake

Optionally you can save this with ``module save chm``.


Then build CHM

::

   git clone https://github.com/Chrismarsh/CHM  # get CHM source code
   mkdir ~/chm-build && cd ~/chm-build # make build directory
   cmake ../CHM -DENABLE_SAFE_CHECKS=ON -DBoost_NO_BOOST_CMAKE=ON -DUSE_TCMALLOC=FALSE -DUSE_JEMALLOC=TRUE -DCMAKE_BUILD_TYPE=Release
   make -j10


