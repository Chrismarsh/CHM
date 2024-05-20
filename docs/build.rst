Compilation
============

.. warning::
    Conan is no longer used to build CHM. Spack is now used

.. note::
   Building CHM without MPI support is now deprecated. MPI is now required.

CHM uses `spack <https://spack.readthedocs.io/>`__ to manage and build all
dependencies. Because of the various requirements on build
configuration, versions, and inter-dependencies, using system libraries (apt/yum/brew/&c)
it not recommended. Take care when compiling against homebrew libraries. Because homebrew releases
new versions of libraries often, it often results in having to frequently recompile CHM.

However, as the build system uses cmake to locate libraries, there are no assumptions about using spack, so any
library provider will work, such as the above noted system libraries or other dependency management tools
like easy_build.

Environment requirements
**************************

Linux (x86_64) and Macos (arm64) are the only supported environments.

Build env requirements:
   - cmake >=3.21
   - C++20 compiler (e.g., gcc 9.3.0+)
   - Fortran 90+ compiler (e.g., gfortran)
   - OpenMPI or IntelMPI

.. warning::
   Unfortunately the Intel compiler doesn't currently work with applications that also
   link against GSL. This is being investigated. For now, please do no build CHM with Intel Compilers.


Spack will build all required libraries and their dependencies, including compilers and MPI as required.
This is the recommended approach.

CHM source
*************

An out of source build should be used. That is, build in a separate folder outside of the CHM source.
This makes it easier to clean up and start from scratch and to keep seperate release and debug builds.

An example is given below:

::

   cd ~/
   git clone --recurse-submodules https://github.com/Chrismarsh/CHM
   mkdir ~/build-CHM
   cd ~/build-CHM
   # This is where the build configuration will occur in the next steps



CHM with spack
***************

Install spack
+++++++++++++++
Install `spack <https://spack-tutorial.readthedocs.io/en/latest/tutorial_basics.html>`__

Use the git repository and use the develop branch, as significant bug fixes to packages CHM uses have been made in
this branch.

Configure Spack
+++++++++++++++++++
It is critical to ensure spack is correctly
configured, as described in the `Spack Getting Started <https://spack.readthedocs.io/en/latest/getting_started.html>`__
guide.

If you need to build a compiler via spack to use to build CHM and the spack libraries, this is the time to do it.
Otherwise, ensure
the `external compiler is found by
spack <https://spack.readthedocs.io/en/latest/getting_started.html#spack-compiler-find>`__ and correctly configured.

If you use a system MPI or intel-oneapi-(mkl|tbb) (i.e., a not-spack built version), this is when it should be configured
`as a spack external <https://spack.readthedocs.io/en/latest/packages_yaml.html#external-packages>`__ package.

On macos, apple-clang does not have a fortran compiler. The current suggestion is to install gfortran via gcc in brew
and it will be detected as the fortran compiler for apple-clang in spack. Note th at when gcc updates, you'll need to
rebuild impacted packages.

CHM uses some libraries that are not currently in mainline spack. Until they have been accepted, please clone
the CHM spack-repo

::

   git clone https://github.com/Chrismarsh/spack-repo.git /some/path/here/


then create ``repos.yml`` in ``~/.spack`` and add the path to the above cloned ``spack-repo``.
It will look like this

::

    $ cat ~/.spack/repos.yaml
    	repos:
    	  - /some/path/here/spack-repo
    	  - $spack/var/spack/repos/builtin


Build dependencies
+++++++++++++++++++++

This step will build and install the dependencies via spack.

::

    spack env create chm ~/CHM/spack.yaml
    spack env activate chm
    spack concretize
    spack install -j 8 # number of parallel builds, adjust accordingly


CHM with easy_build
**********************
When targetting the Digital Alliance Canada stack, `this repository <https://github.com/Chrismarsh/easy_build>`__ hosts
the easy_build scripts needed for missing libraries. They can be installed in a depdency-preserving order with ``install-all.sh``.

Digitial Alliance Canada
-----------------------------

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
   module load trilinos/15.0.0
   module load netcdf-c++4/4.3.1
   module load vtk/9.0.1
   module load proj/9.0.1
   module load jemalloc/5.3.0
   module load cmake

Optionally you can save this with ``module save chm``.


CHM with system packages
**************************

Please install the following libraries using the package manager of your choice:
  - boost >= 1.74.0 with: system, filesystem, date_time, thread, regex, iostreams, program_options, mpi, serialization
  - cgal (header-only)
  - hdf5 with c++ bindings
  - netcdf
  - netcdf-cxx4 >= 4.3
  - gdal >= 3.6
  - proj >=9
  - sparsehash
  - gperftools (only needed if tcmalloc is enabled)
  - gsl
  - armadillo
  - intel-oneapi-tbb
  - eigen
  - meteoio
  - func
  - trilinos@15.0.0 with mpi, optionally openmp & threadsafe if CHM is built with omp
  - jemalloc (only needed if jemalloc is enabled)
  - vtk >= 9.2
  - spdlog
  - openblas
  - MPI

.. warning::
    apple-clang doesn't ship with an OpenMP library,
    so OpenMP `should be installed <https://mac.r-project.org/openmp/>`__.
    Doing so via homebrew (``brew install libomp``) is likely the easiest.

Compile CHM
**************

Regardless of what method was used to build the libraries, the configuration of CHM is the same.

.. note::
   The follow instructions assume that they are invoked from within ``~/build-CHM`` (or your equivalent).


Run cmake
---------

This guide assumes you are building CHM to build and debug it. However, you can set the install prefix to be anywhere,
such as shown in the example below

::

   cmake ~/CHM -DCMAKE_INSTALL_PREFIX=/opt/chm-install

Both ``ninja`` and ``make`` (this is the default) are supported. To use ``ninja``, add

::

   cmake ~/CHM -DCMAKE_INSTALL_PREFIX=/opt/chm-install -G "Ninja"


The default build option creates an optimizted “release” build. To build
a debug build, use ``-DCMAKE_BUILD_TYPE=Debug``.


Set compiler
~~~~~~~~~~~~~~

CMake does not always detect the most-optimal compiler you wish to use. The compiler can be manually specified to cmake.
For example, if the Intel compiler is used, add the following cmake flags:

::

   -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc -DCMAKE_FORTRAN_COMPILER=ifort

If using spack to build the compiler (e.g., gcc), use ``spack find -p gcc`` to find the spack-built compiler.

High performance allocators
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default jemalloc is used.

``-DUSE_TCMALLOC=FALSE -DUSE_JECMALLOC=TRUE``.


Building
--------

Using make

::

   make -jN

where N is the number of parallel jobs (e.g., total core count).

Using Ninja

::

   ninja -C . 

Run tests
---------

Tests can be enabled with ``-DBUILD_TESTS=TRUE`` and run with
``make check``/ ``ninja check``. These have not been updated and currently fail

Install
-------

``make install``/``ninja install``

Build docs
***********
To build the documentation requires `Doxygen <https://www.doxygen.nl/download.html>`__ and Sphinx+Breathe+Exhale.

.. code::

   pip install sphinx
   pip install sphinx-rtd-theme
   pip install breathe
   pip install exhale@git+https://github.com/svenevs/exhale.git@refs/pull/205/merge

The exhale version requirement: see `issue#200 <https://github.com/svenevs/exhale/pull/200>`__.

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

Only the ``gcc/9.3.0`` compiler in ``StdEnv/2020`` is supported. This can be enabled with

::

   module load gcc/9.3.0



