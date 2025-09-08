Installation
==============

The simplest way to build CHM for usage is to install and configure spack (as described below :ref:`configure-spack`)
and then install CHM:

::

    spack install chm


Environment requirements
**************************

.. note::
    Conan is no longer used to build CHM. Spack is now used

.. note::
   Building CHM without MPI support is now deprecated. MPI is now required.

.. note::
    The OpenMP features are disabled by default for good reason. Unless you know you need OpenMP support,
    use MPI to obtain good parallel performance.

Linux (x86_64) and Macos (arm64) are the only supported environments.

.. warning::
   Unfortunately the Intel compiler doesn't currently work with applications that also
   link against GSL. This is being investigated. For now, please do not build CHM with Intel Compilers.

It is recommended to use `spack <https://spack.readthedocs.io/>`_ to manage and build all
dependencies. Because of the various requirements on build
configuration, versions, and inter-dependencies, using system libraries (apt/yum/brew/&c)
is not recommended. Take care when compiling against homebrew libraries. Because homebrew releases
new versions of libraries often, it often results in having to frequently recompile CHM.

Spack will build all required libraries and their dependencies, including compilers and MPI as required.
This is the recommended approach.

As the build system uses cmake to locate libraries, there are no assumptions about using spack, so any
library provider will work, such as the above noted system libraries or other dependency management tools
like easy_build.


CHM with spack
***************

This is the recommend method to build CHM for HPC usage.


Install spack
+++++++++++++++
Install `spack <https://spack-tutorial.readthedocs.io/en/latest/tutorial_basics.html>`__

Use the git repository and use the develop branch, as significant bug fixes to packages CHM uses have been made in
this branch.

.. _configure-spack:

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

On macos, apple-clang does not have a fortran compiler. Install gcc via spack and gfortran will be detected and
added to apple-clang as the FC provider when using spack compiler find. Using brew's gfortran works, but you're
then at the mercy of brew updates.

CHM uses some libraries that are not currently in mainline spack. Please add the CHM spack repo as:

::

   spack repo add https://github.com/Chrismarsh/spack-repo.git


You may also manually clone spack-repo and add it to ``~/.spack/repos.yaml``:

::

    repos:
      chm: $spack/../spack-repo/spack_repo/chm


Build and install CHM
++++++++++++++++++++++++

Once everything for spack is setup:

::

    spack install chm


(optional) Create develop environment
++++++++++++++++++++++++++++++++++++++++
If you wish to develop CHM, then the following should be done to install the dependencies and build an environment
to work on CHM. If you only want to run CHM, use the method above in :ref:`Installation`.

Create a spack environment and build and install the dependencies.

::

    spack env create chm ~/CHM/spack.yaml
    spack env activate chm
    spack concretize
    spack install -j 8 # number of parallel builds, adjust accordingly


CHM with easy_build
**********************
.. warning::
    This hasn't been updated in a while. YMMV.

When targeting the Digital Alliance Canada stack, `this repository <https://github.com/Chrismarsh/easy_build>`__ hosts
the easy_build scripts needed for missing libraries. They can be installed in a dependency-preserving order with ``install-all.sh``.

Digital Alliance Canada
++++++++++++++++++++++++

To build on Compute Canada stack machines, such as Graham, all dependencies must be built
from source to ensure the correct optimizations are used.

Only the ``gcc/9.3.0`` compiler in ``StdEnv/2020`` is tested. This can be enabled with

::

   module load gcc/9.3.0



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

.. warning::
    This often doesn't work well due to system libraries not being compiled with the features needed. YMMV.

Please install the the libraries from your package manager of choice following the version constraints listed
in the `CHM spack.yaml <https://github.com/Chrismarsh/CHM/blob/develop/spack.yaml>`_ file.

.. warning::
    apple-clang doesn't ship with an OpenMP library,
    so OpenMP `should be installed <https://mac.r-project.org/openmp/>`__.
    Doing so via homebrew (``brew install libomp``) is likely the easiest.

and then follow the instructions for CHM standalone.

CHM standalone
****************

Regardless of what method was used to build the libraries, the configuration of CHM is the same.

Build env requirements:
   - cmake >=3.30
   - C++20 compiler (e.g., gcc 9.3.0+)
   - Fortran 90+ compiler (e.g., gfortran)
   - OpenMPI or IntelMPI

Source code
+++++++++++++++

An out of source build should be used. That is, build in a separate folder outside of the CHM source.
This makes it easier to clean up and start from scratch and to keep separate release and debug builds.

An example is given below:

::

   cd ~/
   git clone https://github.com/Chrismarsh/CHM
   mkdir ~/build-CHM
   cd ~/build-CHM
   # This is where the build configuration will occur in the next steps


.. note::
   The following instructions assume that they are invoked from within ``~/build-CHM`` (or your equivalent).


Run cmake
++++++++++++

This guide assumes you are building CHM to build and debug it. However, you can set the install prefix to be anywhere,
such as shown in the example below

::

   cmake ~/CHM -DCMAKE_INSTALL_PREFIX=/opt/chm-install

Both ``ninja`` and ``make`` (this is the default) are supported. To use ``ninja``, add

::

   cmake ~/CHM -DCMAKE_INSTALL_PREFIX=/opt/chm-install -G "Ninja"


The default build option creates an optimized “release” build. To build
a debug build, use ``-DCMAKE_BUILD_TYPE=Debug``.


Set compiler
+++++++++++++++

CMake does not always detect the most-optimal compiler you wish to use. The compiler can be manually specified to cmake.
For example, if the Intel compiler is used, add the following cmake flags:

::

   -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc -DCMAKE_FORTRAN_COMPILER=ifort

If using spack to build the compiler (e.g., gcc), use ``spack find -p gcc`` to find the spack-built compiler.

High performance allocators
++++++++++++++++++++++++++++++++

By default jemalloc is used.

``-DUSE_TCMALLOC=FALSE -DUSE_JECMALLOC=TRUE``.


Building
+++++++++++

Using make

::

   make -jN

where N is the number of parallel jobs (e.g., total core count).

Using Ninja

::

   ninja -C . 

Run tests
++++++++++

Tests can be enabled with ``-DBUILD_TESTS=TRUE`` and run with
``make check``/ ``ninja check``. These have not been updated and currently fail

Install
+++++++++++

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
+++++++++++++++++++++

The high performance allocators may need to be disabled and can be done via
``-DUSE_TCMALLOC=FALSE -DUSE_JEMALLOC=FALSE``





