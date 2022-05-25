Compilation
============

CHM uses `conan <https://conan.io/>`__ to manage and build all
dependencies. Because of the various requirements on build
configuration, versions, and interdependencies, using system libraries
it not supported.

All of the CHM dependencies are built on Travis-CI and uploaded to the an Artifactory repository to serve
prebuilt binaries and build scripts. This means that *if* the CHM build is done with
supported compilers and operating system (described later), the
dependencies do not need to be built by the end user. However it is generally recommended to build all dependencies
from source and this is especially the case for high-performance environments and MPI.

Build requirements
*******************

Linux and MacOS are the only supported environments. The following have been tested

=======  =====  ========  =====
   Linux          MacOS
--------------  ---------------
Ubuntu   18.04  Moajave   10.x
  -      20.04  Catalina  15.15
=======  =====  ========  =====        

Only gcc is currently support:
gcc (libc 2.27+): 7.x, 8.x, 9.x, 10.x

.. warning::
   It is best to use gcc/7.x + as earlier gcc versions do not evaluate the constexpr hashes at compile time, leading to lower performance.
   Example is here https://www.godbolt.org/z/PHZ4P4

.. warning::
   Below, references to building with the Intel compiler are made. Unfortunately the Intel compiler doesn't currently work with applications that also
   link against GSL. This is being investigated. For now, please do no build CHM with Intel Compilers.

If using conan to build the dependencies, the only requirements are:

   - conan (via Python pip)
   - cmake >=3.17  (via apt-get/brew)
   - C++14 compiler (gcc 7.2+) (via apt-get/brew)
   - Fortran 90+ compiler (gfortran) (via apt-get/brew)
   - m4 (via apt-get/brew)
   - autotools (although this is generally installed in most environments) (via apt-get/brew)
   - BLAS library (via apt-get/brew) e.g., ``libopenblas-dev``

Using system sqlite3, curl, and libtiff enables the compilation of gdal to proceed smoothly.

On Ubuntu 20.04 these can be installed as:

::

   libopenblas-dev
   libtiff-dev
   libcurl4-gnutls-dev
   libsqlite3-dev
   sqlite3

.. note::

   For distributed MPI support, optionally ensure MPI is installed:
   ``libopenmpi-dev`` and ``openmpi-bin``



On MacOS, `homebrew <https://brew.sh/>`__ should be used to install
cmake and optionally conan. Macport based installs likely work, but have not been
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

The required 3rd party libraries needed for CHM can be built using conan.

Throughout, this section assumes a working development environment, but
a blank conan environment. 

Setup conan
-------------

::

   conan profile new default --detect
   conan profile update settings.compiler.cppstd=14 default

conan needs to be told to use new C++11 ABI. If using clang (e.g.,
Macos), do

::

   conan profile update settings.compiler.libcxx=libc++ default  #with clang

and if using gcc, do

::

   conan profile update settings.compiler.libcxx=libstdc++11 default  #with gcc

If you change compilers, such as on a cluster with a modules system, you
can rerun

::

   conan profile new default --detect --force

to detect the new compiler settings. The ``cppstd`` and ``libcxx``
settings need to be reapplied once this is done.

The next few steps can be avoided by installing the following config,

::

   conan config install https://github.com/Chrismarsh/conan-config.git

This adds the CHM remote, enables revisions, and applies a ``settings.yml`` file that allows for distinguishing between the libc
versions of ubuntu 18.04 and 20.04.

Add Conan remote
-----------------

Add the CHM artifactory repository as a conan remote -- this is where the conan scripts to build the dependencies reside.

::

   conan remote add chm http://conan.snowcast.ca/artifactory/api/conan/chm


.. note::

   If the above Conan remote is not working, you can use the ``conan-`` submodules to initialize the local conan build.

   Initialize the submodules that contain the conan recipes

   ::

      cd CHM && git submodule update --init --recursive  # get recipes for dependency builds
      ./conan_export_deps.sh  # tell conan which versions are needed


Enable revisions
-----------------
Enable conan `revisions <https://docs.conan.io/en/latest/versioning/revisions.html#how-to-activate-the-revisions>`__ by
adding ``revisions_enabled=1`` in the ``[general]`` section of your conan.conf file.

Build
--------
This step will install the dependencies into your local conan cache (``~/.conan/data``).
Further, this command will produce the ``FindXXX.cmake`` files required for the
CHM build.

.. note::

   If something goes wrong, you can remove this directory (``~/.conan/data``) or a specific package (``~/.conan/data/package``) to "start fresh".

Without MPI
~~~~~~~~~~~~~~

To build without MPI support:

::

   cd ~/build-CHM
   conan install ~/CHM -if=. --build missing


With MPI support
~~~~~~~~~~~~~~~~~~

If MPI is to be used, then include the following ``-o`` switches:

::

   conan install ~/CHM -if=. -o boost:without_mpi=False -o trilinos:with_mpi=True --build missing

During the CHM cmake configure step, ensure you enable MPI!

Various gotchas
-----------------

Note that custom options can be specified for any of the dependencies using ``-o package:option=value`` at the ``conan install`` stage.

Trilinos
~~~~~~~~~

Trilinos is the only dependency that is not obvious to setup. Because of the tuned nature of BLAS and LAPACK libraries,
only system BLAS and LAPACK are used in compilation.


Intel MKL
~~~~~~~~~

.. warning::
   Using MKL with Trilinos is not supported as the final CHM link will conflict with the internal BLAS in GSL.


OpenBLAS
~~~~~~~~~

Linking Trilinos against OpenBLAS is the best option as it has the LAPACK API.

Set the conan option ```-o trilinos:with_openblas=True`` to change the link library name to ``openblas``.
This may only be useful on some systems. E.g., homebrew openblas has a ``lblas`` symlink.

Custom BLAS location
~~~~~~~~~~~~~~~~~~~~~~

The Trilinos dependencies look for the BLAS libraries in a standard location.
On HPC machines this will almost certainly fail, so the location of the library direction may be set via the env var
``$BLASROOT``. LAPACK search will be set to the same path.

If a custom BLAS location is specified to build Trilinos, this will be automatically detected for the final CHM link.

MacOS
~~~~~~

Homebrew should be used to install -- ``brew install openblas``. A homebrew installed ``openblas`` will be automatically detected and used.
This is prefered over the system default Accelerate framework.


OpenMP
~~~~~~

On MacOS, the openmp library should be installed via homebrew:

::

   brew install libomp


.. warning::
   The Trilinos openmp implementation is not compatible with homebrew omp. It is automatically disabled. It can be explicitly disabled via
   ``-o trilinos:with_openmp=False``



Build CHM
***********

Setup CHM source folders
------------------------

An out of source build should be used. That is, build in a separate folder removed from the CHM source. This makes it easier to clean up
and start from scratch. An example is given below:

::

   cd ~/
   git clone https://github.com/Chrismarsh/CHM

   mkdir ~/build-CHM

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


To use MPI, pass the following to cmake

::

   cmake ~/CHM <other args here> -DUSE_MPI=TRUE


Intel compiler
~~~~~~~~~~~~~~

If the Intel compiler is used, add the following cmake flags:

::

   -DCMAKE_CXX_COMPILER=icpc -DCMAKE_C_COMPILER=icc -DCMAKE_FORTRAN_COMPILER=ifort

High performance allocators
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default tcmalloc is used. Optionally, if system `jemalloc` is available it can be enabled with
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

TCMALLOC
--------

TCmalloc may need to be disabled and can be done via
``-DUSE_TCMALLOC=FALSE``

gepertool heap profiler & libunwnd
----------------------------------

Some machines do not build gperftools with the heap profiling correctly.
This can be disabled when building gperftools

::

   conan install ~/code/CHM/ -if=. --build missing -o gperftools:heapprof=False

Full build including dependencies (summary)
***********************************************

In summary a full MPI Release build of CHM (this assumes conan is setup correctly)

::

   cd ~/
   git clone https://github.com/Chrismarsh/CHM  # get CHM source code
   mkdir ~/build-CHM && cd ~/build-CHM  # create a build directory
   conan install ~/CHM -if=. -o boost:without_mpi=False -o trilinos:with_mpi=True --build missing  # build dependencies that haven't been built, produce custom FindXXX.cmake for all dependencies
   cmake ~/CHM -DUSE_MPI=ON # run cmake configuration
   make -j   # build the CHM executable using all build threads



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
   module load boost-mpi
   module load openblas
   module load gsl
   module load eigen/3.3.7
   module load sparsehash
   module load tbb
   module load trilinos/chm
   module load netcdf-c++4
   module load vtk
   module load proj
   module load jemalloc
   module load cmake

Optionally you can save this with ``module save chm``.


Then build CHM

::

   git clone https://github.com/Chrismarsh/CHM  # get CHM source code
   mkdir ~/chm-build && cd ~/chm-build # make build directory
   cmake ../CHM -DBUILD_WITH_CONAN=FALSE -DUSE_MPI=TRUE -DENABLE_SAFE_CHECKS=ON -DBoost_NO_BOOST_CMAKE=ON -DUSE_TCMALLOC=FALSE -DUSE_JEMALLOC=TRUE -DCMAKE_BUILD_TYPE=Release
   make -j10


