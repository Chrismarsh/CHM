Compilation
============

CHM uses `conan <https://conan.io/>`__ to manage and build all
dependencies. Because of the various requirements on build
configuration, versions, and interdependencies, using system libraries
it not supported.

All of the CHM dependencies are built on Travis-CI and uploaded to the
`bintray <https://bintray.com/chrismarsh/CHM>`__ repository to serve
prebuilt binaries. This means that *if* the CHM build is done with
supported compilers and operating system (described later), the
dependencies do not need to be built by the end user.

Build requirements
*******************

Linux and MacOS are the only supported environments. The following have been tested

=======  =====  ========  =====
   Linux          MacOS
--------------  ---------------
Ubuntu   18.04  Moajave   10.x
  -      20.04  Catalina  15.15
=======  =====  ========  =====        


gcc (libc 2.27+): 7.x, 8.x

.. warning::
   It is best to use gcc/7.x + as earlier gcc versions do not evaluate the constexpr hashes at compile time, leading to lower performance.
   Example is here https://www.godbolt.org/z/PHZ4P4

If using conan to build the dependencies, the only requirements are:

   - conan >=1.21 (via Python pip)
   - cmake >=3.16  (via apt-get/brew)
   - C++14 compiler (gcc 7.2+) (via apt-get/brew)
   - Fortran 90+ compiler (gfortran) (via apt-get/brew)
   - m4 (via apt-get/brew)
   - autotools (although this is generally installed in most environments) (via apt-get/brew)
   - BLAS library (via apt-get/brew) e.g., `libopenblas-dev`

.. note::
   For distributed MPI support, optionally ensure MPI is installed


On MacOS, `homebrew <https://brew.sh/>`__ should be used to install
cmake and optionally conan. Macport based installs likely work, but have not been
tested.


Compile
********

Throughout, this document assumes a working development environment, but
a blank conan environment. 

Setup conan
-----------

::

   conan profile new default --detect
   conan profile update settings.compiler.cppstd=14 default

conan needs to be told to use new C++11 ABI. If using clang (e.g.,
MacOs), do

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

Intel compiler
~~~~~~~~~~~~~~

If the Intel compiler is being used (this is optional), ensure the Intel compilervars is sourced, e.g.,

::

   source /opt/intel/bin/compilervars.sh intel64

prior to running the conan. Use the above gcc settings for conan.

.. warning::

   Intel compiler is not actively tested. There is a known issue using the MKL at the moment.


Setup CHM source folders
------------------------

An out of source build should be used. That is, build in a seperate folder removed from the CHM source. This makes it easier to clean up
and start from scratch. An example is given below:

::

   cd ~/
   git clone https://github.com/Chrismarsh/CHM
   (cd CHM && git submodule update --init --recursive)

   mkdir ~/build-CHM

.. note::
   The follow instructions assume that they are invoked from within ``~/build-CHM`` (or your equivalent).

Setup dependencies
------------------

Initialize the submodules that contain the conan recepies

::

   cd CHM && git submodule update --init --recursive  # get recipes for dependency builds
   ./conan_export_deps.sh  # tell conan which versions are needed

Install the dependencies into your local conan cache (``~/.conan/data``)

::

   cd ~/build-CHM
   conan install ~/CHM -if=. --build missing

Further, this command will produce the ``FindXXX.cmake`` files required for the
CHM build.

Enabling MPI
-------------

If MPI is to be used:

::

   conan install ~/CHM -if=. -o boost:without_mpi=False -o trilinos:with_mpi=True --build missing


Full build including dependencies (summary)
------------------------------------

In summary:

::

   cd ~/
   git clone https://github.com/Chrismarsh/CHM  # get CHM source code
   cd CHM && git submodule update --init --recursive  # get recipes for dependency builds
   ./conan_export_deps.sh  # tell conan which versions are needed

   mkdir ~/build-CHM && cd ~/build-CHM  # create a build directory
   conan install ~/CHM -if=. --build missing  # build dependencies that haven't been built, produce custom FindXXX.cmake for all dependencies
   cmake ~/CHM  # run cmake configuration
   make -j  # build the CHM executable using all build threads

Additionally, configuration can be setup and built with MPI using:

::

   mkdir ~/build-CHM-mpi && cd ~/build-CHM-mpi
   conan install ~/CHM -if=. -o boost:without_mpi=False -o trilinos:with_mpi=True --build missing
   cmake -DUSE_MPI=ON ~/CHM
   make -j

Note that custom options can be specified for any of the dependencies using `-o package:option=value` at the `conan install` stage.

Trilinos
~~~~~~~~~

Trilinos is the only dependency that is not obvious to setup. Because of the tuned nature of BLAS and LAPACK libraries,
only system BLAS and LAPACK are used in compilation.


Intel MKL
++++++++++

.. warning::
   Using MKL with Trilinos is not supported as the final CHM link will conflict with the internal BLAS in GSL.


OpenBLAS
+++++++++

Linking Trilinos against OpenBLAS is the best option as it has the LAPACK API.

Set the conan option ```-o trilinos:with_openblas=True`` to change the link library name to ``openblas``.
This may only be useful on some systems. E.g., homebrew openblas has a ``lblas`` symlink.

Custom BLAS location
++++++++++++++++++++++

The Trilinos dependencies look for the BLAS libraries in a standard location.
On HPC machines this will almost certainly fail, so the location of the library direction may be set via the env var
``$BLASROOT``. LAPACK search will be set to the same path.

If a custom BLAS location is specified to build Trilinos, this will be automatically detected for the final CHM link.

MacOS
+++++++

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

Run cmake
---------

You can set the install prefix to be anywhere, such as shown in the
example below

::

   cmake ~/CHM -DCMAKE_INSTALL_PREFIX=/opt/chm-install

This should complete without any errors. Both ``ninja`` and ``make``
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

The documentation can be built out of source with:

.. code::

   make docs

or it can be built in source tree with

.. code::

   cd docs
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

Matlab
------

OSX
~~~

-  Create a symbolic link from /usr/bin to the matlab install
-  ``sudo ln -s /Applications/MATLAB_R2013a.app/bin/matlab /usr/bin/matlab``

Linux:
~~~~~~

Usage of the matlab engine requires installing ``csh``


Building on WestGrid
*********************

To build on WestGrid’s Graham machine, all dependencies must be built
from source to ensure the correct optimizations are used. As well, Conan
detects libc versions via compiler version, however on the CentOS 7
system on Graham, the libc is much older than the compiler would
suggest, thus the prebuilt libraries will not link correctly.


Setup Conan
-----------

::

   module load gcc/8.3.0
   module load python/3.7.4
   module load cmake/3.16

   virtualenv ~/conan_env
   source ~/conan_env/bin/activate
   pip install conan
   conan profile new default --detect
   conan remote add bincrafters https://api.bintray.com/conan/bincrafters/public-conan
   conan remote add CHM https://api.bintray.com/conan/chrismarsh/CHM
   conan profile update settings.compiler.cppstd=14 default  
   conan profile update settings.compiler.libcxx=libstdc++11 default  #with gcc

If a different gcc version is used,

::

   conan profile new default --detect --force 
   conan profile update settings.compiler.cppstd=14 default  
   conan profile update settings.compiler.libcxx=libstdc++11 default  #with gcc

Needs to be re-run. Doing so will require a full rebuilt of all
dependencies.

Building CHM
------------

Ensure the environment is correctly setup

::

   module load gcc/8.3.0
   module load python/3.7.4
   module load cmake/3.16
   module load openblas/0.3.6
   module load openmpi/4.0.1
   source ~/conan_env/bin/activate

then build dependencies and CHM

::

   mkdir ~/chm-build && cd ~/chm-build
   conan install ~/CHM -if=.  -o boost:without_mpi=False  -o trilinos:with_mpi=True -o trilinos:with_openblas=True -o trilinos:blas_root=$EBROOTOPENBLAS/lib --build
   cmake ~/CHM -DCMAKE_INSTALL_PREFIX=~/bin/CHM -DUSE_MPI=TRUE
   make -j10 install


