from conans import ConanFile, CMake, tools
import os

class CHMConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"

    name = "CHM"
    version = "1.1"
    license = "https://github.com/Chrismarsh/CHM/blob/master/LICENSE"
    author = "Chris Marsh"
    url = "https://github.com/Chrismarsh/CHM"
    description = "Canadian hydrological model"
    generators = "cmake_find_package"

    options = {
        "verbose_cmake":[True,False],
        "build_tests":[True,False],
        "with_mpi": [True, False],
        "with_omp": [True,False]
    }

    default_options = {
       "verbose_cmake":False,
       "build_tests":False,
        #default without openmp or mpi
        "with_omp": True,
        "with_mpi": False,

        #dependency options
        "*:shared": True,
        "gdal:with_curl": True,
        "gdal:with_netcdf": True,

        # "netcdf:dap": False,

        "proj:with_curl": False,

        "hdf5:enable_cxx": True

    }

    def source(self):

        #master default
        branch = None

        try:
            is_ci = os.environ["CI"]
        except KeyError as e:
            raise Exception('This conanfile is intended to be called from within a CI environment. Please use CMake '
                            'to compile.')

        try:
            branch = os.environ["GITHUB_SHA"]
        except KeyError as e:
            try:
                if os.environ["CI"]:
                    self.output.error('When running under CI, $GITHUB_SHA should be available.')
            except KeyError as e:
                pass


        git = tools.Git()
        git.clone("https://github.com/Chrismarsh/CHM.git")
        if branch is None:
            raise Exception('No branch specified')

        git.run(f'checkout {branch}')
        git.run("submodule update --init --recursive")




    def requirements(self):
        self.requires( "cgal/[>=5.2]@CHM/stable" )
        self.requires( "boost/[>=1.75]@CHM/stable" )
        self.requires( "vtk/[>=9.0.1]@CHM/stable" )

        self.requires("netcdf/4.8.1" )
        self.requires("netcdf-cxx/4.3.1@CHM/stable" )

        # Guide libtiff (via gdal) to make the right version selection, but this is not a depednency we explicitly declare
        # self.requires( "libdeflate/1.12", override=True)
        self.requires( "proj/9.0.1" )
        self.requires( "gdal/3.5.2" )

        #because libproj and gdal conflict on what libtiff to use
        self.requires( "zlib/1.2.13", override=True)
        self.requires( "libtiff/4.4.0", override=True)
        self.requires( "sqlite3/3.39.3", override=True)
        self.requires( "libcurl/7.85.0", override=True)

        self.requires( "sparsehash/[>=2.0.3]@CHM/stable" )
        self.requires( "gperftools/[>=2.7]@CHM/stable" )
        self.requires( "gsl/[>=2.6]@CHM/stable" )
        self.requires( "armadillo/[>=10.2.0]@CHM/stable" )
        self.requires( "onetbb/[>=2021.5.0]@CHM/stable" )
        self.requires( "eigen/[>=3.3.9]" )

        self.requires( "meteoio/2.8.0@CHM/stable")
        self.requires( "func/0.1@CHM/stable")
        self.requires( "trilinos/13.4.0@CHM/stable")


    def _configure_cmake(self):
        cmake = CMake(self)

        if self.options['with_mpi']:
        #default to no MPI
            self.options["boost:without_mpi"] = False
            self.options["trilinos:with_mpi"] = True

        if self.options["with_omp"]:
            # trilinos does not support omp on macos
            if not tools.os_info.is_macos:
                self.options["trilinos:with_openmp"] = True


        if self.options.build_tests:
            cmake.definitions["BUILD_TESTS"] = "ON"

        if self.options.verbose_cmake:
            cmake.verbose = True
            cmake.definitions["CMAKE_FIND_DEBUG_MODE"]=1

        if self.options["with_omp"]:
            cmake.definitions["USE_OMP"] = "ON"

        if self.options["with_mpi"]:
            cmake.definitions["USE_MPI"] = "ON"

        cmake.configure(source_folder=self.source_folder)

        return cmake

    def build(self):
        cmake = self._configure_cmake()
        cmake.build()
        # cmake.test(target="check")

    def package(self):
        cmake = self._configure_cmake()
        cmake.install()


    def imports(self):
        self.copy("*.so*", dst="lib", src="lib")  # From bin to bin
        self.copy("*.dylib*", dst="lib", src="lib")  # From lib to bin
