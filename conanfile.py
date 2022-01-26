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
        "gdal:libcurl": True,
        "gdal:netcdf": True
        # "gperftools:heapprof":True
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
        self.requires( "netcdf-cxx/[>=4.3]@CHM/stable" )
        self.requires( "proj/[>=7.2.1]@CHM/stable" )
        self.requires( "gdal/[>=3.2.1]@CHM/stable" )
        self.requires( "sparsehash/[>=2.0.3]@CHM/stable" )
        self.requires( "gperftools/[>=2.7]@CHM/stable" )
        self.requires( "gsl/[>=2.6]@CHM/stable" )
        self.requires( "armadillo/[>=10.2.0]@CHM/stable" )
        self.requires( "onetbb/[>=2021.3.0]" )
        self.requires( "eigen3/[>=3.3.9]@CHM/stable" )
        self.requires( "meteoio/2.8.0@CHM/stable")
        self.requires( "func/0.1@CHM/stable")
        self.requires( "trilinos/chm@CHM/stable")

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
