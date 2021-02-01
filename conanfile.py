from conans import ConanFile, CMake, tools
import os

class CHMConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"

    name = "CHM"
    version = "1.0"
    license = "https://github.com/Chrismarsh/CHM/blob/master/LICENSE"
    author = "Chris Marsh"
    url = "https://github.com/Chrismarsh/CHM"
    description = "Canadian hydrological model"
    generators = "cmake_find_package"

    options = {
        "verbose_cmake":[True,False],
        "build_tests":[True,False],
        "with_mpi":[True, False],
        "with_omp": [True,False]
    }

    default_options = {
       "verbose_cmake":False,
       "build_tests":True,
        #default without openmp or mpi
        "with_omp": False,
        "with_mpi": False,


        #dependency options
        "gdal:libcurl": True,
        "gdal:netcdf": True
        # "gperftools:heapprof":True
    }

    def source(self):
        try:
            branch = os.environ["GITHUB_REF"]
        except KeyError as e:
            branch = "github-actions"

        git = tools.Git()
        git.clone("https://github.com/Chrismarsh/CHM.git",branch=branch)
        git.run("submodule update --init --recursive")

        if self.options['with_mpi']:
            #default to no MPI
            self.options["boost:without_mpi"] = False
            self.options["trilinos:with_mpi"] = True

        if self.options["with_omp"]:
            self.options["trilinos:with_openmp"] = False

    def requirements(self):

        self.requires( "cgal/[>=5.2]@CHM/stable" )
        self.requires( "boost/[>=1.71]@CHM/stable" )
        self.requires( "vtk/8.2.0@CHM/stable" )
        self.requires( "netcdf-cxx/[>=4.3]@CHM/stable" )
        self.requires( "proj/[>=7.2.1]@CHM/stable" )
        self.requires( "gdal/[>=3.2.1]@CHM/stable" )
        self.requires( "sparsehash/[>=2.0.3]@CHM/stable" )
        self.requires( "gperftools/[>=2.7]@CHM/stable" )
        self.requires( "gsl/[>=2.6]@CHM/stable" )
        self.requires( "armadillo/[>=10.2.0]@CHM/stable" )
        self.requires( "tbb/[>=2020.3]@CHM/stable" )
        self.requires( "eigen3/[>=3.3.9]@CHM/stable" )
        self.requires( "meteoio/2.8.0@CHM/stable")
        self.requires( "func/0.1@CHM/stable")
        self.requires( "trilinos/12.18.1@CHM/stable")

    def _configure_cmake(self):
        cmake = CMake(self)

        if self.options.build_tests:
            cmake.definitions["BUILD_TESTS"] = True

        if self.options.verbose_cmake:
            cmake.verbose = True
            cmake.definitions["CMAKE_FIND_DEBUG_MODE"]=1

        if self.options["with_omp"]:
            cmake.definitions["USE_OMP"] = True

        if self.options["with_mpi"]:
            cmake.definitions["USE_MPI"] = True

        cmake.configure(source_folder=self.source_folder)

        return cmake

    def build(self):
        cmake = self._configure_cmake()
        cmake.build()
        cmake.test(target="check")

    def package(self):
        cmake = self._configure_cmake()
        cmake.install()


    def imports(self):
        self.copy("*.so*", dst="lib", src="lib")  # From bin to bin
        self.copy("*.dylib*", dst="lib", src="lib")  # From lib to bin
