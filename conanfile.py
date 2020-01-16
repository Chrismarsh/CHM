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
    # default_options = {"boost:without_python": True,
    #                    "boost:without_mpi": True}

    default_options = {"gperftools:heapprof":True}
    # [options]
    # boost:without_python=True
    # boost:without_mpi=False

    # cgal:with_tbb=True
    # cgal:with_gmp=True

    # netcdf-c:parallel4=False

    def source(self):

        # branch = os.environ.get("TRAVIS_BRANCH","master")
        branch = os.environ["CONAN_TEST_BRANCH"]
        git = tools.Git()
        git.clone("https://github.com/Chrismarsh/CHM.git",branch=branch)
        git.run("submodule update --init --recursive")

        # git.run("clone https://github.com/Chrismarsh/CHM.git")
        # git.run("-C CHM checkout %s" %branch)
        # git.run("-C CHM submodule update --init --recursive")


    def requirements(self):
        self.requires( "cgal/5.0.0@CHM/stable" )
        self.requires( "boost/1.71.0@CHM/stable" )
        self.requires( "vtk/8.2.0@CHM/stable" )
        self.requires( "netcdf-cxx/4.3.1@CHM/stable" )
        self.requires( "proj/4.9.3@CHM/stable" )
        self.requires( "gdal/2.4.1@CHM/stable" )
        self.requires( "sparsehash/2.0.3@CHM/stable" )
        self.requires( "gperftools/2.7@CHM/stable" )
        self.requires( "gsl/2.6@CHM/stable" )
        self.requires( "armadillo/9.800.2@CHM/stable" )
        self.requires( "viennacl/1.7.1@CHM/stable" )
        self.requires( "tbb/2019_u9@CHM/stable" )
        self.requires( "eigen3/3.3.7@CHM/stable" )
        self.requires( "meteoio/2.8.0@CHM/stable")
        
    def build(self):
        cmake = CMake(self)
        cmake.definitions["BUILD_TESTS"] = True
        cmake.configure(source_folder=self.source_folder)
        cmake.build()
        cmake.install()
        self.run('install/CHM -v')

        #cmake.test(target="check")
    # def test(self):
        

    def imports(self):
        self.copy("*.so*", dst="lib", src="lib")  # From bin to bin
        self.copy("*.dylib*", dst="lib", src="lib")  # From lib to bin
