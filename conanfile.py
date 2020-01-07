from conans import ConanFile, CMake


class CHMConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    requires = "boost/1.71.0@CHM/stable",
    "cgal/5.0.0@CHM/stable",
    "vtk/8.2.0@CHM/stable",
    "netcdf-cxx/4.3.1@CHM/stable",
    "proj/4.9.3@CHM/stable",
    "gdal/2.4.1@CHM/stable",
    "sparsehash/2.0.3@CHM/stable",
    "gperftools/2.7@CHM/stable",
    "gsl/2.4@CHM/stable",
    "armadillo/9.800.2@CHM/stable",
    "viennacl/1.7.1@CHM/stable",
    "tbb/2019_u9@CHM/stable",
    "eigen3/3.3.7@CHM/stable",
    "meteoio/2.8.0@CHM/stable"


generators = "cmake_find_package"
default_options = {"boost:without_python": True,
                   "boost:without_mpi": True}


# [options]
# boost:without_python=True
# boost:without_mpi=False

# cgal:with_tbb=True
# cgal:with_gmp=True

# netcdf-c:parallel4=False

def imports(self):
    self.copy("*.dll", dst="bin", src="bin")  # From bin to bin
    self.copy("*.dylib*", dst="bin", src="lib")  # From lib to bin
