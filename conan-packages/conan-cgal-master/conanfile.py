from conans import ConanFile, CMake, tools


class CgalConan(ConanFile):
    name = "cgal"
    version = "5.0.0"
    license = "GPL/LGPL"
    url = "git@github.com:garlyon/conan-cgal.git"
    description = "Computational Geometry Algorithms Library"
    no_copy_source = True
    settings = "os", "compiler", "build_type", "arch"
    options = {
                "shared": [True, False], 
                "with_gmp": [True, False],
                "with_qt5" : [True, False],
                "with_imageio":[True,False],
                "header_only":[True,False],
                "with_tbb":[True,False]
                }                
    default_options = "shared=False", "with_gmp=False", "with_qt5=False", "with_imageio=False","header_only=True","with_tbb=True"

    scm = {
        "type": "git",
        "url": "https://github.com/cgal/cgal.git",
        "revision": "releases/CGAL-5.0"
    }

    def requirements(self):
        self.requires("boost/[>=1.67]")
        if self.options.with_gmp:
            self.requires("gmp/[>=5.0]@grif/dev")
            self.requires("mpfr/[>=3.0]@grif/dev")
        if self.options.with_tbb:
            self.requires("tbb/2019_u9")

    def build(self):
        with tools.environment_append(self.cmake_env_vars):
            cmake = CMake(self)
            cmake.definitions["BOOST_ROOT"] = self.deps_cpp_info["boost"].rootpath
            cmake.definitions["CGAL_DISABLE_GMP"] = "OFF" if self.options.with_gmp else "ON"
            cmake.definitions["WITH_CGAL_Qt5"] = "OFF" if self.options.with_qt5 else "ON"
            cmake.definitions["WITH_CGAL_ImageIO"] = "OFF" if self.options.with_imageio else "ON"
            cmake.definitions["CGAL_HEADER_ONLY"] = "ON" if self.options.header_only else "OFF"
            cmake.definitions["CMAKE_POSITION_INDEPENDENT_CODE"] = "True"
            cmake.configure()
            cmake.build()
            cmake.install()

    def package(self):
        pass

    def package_info(self):
        self.cpp_info.libs = tools.collect_libs(self)

    @property
    def cmake_env_vars(self):
        env = {}
        if self.options.with_gmp:
            env["GMP_DIR"] = self.deps_cpp_info["gmp"].rootpath
            env["MPFR_DIR"] = self.deps_cpp_info["mpfr"].rootpath
        return env


