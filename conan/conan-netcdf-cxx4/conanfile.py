from conans import ConanFile, CMake, tools
from conans.errors import ConanInvalidConfiguration

class NetcdfcConan(ConanFile):
    name = "netcdf-cxx"
    version = "4.3.1"
    license = "MIT"
    author = "Lars Bilke, lars.bilke@ufz.de"
    url = "https://github.com/bilke/conan-netcdf-cxx"
    description = "Unidata network Common Data Form cxx"
    settings = "os", "compiler", "build_type", "arch"
    options = {"shared": [True, False], "fPIC": [True, False]}
    default_options = "shared=True", "fPIC=True"
    generators = "cmake"

    def source(self):
        self.run("git clone --depth=1 https://github.com/Unidata/netcdf-cxx4.git")

        if tools.os_info.is_macos:
            tools.replace_in_file("netcdf-cxx4/CMakeLists.txt", "PROJECT(NCXX C CXX)",
                                  '''PROJECT(NCXX C CXX)
                                    include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
                                    conan_basic_setup(KEEP_RPATHS)''')
        else:
            tools.replace_in_file("netcdf-cxx4/CMakeLists.txt", "PROJECT(NCXX C CXX)",
                                  '''PROJECT(NCXX C CXX)
                                    include(${CMAKE_BINARY_DIR}/conanbuildinfo.cmake)
                                    conan_basic_setup()''')

    def requirements(self):
        self.requires("netcdf-c/4.6.2@CHM/stable")

    # def config_options(self):
        # if self.settings.os == "Windows":
            # del self.options.fPIC

    def configure(self):
        del self.settings.compiler.libcxx
        del self.settings.compiler.cppstd
        # if self.settings.os == "Windows" and self.options.shared:
            # raise ConanInvalidConfiguration("Windows shared builds are not supported right now")

    def configure_cmake(self):
        cmake = CMake(self)
        cmake.definitions["NCXX_ENABLE_TESTS"] = False
        cmake.definitions["ENABLE_CONVERSION_WARNINGS"] = False
        cmake.definitions["BUILD_SHARED_LIBS"] = self.options.shared

        if tools.os_info.is_macos:
            cmake.definitions["CMAKE_INSTALL_NAME_DIR"] = "@rpath"
        cmake.configure(source_folder="netcdf-cxx4")
        return cmake

    def build(self):
        cmake = self.configure_cmake()
        cmake.build()

    def package(self):
        cmake = self.configure_cmake()
        cmake.install()

    def package_info(self):
        self.cpp_info.libs = ["netcdf-cxx4"]
