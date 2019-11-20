from conans import ConanFile, CMake
from pathlib import Path


class CgalTestConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    generators = "cmake"

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def imports(self):
        self.copy("*.dll", dst="bin", src="bin")
        self.copy("*.so*", dst="bin", src="lib")

    def test(self):
        self.run(str(Path("bin").joinpath("example")))
