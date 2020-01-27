from conans import ConanFile, CMake
import os


class CHMTestConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    requires = "CHM/10@CHM/stable"
    generators = "cmake"

    def build(self):
        pass

    def imports(self):
        self.copy("CHM", "bin", "bin") 
        self.copy("*.so*", "lib", "lib") 
        self.copy("*.dylib*", "lib", "lib") 


    def test(self):
        cmd = "%s -v" % os.path.join(".", "bin", "CHM")
        print(cmd)
        self.run(cmd)
        return 0