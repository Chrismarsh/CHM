from conans import ConanFile, CMake
from six import StringIO
import os


class CHMTestConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    requires = "CHM/1.1@CHM/testing"
    generators = "cmake"

    def build(self):
        pass

    def imports(self):
        self.copy("CHM", "bin", "bin") 
        self.copy("*.so*", "lib", "lib") 
        self.copy("*.dylib*", "lib", "lib") 


    def test(self):

        # self.run("""patchelf --set-rpath '\$ORIGIN/../lib' ./bin/CHM""")

        cmd = "%s -v" % os.path.join(".", "bin", "CHM")
        print(cmd)
        self.run(cmd)
        return 0