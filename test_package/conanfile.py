from conans import ConanFile, CMake
import os


class CHMTestConan(ConanFile):
    settings = "os", "compiler", "build_type", "arch"
    requires = "CHM/10@CHM/stable"
    generators = "cmake"

    def build(self):
        pass
        # self.run('%s ../../message.proto --proto_path=../.. --cpp_out="."'
        #          % os.path.join('.', 'bin', 'protoc'))
        # cmake = CMake(self.settings)
        # self.run('cmake "%s" %s' % (self.conanfile_directory, cmake.command_line))
        # self.run("cmake --build . %s" % cmake.build_config)
        # if self.settings.os == "Macos":
        #     self.run("cd bin; for LINK_DESTINATION in $(otool -L client | grep libproto | cut -f 1 -d' '); do install_name_tool -change \"$LINK_DESTINATION\" \"@executable_path/$(basename $LINK_DESTINATION)\" client; done")

    def imports(self):
        self.copy("CHM", "bin", "bin") 
        self.copy("*.so*", "lib", "lib") 
        self.copy("*.dylib*", "lib", "lib") 


    def test(self):
        cmd = "%s -v" % os.path.join(".", "bin", "CHM")
        print(cmd)
        self.run(cmd)
        return 0