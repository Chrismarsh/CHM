from cpt.packager import ConanMultiPackager
from collections import defaultdict
from sys import platform

if __name__ == "__main__":

    if platform == "linux":
        command = "sudo apt-get -qq update && sudo apt-get -qq install -y patchelf && sudo apt-get -qq install -y gfortran"

    builder = ConanMultiPackager(cppstds=[14],
                                archs=["x86_64"],
                                build_types=["Release"],
                                docker_entry_script = command)
                              
    builder.add_common_builds(pure_c=False)

    builder.remove_build_if(lambda build: build.settings["compiler.libcxx"] == "libstdc++")

    named_builds = defaultdict(list)
    for settings, options, env_vars, build_requires, reference in builder.items:

        shared="shared"

        named_builds[settings['compiler'] +"_"+shared].append([settings, options, env_vars, build_requires, reference])

    builder.named_builds = named_builds

    builder.run()

