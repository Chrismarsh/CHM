from cpt.packager import ConanMultiPackager
from collections import defaultdict

if __name__ == "__main__":
    builder = ConanMultiPackager(cppstds=[14],
                                archs=["x86_64"],
                                build_types=["Release"])
                              
    builder.add_common_builds(pure_c=False,shared_option_name="netcdf-cxx:shared")

    builder.remove_build_if(lambda build: build.settings["compiler.libcxx"] == "libstdc++")

    named_builds = defaultdict(list)
    for settings, options, env_vars, build_requires, reference in builder.items:

        shared="shared"

        if not options['netcdf-cxx:shared']:
            shared = "static" 

        named_builds[settings['compiler'] +"_"+shared].append([settings, options, env_vars, build_requires, reference])

    builder.named_builds = named_builds

    builder.run()

