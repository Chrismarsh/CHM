from cpt.packager import ConanMultiPackager
from collections import defaultdict
import os


if __name__ == "__main__":
    builder = ConanMultiPackager(cppstds=[14],
                                archs=["x86_64"],
                                build_types=["Release"],
                                build_policy='missing')
                        
    builder.add_common_builds(pure_c=False, shared_option_name=False)
    builder.remove_build_if(lambda build: build.settings["compiler.libcxx"] == "libstdc++") # old ABI

    builder.update_build_if(lambda build: os.environ['USE_MPI'] == 'with-mpi',
                        new_options={'CHM:with_mpi': False})

    builder.run()
