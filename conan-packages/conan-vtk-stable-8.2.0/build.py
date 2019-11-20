from conan.packager import ConanMultiPackager
import copy
import platform

if __name__ == "__main__":
    builder = ConanMultiPackager(archs = ["x86_64"])
    builder.add_common_builds(pure_c=False)
    items = []
    for item in builder.items:
        if item.settings["compiler"] == "Visual Studio":
            if item.settings["compiler.runtime"] == "MT" or item.settings["compiler.runtime"] == "MTd":
                # Ignore MT runtime
                continue
        # Build static only
        if item.options["vtk:shared"]:
            continue

        new_options = copy.copy(item.options)
        new_options["vtk:qt"] = True
        new_options["vtk:ioxml"] = True
        items.append([item.settings, new_options, item.env_vars, item.build_requires])

        new_options = copy.copy(item.options)
        new_options["vtk:minimal"] = True
        new_options["vtk:ioxml"] = True
        items.append([item.settings, new_options, item.env_vars, item.build_requires])

    builder.items = items
    builder.run()
