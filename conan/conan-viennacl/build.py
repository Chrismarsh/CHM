from conan.packager import ConanMultiPackager

if __name__ == "__main__":
    builder = ConanMultiPackager(username="bilke")
    builder.add_common_builds()
    builder.run()
