Conan build script for the GNU Multiple Precision Arithmetic library (GMP)

Build
--

`conan create mpfr/4.0.1@ntc/stable -pr <profile>`

Notes
--

- This was created in order to build CGAL
- This should build directly with g++
- For MSVC, this is cross-compiled using the `msys2_mingw` profile (see
  `profiles/msys2_mingw`, or better yet, the [Windows
  Subsystem](http://docs.conan.io/en/latest/systems_cross_building/windows_subsystems.html?highlight=msys2_mingw)
  Conan docs.)  The `package_id()` function then makes the package appear as if
  it was build with Microsoft Visual Studio.

See also:
- [conan-gmp](https://github.com/kheaactua/conan-gmp).
- [conan-build-helpers](https://github.com/kheaactua/conan-build-helpers).
