# conan-cgal

## Create package localy

Without GMP/MPFR

`conan create . grif/dev`

With GMP/MPFR

`conan create . grif/dev -o cgal:with_gmp=True`

Note: only shared gmp/mpfr working on Windows so far

`conan create . grif/dev -o cgal:with_gmp=True -o gmp:shared=True -o mpfr:shared=True`