#from https://github.com/vice87/gam-ngs/blob/master/cmake/Modules/FindSparsehash.cmake

# - Try to find Sparsehash
# Once done this will define
#  SPARSEHASH_FOUND - System has Sparsehash
#  SPARSEHASH_INCLUDE_DIRS - The Sparsehash include directories

find_path(SPARSEHASH_INCLUDE_DIR NAMES google/sparse_hash_map
        PATHS         ${SPARSEHASH_ROOT} $ENV{SPARSEHASH_ROOT} /usr/local/ /usr/ /sw/ /opt/local /opt/csw/ /opt/ ENV CPLUS_INCLUDE_PATH
        PATH_SUFFIXES include/
        )

set(SPARSEHASH_INCLUDE_DIRS ${SPARSEHASH_INCLUDE_DIR})

# handle the QUIETLY and REQUIRED arguments and set SPARSEHASH_FOUND to TRUE if
# all listed variables are TRUE
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(Sparsehash DEFAULT_MSG SPARSEHASH_INCLUDE_DIR)

mark_as_advanced(SPARSEHASH_INCLUDE_DIR)