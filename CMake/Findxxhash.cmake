###############################################################################
# CMake module to search for XXHASH library
#
# XXHASH_ROOT = install prefix to search

# On success, the macro sets the following variables:
# XXHASH_FOUND       = if the library found
# XXHASH_LIBRARY     = full path to the library
# XXHASH_INCLUDE_DIR = where to find the library headers
# also defined, but not for general use are
# XXHASH_LIBRARY, where to find the XXHASH library.
#
###############################################################################

FIND_PATH(XXHASH_INCLUDE_DIR xxhash.h
        PATHS ${XXHASH_ROOT}/include
        NO_DEFAULT_PATH
        DOC "Path to XXHASH library include directory")

SET(XXHASH_NAMES ${XXHASH_NAMES} xxhash)

FIND_LIBRARY(XXHASH_LIBRARY
        NAMES ${XXHASH_NAMES}
        PATHS ${XXHASH_ROOT}/lib
        PATHS ${XXHASH_ROOT}/lib64
        NO_DEFAULT_PATH
        DOC "Path to XXHASH library file")

# Handle the QUIETLY and REQUIRED arguments
# if all listed variables are TRUE
INCLUDE(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(XXHASH DEFAULT_MSG XXHASH_LIBRARY XXHASH_INCLUDE_DIR)

IF(XXHASH_FOUND)
    SET(XXHASH_LIBRARIES ${XXHASH_LIBRARY})
ENDIF()