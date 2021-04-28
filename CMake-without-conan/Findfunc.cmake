###############################################################################
# CMake module to search for func library
#
# FUNC_ROOT = install prefix to search

# On success, the macro sets the following variables:
# FUNC_FOUND       = if the library found
# FUNC_LIBRARY     = full path to the library
# FUNC_INCLUDE_DIR = where to find the library headers
# also defined, but not for general use are
# FUNC_LIBRARY, where to find the func library.
#
# Copyright (c) 2009 Mateusz Loskot <mateusz@loskot.net>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#
###############################################################################
SET(FUNC_INCLUDE_SEARCH_PATHS
	${FUNC_DIR}/include
    )

SET(FUNC_FOUND ON)

FIND_PATH(FUNC_INCLUDE_DIR NAMES func.hpp PATHS "${FUNC_INCLUDE_SEARCH_PATHS}/func" DOC "Include for func")

IF(NOT FUNC_INCLUDE_DIR)
    SET(FUNC_FOUND OFF)
ENDIF()

IF (FUNC_FOUND)
#    IF (NOT VIENNACL_FIND_QUIETLY)
        MESSAGE(STATUS "Found func include: ${FUNC_INCLUDE_DIR}")
#    ENDIF (NOT VIENNACL_FIND_QUIETLY)
ELSE (FUNC_FOUND)
    IF (FUNC_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Could not find func")
    ENDIF (FUNC_FIND_REQUIRED)
ENDIF (FUNC_FOUND)

LIST( REMOVE_DUPLICATES FUNC_INCLUDE_DIR )

SET( HAVE_FUNC TRUE )
message(STATUS "Func detected: " ${FUNC_INCLUDE_DIR})

find_library(FUNC_LIBRARY
        NAMES func
        PATHS ${FUNC_DIR}/lib/
        )
find_library(FUNC_IMPLS_LIBRARY
        NAMES func_impls
        PATHS ${FUNC_DIR}/lib/
        )

MARK_AS_ADVANCED(
        FUNC_INCLUDE_DIR
        FUNC_LIBRARY
        FUNC_IMPLS_LIBRARY
)

