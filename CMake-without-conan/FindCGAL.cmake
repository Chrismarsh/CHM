###############################################################################
# CMake module to search for CGAL library
#
# CGAL_ROOT = install prefix to search

# On success, the macro sets the following variables:
# CGAL_FOUND       = if the library found
# CGAL_LIBRARY     = full path to the library
# CGAL_INCLUDE_DIR = where to find the library headers
# also defined, but not for general use are
# CGAL_LIBRARY, where to find the CGAL library.
#
# Copyright (c) 2009 Mateusz Loskot <mateusz@loskot.net>
#
# Redistribution and use is allowed according to the terms of the BSD license.
# For details see the accompanying COPYING-CMAKE-SCRIPTS file.
#
###############################################################################

find_package_handle_standard_args(CGAL CGAL_INCLUDE_DIR)

set(CGAL_INCLUDE_SEARCH_PATHS
	${CGAL_DIR}/include
    )


find_path(CGAL_INCLUDE_DIR NAMES version.h PATHS "${CGAL_INCLUDE_SEARCH_PATHS}/CGAL" DOC "Include for CGAL")

list( REMOVE_DUPLICATES CGAL_INCLUDE_DIR )

message(STATUS "CGAL detected: " ${CGAL_INCLUDE_DIR})


MARK_AS_ADVANCED(
        CGAL_INCLUDE_DIR
)

