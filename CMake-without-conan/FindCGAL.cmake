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

IF( DEFINED ENV{CGAL_DIR} )
	SET( CGAL_DIR "$ENV{CGAL_DIR}" )
ENDIF()

find_path(CGAL_INCLUDE_DIR
		 include/CGAL/version.h
		 HINTS ${CGAL_DIR}

		DOC "Include for CGAL")

find_package_handle_standard_args(CGAL DEFAULT_MSG CGAL_INCLUDE_DIR)

if(CGAL_FOUND)
	set(CGAL_INCLUDE_DIRS ${CGAL_INCLUDE_DIR})
	MARK_AS_ADVANCED(
			CGAL_INCLUDE_DIR
			CGAL_DIR
	)

	message(STATUS "CGAL found: " ${CGAL_INCLUDE_DIR})

	add_library(CGAL::CGAL INTERFACE IMPORTED)
	set_target_properties(CGAL::CGAL PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${CGAL_INCLUDE_DIRS})
endif()