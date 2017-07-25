# Copyright 2012 by Kitware, Inc. All Rights Reserved. Please refer to
# KITWARE_LICENSE.TXT for licensing information, or contact General Counsel,
# Kitware, Inc., 28 Corporate Drive, Clifton Park, NY 12065.

# Locate the system installed PROJ
#
# The following variables will guide the build:
#
# PROJ4_ROOT        - Set to the install prefix of the PROJ library
#
# The following variables will be set:
#
# PROJ4_FOUND       - Set to true if PROJ can be found
# PROJ4_INCLUDE_DIR - The path to the PROJ header files
# PROJ4_LIBRARY     - The full path to the PROJ library

if( PROJ4_DIR )
    find_package( PROJ NO_MODULE )
elseif( NOT PROJ4_FOUND )

    # Backup the previous root path
    if(PROJ4_ROOT)
	set(_CMAKE_FIND_ROOT_PATH ${CMAKE_FIND_ROOT_PATH})
	set(CMAKE_FIND_ROOT_PATH ${PROJ4_ROOT})
	set(_PROJ4_ROOT_OPTS ONLY_CMAKE_FIND_ROOT_PATH)
    endif()

    find_path( PROJ4_INCLUDE_DIR proj_api.h ${_PROJ4_ROOT_OPTS})
    find_library( PROJ4_LIBRARY proj ${_PROJ4_ROOT_OPTS})

    # Restore the original root path
    if(PROJ4_ROOT)
	set(CMAKE_FIND_ROOT_PATH ${CMAKE_FIND_ROOT_PATH})
    endif()

    include( FindPackageHandleStandardArgs )
    FIND_PACKAGE_HANDLE_STANDARD_ARGS( PROJ4 PROJ4_INCLUDE_DIR PROJ4_LIBRARY )
endif()
