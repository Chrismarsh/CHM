###############################################################################
# CMake module to search for PROJ library
#
# PROJ_ROOT = install prefix to search
#
# On success, the macro sets the following variables:
# PROJ_FOUND       = if the library found
# PROJ_LIBRARY     = full path to the library
# PROJ_INCLUDE_DIR = where to find the library headers
#
# Redistribution and use is allowed according to the terms of the BSD license.
#
###############################################################################

IF( DEFINED ENV{PROJ_DIR} )
    SET( PROJ_DIR "$ENV{PROJ_DIR}" )
ENDIF()

find_path(PROJ_INCLUDE_DIR
        proj.h
        PATHS ${PROJ_DIR}/include
        DOC "Include for PROJ")

find_library(PROJ_LIBRARY
        NAMES proj
        PATHS ${PROJ_DIR}/lib ${PROJ_DIR}/lib64
        )

find_package_handle_standard_args(PROJ DEFAULT_MSG
        PROJ_INCLUDE_DIR PROJ_LIBRARY)

if(PROJ_FOUND)
    set( PROJ_INCLUDE_DIRS ${PROJ_INCLUDE_DIR})
    set( PROJ_LIBRARIES ${PROJ_LIBRARY})

    mark_as_advanced(
            PROJ_INCLUDE_DIR
            PROJ_LIBRARY

    )

    add_library(PROJ::PROJ INTERFACE IMPORTED)

    set_target_properties(PROJ::PROJ PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${PROJ_INCLUDE_DIRS}")

    set_property(TARGET PROJ::PROJ PROPERTY INTERFACE_LINK_LIBRARIES "${PROJ_LIBRARIES}")

endif()

