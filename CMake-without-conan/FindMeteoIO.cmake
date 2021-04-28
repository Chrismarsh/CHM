###############################################################################
# CMake module to search for METEOIO library
#
# MeteoIO_ROOT = install prefix to search
#
# On success, the macro sets the following variables:
# MeteoIO_FOUND       = if the library found
# MeteoIO_LIBRARY     = full path to the library
# MeteoIO_INCLUDE_DIR = where to find the library headers
#
# Redistribution and use is allowed according to the terms of the BSD license.
#
###############################################################################

IF( DEFINED ENV{MeteoIO_DIR} )
    SET( MeteoIO_DIR "$ENV{MeteoIO_DIR}" )
ENDIF()

find_path(MeteoIO_INCLUDE_DIR
        include/meteoio/MeteoIO.h
        HINTS ${MeteoIO_DIR}
        DOC "Include for meteoio")

find_library(MeteoIO_LIBRARY
        NAMES meteoio
        HINTS ${MeteoIO_DIR}
        )

find_package_handle_standard_args(MeteoIO DEFAULT_MSG
        MeteoIO_INCLUDE_DIR MeteoIO_LIBRARY)

message(STATUS "METEOIO detected: " ${MeteoIO_INCLUDE_DIR})

if(MeteoIO_FOUND)
    set( MeteoIO_INCLUDE_DIRS ${MeteoIO_INCLUDE_DIR})
    set( MeteoIO_LIBRARIES ${MeteoIO_LIBRARY})

    mark_as_advanced(
            MeteoIO_INCLUDE_DIR
            MeteoIO_LIBRARY

    )

    add_library(MeteoIO::MeteoIO INTERFACE IMPORTED)

    set_target_properties(MeteoIO::MeteoIO PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${MeteoIO_INCLUDE_DIRS}")

    set_property(TARGET MeteoIO::MeteoIO PROPERTY INTERFACE_LINK_LIBRARIES "${MeteoIO_LIBRARIES}")

endif()

