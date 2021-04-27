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

find_package_handle_standard_args(MeteoIO MeteoIO_INCLUDE_DIR MeteoIO_LIBRARY)

SET(MeteoIO_INCLUDE_SEARCH_PATHS
        ${MeteoIO_DIR}/include
        )


find_path(MeteoIO_INCLUDE_DIR NAMES MeteoIO.h PATHS "${MeteoIO_INCLUDE_SEARCH_PATHS}/meteoio" DOC "Include for meteoio")

list(REMOVE_DUPLICATES MeteoIO_INCLUDE_DIR )

set( HAVE_METEOIO TRUE )
message(STATUS "METEOIO detected: " ${MeteoIO_INCLUDE_DIR})

find_library(MeteoIO_LIBRARY
        NAMES meteoio
        PATHS ${MeteoIO_DIR}/lib/
        )


add_library(MeteoIO::MeteoIO INTERFACE IMPORTED)

set_target_properties(MeteoIO::MeteoIO PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${MeteoIO_INCLUDE_DIR}")

set_property(TARGET MeteoIO::MeteoIO PROPERTY INTERFACE_LINK_LIBRARIES "${MeteoIO_LIBRARY}")


mark_as_advanced(
        MeteoIO_INCLUDE_DIR
        MeteoIO_LIBRARY

)
