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
include (FindPackageHandleStandardArgs)

if( DEFINED ENV{Func_DIR} )
    set( Func_DIR "$ENV{Func_DIR}" )
endif()

set(FUNC_FOUND ON)

find_path(FUNC_INCLUDE_DIR
        func/func.hpp
        PATHS ${Func_DIR}/include
        DOC "Include for func"
        )
find_library(FUNC_LIBRARY
        NAMES func
        PATHS ${Func_DIR}/lib
        )
find_library(FUNC_IMPLS_LIBRARY
        NAMES func_impls
        PATHS ${Func_DIR}/lib
        )

find_package_handle_standard_args(Func DEFAULT_MSG
        FUNC_INCLUDE_DIR FUNC_LIBRARY FUNC_IMPLS_LIBRARY)

if(Func_FOUND)
    set( Func_INCLUDE_DIRS ${FUNC_INCLUDE_DIR})
    set( Func_LIBRARIES ${FUNC_LIBRARY} ${FUNC_IMPLS_LIBRARY})

    mark_as_advanced(
            FUNC_INCLUDE_DIR
            FUNC_LIBRARY
            FUNC_IMPLS_LIBRARY
    )

    add_library(Func::Func INTERFACE IMPORTED)
    set_target_properties(Func::Func PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${Func_INCLUDE_DIRS}")
    set_property(TARGET Func::Func PROPERTY INTERFACE_LINK_LIBRARIES "${Func_LIBRARIES}")
endif()



