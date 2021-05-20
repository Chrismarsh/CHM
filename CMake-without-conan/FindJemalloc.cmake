# - Find Jemalloc
# Find the native Jemalloc includes and library
#
#  Jemalloc_INCLUDE_DIR - where to find incl, etc.
#  Jemalloc_LIBRARIES   - List of libraries when using Jemalloc.
#  Jemalloc_FOUND       - True if Jemalloc found.


IF( DEFINED ENV{Jemalloc_DIR} )
    SET( Jemalloc_DIR "$ENV{Jemalloc_DIR}" )
ENDIF()


find_path(Jemalloc_INCLUDE_DIR
        jemalloc/jemalloc.h
        PATHS ${Jemalloc_DIR}/include
        )

find_library(Jemalloc_LIBRARY
        NAMES jemalloc libjemalloc JEMALLOC
        PATHS ${Jemalloc_DIR}/lib
        )

find_package_handle_standard_args(Jemalloc DEFAULT_MSG
        Jemalloc_INCLUDE_DIR Jemalloc_LIBRARY)


if(Jemalloc_FOUND)

    set(Jemalloc_INCLUDE_DIRS ${Jemalloc_INCLUDE_DIR} )
    set(Jemalloc_LIBRARIES ${Jemalloc_LIBRARY} )

    mark_as_advanced(
            Jemalloc_LIBRARY
            Jemalloc_INCLUDE_DIR
    )

    add_library(Jemalloc::Jemalloc INTERFACE IMPORTED)

    set_target_properties(Jemalloc::Jemalloc PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${Jemalloc_INCLUDE_DIRS}")
    set_property(TARGET Jemalloc::Jemalloc PROPERTY INTERFACE_LINK_LIBRARIES "${Jemalloc_LIBRARIES}")

endif()