# - Find Gperftools
# Find the native Gperftools includes and library
#
#  Gperftools_INCLUDE_DIR - where to find Gperftools.h, etc.
#  Gperftools_LIBRARIES   - List of libraries when using Gperftools.
#  Gperftools_FOUND       - True if Gperftools found.

include (FindPackageHandleStandardArgs)

IF( DEFINED ENV{GperftoolsF_DIR} )
    SET( GperftoolsF_DIR "$ENV{GperftoolsF_DIR}" )
ENDIF()


find_path(Gperftools_INCLUDE_DIR
        google/tcmalloc.h
        PATHS ${GperftoolsF_DIR}/include
        )

find_library(Gperftools_LIBRARY
        NAMES tcmalloc_minimal
        PATHS ${GperftoolsF_DIR}/lib
        )

find_package_handle_standard_args(Gperftools DEFAULT_MSG
        Gperftools_INCLUDE_DIR Gperftools_LIBRARY)


if(Gperftools_FOUND)

    set(Gperftools_INCLUDE_DIRS ${Gperftools_INCLUDE_DIR} )
    set(Gperftools_LIBRARIES ${Gperftools_LIBRARY} )

    mark_as_advanced(
            Gperftools_LIBRARY
            Gperftools_INCLUDE_DIR
    )

    add_library(Gperftools::Gperftools INTERFACE IMPORTED)

    set_target_properties(Gperftools::Gperftools PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${Gperftools_INCLUDE_DIRS}")
    set_property(TARGET Gperftools::Gperftools PROPERTY INTERFACE_LINK_LIBRARIES "${Gperftools_LIBRARIES}")

endif()