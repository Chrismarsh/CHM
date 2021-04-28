# - Find Gperftools
# Find the native Gperftools includes and library
#
#  Gperftools_INCLUDE_DIR - where to find Gperftools.h, etc.
#  Gperftools_LIBRARIES   - List of libraries when using Gperftools.
#  Gperftools_FOUND       - True if Gperftools found.

find_package_handle_standard_args(Gperftools Gperftools_INCLUDE_DIR Gperftools_LIBRARY)

SET(Gperftools_INCLUDE_SEARCH_PATHS
        ${Gperftools_DIR}/include
        )

find_path(Gperftools_INCLUDE_DIR NAMES tcmalloc.h PATHS "${Gperftools_INCLUDE_SEARCH_PATHS}/google/tcmalloc.h")



#if (USE_Gperftools)
#    set(Gperftools_NAMES Gperftools)
#else ()
set(Gperftools_NAMES tcmalloc_minimal)
#endif ()

message(STATUS "${Gperftools_INCLUDE_DIR}/../")
find_library(Gperftools_LIBRARY
        NAMES ${Gperftools_NAMES}
        PATHS ${Gperftools_INCLUDE_DIR}/../lib/
        )


mark_as_advanced(
        Gperftools_LIBRARY
        Gperftools_INCLUDE_DIR
)