#from https://github.com/vice87/gam-ngs/blob/master/cmake/Modules/FindSparsehash.cmake

# - Try to find Sparsehash
# Once done this will define
#  SPARSEHASH_FOUND - System has Sparsehash
#  SPARSEHASH_INCLUDE_DIRS - The Sparsehash include directories

IF( DEFINED ENV{SPARSEHASH_DIR} )
    SET( SPARSEHASH_DIR "$ENV{SPARSEHASH_DIR}" )
ENDIF()

find_path(SPARSEHASH_INCLUDE_DIR
        include/google/sparse_hash_map
        HINTS ${SPARSEHASH_DIR}

        DOC "Include for SPARSEHASH")

find_package_handle_standard_args(SPARSEHASH DEFAULT_MSG SPARSEHASH_INCLUDE_DIR)

if(SPARSEHASH_FOUND)
    set(SPARSEHASH_INCLUDE_DIRS ${SPARSEHASH_INCLUDE_DIR})
    MARK_AS_ADVANCED(
            SPARSEHASH_INCLUDE_DIR
            SPARSEHASH_DIR
    )

    message(STATUS "SPARSEHASH found: " ${SPARSEHASH_INCLUDE_DIR})

    add_library(Sparsehash::Sparsehash INTERFACE IMPORTED)
    set_target_properties(Sparsehash::Sparsehash PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${SPARSEHASH_INCLUDE_DIRS})
endif()

