#from https://github.com/vice87/gam-ngs/blob/master/cmake/Modules/FindSparsehash.cmake

# - Try to find Sparsehash
# Once done this will define
#  Sparsehash_FOUND - System has Sparsehash
#  Sparsehash_INCLUDE_DIRS - The Sparsehash include directories

IF( DEFINED ENV{Sparsehash_DIR} )
    SET( Sparsehash_DIR "$ENV{Sparsehash_DIR}" )
ENDIF()

find_path(Sparsehash_INCLUDE_DIR
        google/sparse_hash_map
        PATHS ${Sparsehash_DIR}/include

        DOC "Include for Sparsehash")

find_package_handle_standard_args(Sparsehash DEFAULT_MSG Sparsehash_INCLUDE_DIR)

if(Sparsehash_FOUND)
    set(Sparsehash_INCLUDE_DIRS ${Sparsehash_INCLUDE_DIR})
    MARK_AS_ADVANCED(
            Sparsehash_INCLUDE_DIR
            Sparsehash_DIR
    )

    message(STATUS "Sparsehash found: " ${Sparsehash_INCLUDE_DIR})

    add_library(Sparsehash::Sparsehash INTERFACE IMPORTED)
    set_target_properties(Sparsehash::Sparsehash PROPERTIES INTERFACE_INCLUDE_DIRECTORIES ${Sparsehash_INCLUDE_DIRS})
endif()

