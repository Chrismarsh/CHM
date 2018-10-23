#original source for this Cmake find file is
#https://github.com/naibaf7/libdnn/blob/master/cmake/Modules/FindViennaCL.cmake
SET(VIENNACL_INCLUDE_SEARCH_PATHS
	${VCL_ROOT}/include
    .
    ..
    ../ViennaCL
    ../viennacl-dev
    /usr/include
    /usr/local/include
    /opt/viennacl/include
    /opt/Viennacl/include
    /opt/ViennaCL/include
    )

SET(VIENNACL_FOUND ON)

FIND_PATH(VIENNACL_INCLUDE_DIR NAMES viennacl/version.hpp PATHS ${VIENNACL_INCLUDE_SEARCH_PATHS} DOC "Include for ViennaCL")

IF(NOT VIENNACL_INCLUDE_DIR)
    SET(VIENNACL_FOUND OFF)
ENDIF()

IF (VIENNACL_FOUND)
#    IF (NOT VIENNACL_FIND_QUIETLY)
        MESSAGE(STATUS "Found ViennaCL include: ${VIENNACL_INCLUDE_DIR}")
#    ENDIF (NOT VIENNACL_FIND_QUIETLY)
ELSE (VIENNACL_FOUND)
    IF (VIENNACL_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Could not find VIENNACL")
    ENDIF (VIENNACL_FIND_REQUIRED)
ENDIF (VIENNACL_FOUND)

LIST( REMOVE_DUPLICATES VIENNACL_INCLUDE_DIR )

SET( HAVE_VIENNACL TRUE )
message(STATUS "ViennaCL detected: " ${VIENNACL_INCLUDE_DIR})

MARK_AS_ADVANCED(
        VIENNACL_INCLUDE_DIR
        VIENNACL_LIBRARIES
)
