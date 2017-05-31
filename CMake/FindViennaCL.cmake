SET(VIENNACL_WITH_OPENCL FALSE)

#original source for this Cmake find file is
#https://github.com/naibaf7/libdnn/blob/master/cmake/Modules/FindViennaCL.cmake
SET(VIENNACL_INCLUDE_SEARCH_PATHS
        .
        ..
        ../ViennaCL
        ../viennacl-dev
        /usr/include
        /usr/local/include
        /opt/viennacl/include
        /opt/Viennacl/include
        /opt/ViennaCL/include
        $ENV{VIENNACL_HOME}
        )

SET(VIENNACL_FOUND OFF)

FIND_PATH(VIENNACL_INCLUDE_DIR NAMES viennacl/version.hpp viennacl/linalg/cg.hpp PATHS ${VIENNACL_INCLUDE_SEARCH_PATHS} DOC "Include for ViennaCL")

SET(VIENNACL_FOUND ON)

#    Check include files
IF(NOT VIENNACL_INCLUDE_DIR)
    MESSAGE(STATUS "Could not find VIENNACL include. Turning VIENNACL_FOUND off")
ENDIF()

IF (VIENNACL_FOUND)
    IF (NOT VIENNACL_FIND_QUIETLY)
        MESSAGE(STATUS "Found ViennaCL include: ${VIENNACL_INCLUDE_DIR}")
    ENDIF (NOT VIENNACL_FIND_QUIETLY)
ELSE (VIENNACL_FOUND)
    IF (VIENNACL_FIND_REQUIRED)
        MESSAGE(FATAL_ERROR "Could not find VIENNACL")
    ENDIF (VIENNACL_FIND_REQUIRED)
ENDIF (VIENNACL_FOUND)

IF(VIENNACL_WITH_OPENCL)
    find_package(OpenCL)
ENDIF(VIENNACL_WITH_OPENCL)

LIST( APPEND VIENNACL_INCLUDE_DIRS ${VIENNACL_INCLUDE_DIR} ${OPENCL_INCLUDE_DIRS} )

IF(VIENNACL_WITH_OPENCL)
    LIST( APPEND VIENNACL_LIBRARIES ${OPENCL_LIBRARIES} )
    LIST( REMOVE_DUPLICATES VIENNACL_LIBRARIES )
ENDIF(VIENNACL_WITH_OPENCL)

LIST( REMOVE_DUPLICATES VIENNACL_INCLUDE_DIRS )


SET( HAVE_VIENNACL TRUE )
message(STATUS "ViennaCL detected: " ${VIENNACL_INCLUDE_DIRS})

MARK_AS_ADVANCED(
        VIENNACL_INCLUDE_DIR
        VIENNACL_INCLUDE_DIRS
        VIENNACL_LIBRARIES
)