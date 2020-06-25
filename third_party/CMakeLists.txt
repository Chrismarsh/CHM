
#snowpack specific
SET(_versionString "2.8.0")
SET(_VERSION "-D_VERSION=\\\"${_versionString}\\\"")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${_VERSION}")

INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/third_party/snowpack/CMakeLists.txt)
SET(SNOWPACK_SRCS "")
FOREACH(f ${snowpacklib_sources})
    #	LIST(APPEND SNOWPACK_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/modules/snowpack/${f}")
    LIST(APPEND SNOWPACK_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/third_party/snowpack/${f}")
ENDFOREACH(f)


###################################
set(SNOBAL_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/snobal/sno.cpp
        )
###################################
set(FSM_SRC
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/FSM/FSM2_CHM.f90
        )

###################################
set(CHM_THIRD_PARTY
        ${SNOWPACK_SRCS}
        ${SNOBAL_SRC}
        ${FSM_SRC}
        )

set(THIRD_PARTY_INCL_PATHS
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/
        )