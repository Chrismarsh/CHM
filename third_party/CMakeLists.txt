
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

add_library(SNOWPACK OBJECT
  ${SNOWPACK_SRCS}
  )

target_link_libraries(SNOWPACK PUBLIC
    MeteoIO::MeteoIO
    Boost::headers
    Boost::system
    Boost::filesystem
    Boost::date_time
    Boost::thread
    Boost::regex
    Boost::iostreams
    Boost::program_options
)


target_include_directories(SNOWPACK PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/src/    # needed for logger.hpp... that's a funny dependence, what is meant by "third party"?
)

###################################
add_library(SNOBAL OBJECT
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/snobal/sno.cpp
)
target_link_libraries(SNOBAL PUBLIC
    Func::Func
    Boost::headers
)
if(OpenMP_FOUND)
    target_link_libraries(SNOBAL PUBLIC OpenMP::OpenMP_CXX)
endif()
target_include_directories(SNOBAL PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/src/    # needed for logger.hpp... that's a funny dependence
)


set_target_properties( SNOBAL
        PROPERTIES
        COMPILE_FLAGS ${CHM_BUILD_FLAGS})

###################################
add_library(FSM OBJECT
        ${CMAKE_CURRENT_SOURCE_DIR}/third_party/FSM/FSM2_CHM.f90
	)
target_compile_options(FSM PRIVATE -std=f2003 -fno-underscoring -O3) # may need to be PUBLIC to work properly



###################################

set(THIRD_PARTY_TARGETS
  SNOWPACK
  SNOBAL
  FSM
  )

set_target_properties( SNOWPACK
        PROPERTIES
        COMPILE_FLAGS ${CHM_BUILD_FLAGS})


# third party includes in CHM src files to be specified relative to this directory
FOREACH(target ${THIRD_PARTY_TARGETS})
  target_include_directories(${target} INTERFACE
    ${CMAKE_CURRENT_SOURCE_DIR}/third_party/
    )

ENDFOREACH()
