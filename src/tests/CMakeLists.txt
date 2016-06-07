#
## GTest
#add_subdirectory (gtest)
#include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
#
#set(TEST_SRCS
# #   test_station.cpp
#        test_interpolation.cpp
#        test_timeseries.cpp
#        test_core.cpp
##    test_mesh.cpp
#        test_regexptokenizer.cpp
##    test_daily.cpp
#)
#
#
#set( GTEST_LINK
#    gtest_main
#)
#
##meteoio specific
#SET(_versionString "2.5.0")
#SET(_VERSION "-D_VERSION=\\\"${_versionString}\\\"")
#SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${_VERSION}")
#INCLUDE(meteoio/CMakeLists.txt)
#
#set(SRCS
#        ${CHM_SRCS}
#        ${FILTER_SRCS}
#        ${MODULE_SRCS}
#        ${LIBMAW_SRCS})
#
##SET(ALL_SRC "")
##FOREACH(f ${SRCS})
##	LIST(APPEND ALL_SRC "${CMAKE_SOURCE_DIR}/src/${f}") #${CMAKE_SOURCE_DIR}
##ENDFOREACH(f)
#
##SET(ALL_HEADER "")
##FOREACH(f ${HEADER_FILES})
##    LIST(APPEND ALL_HEADER "${CMAKE_SOURCE_DIR}/src/${f}")
##ENDFOREACH(f)
##
##set (HEADER_FILES
##        ../src/mesh
##        ../src/ui
##        ../src/modules
##        ../src/modules/interp_met
##        ../src/modules/snobal
##        ../src/libmaw
##        ../src/interpolation
##        ../src/timeseries
##        ../src/utility
##        ../src
##        ../src/filters
##        .
##        )
#include_directories(${HEADER_FILES} ${EXT_HEADER_FILES})
#
#message(STATUS  ${HEADER_FILES})
#
#add_executable(
#    runUnitTests
#	${SRCS}
#    ${TEST_SRCS}
#)
#
#target_link_libraries(
#    runUnitTests
#    ${LIB_FILES}
#    ${GTEST_LINK}
#)
#
## You can also omit NAME and COMMAND. The second argument could be some other
#  # test executable.
#  add_test(runUnitTests runUnitTests)
#
##  add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
##                  DEPENDS ${EXECUTABLE_OUTPUT_PATH}/runUnitTests)
# add_custom_target(check COMMAND ${EXECUTABLE_OUTPUT_PATH}/runUnitTests DEPENDS runUnitTests WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
# message(STATUS "Copying files to test dir ${EXECUTABLE_OUTPUT_PATH} from ${CMAKE_SOURCE_DIR}/tests/files/")
#
#add_custom_command(TARGET runUnitTests PRE_BUILD
#                   COMMAND ${CMAKE_COMMAND} -E copy_directory
#                       ${CMAKE_SOURCE_DIR}/tests/files $<TARGET_FILE_DIR:runUnitTests>)
#
