
# GTest
add_subdirectory (gtest)
include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

set(TEST_SRCS
 #   test_station.cpp
#    test_timeseries.cpp
    test_core.cpp
#    test_mesh.cpp
    test_regexptokenizer.cpp
#    test_daily.cpp
)

set(MODULE_SRCS

		../src/modules/interp_met/const_llra_ta.cpp
		../src/modules/interp_met/Liston_monthly_llra_ta.cpp
		../src/modules/interp_met/Kunkel_monthlyTd_rh.cpp
		../src/modules/interp_met/Liston_wind.cpp

		../src/modules/slope_iswr.cpp
		../src/modules/Marsh_shading_iswr.cpp
		../src/modules/Sicart_ilwr.cpp
		../src/modules/PenmanMonteith_evaporation.cpp
		../src/modules/interp_met/Thornton_p.cpp
		../src/modules/Walcek_cloud.cpp
		../src/modules/Harder_precip_phase.cpp
		../src/modules/Burridge_iswr.cpp
		../src/modules/Iqbal_iswr.cpp
		../src/modules/interp_met/iswr_from_obs.cpp
		../src/modules/interp_met/Dodson_NSA_ta.cpp
		../src/modules/interp_met/Thornton_var_p.cpp
		../src/modules/Richard_albedo.cpp
		)

if(MATLAB)
set(LIBMAW_SRCS
    ../src/libmaw/graphics.cpp
    ../src/libmaw/matlab_engine.cpp
)
endif()

set(CHM_SRCS

		../src/core.cpp
		../src/module_factory.cpp
		../src/global.cpp
		../src/station.cpp
		../src/variable_map.cpp

		../src/mesh/triangulation.cpp

		../src/interpolation/inv_dist.cpp
		../src/interpolation/TPSpline.cpp

		../src/timeseries/timestep.cpp
		../src/timeseries/timeseries.cpp
		../src/timeseries/daily.cpp

		../src/utility/crc_hash_compare.cpp
		../src/utility/regex_tokenizer.cpp
		../src/utility/timer.cpp
		)
set (HEADER_FILES
    ../src/mesh
    ../src/modules
	../src/modules/interp_met
    ../src/libmaw
    ../src/interpolation
    ../src/timeseries
    ../src/utility
    ../src
    ${CMAKE_CURRENT_SOURCE_DIR}
    )

set (EXT_HEADER_FILES
	${BOOST_INCLUDE_PATH}
        ${TBB_INCLUDE_DIRS}
        ${GSL_INCLUDE_DIRS}
        ${AMRADILLO_INCLUDE_DIR}
        ${MATLAB_INCLUDE_DIR}
        ${CGAL_INCLUDE_DIRS}
        ${CGAL_3RD_PARTY_INCLUDE_DIRS}
		${METEOIO_INCLUDE_DIR}
    )
set( LIB_FILES
	    ${Boost_LIBRARIES}
        ${TBB_LIBRARIES}
        ${GSL_LIBRARIES}
	    ${CMAKE_THREAD_LIBS_INIT}
        ${ARMADILLO_LIBRARIES} 
	    ${MATLAB_ENG_LIBRARY}
	    ${MATLAB_MX_LIBRARY}
        ${CGAL_LIBRARIES}
        ${CGAL_3RD_PARTY_LIBRARIES}
		${VTK_LIBRARIES}
		${METEOIO_LIBRARY}
		${LAPACK_LIBRARIES}
)

set( GTEST_LINK
    gtest
    gtest_main
)

include_directories(${HEADER_FILES} ${EXT_HEADER_FILES})

add_executable(
    runUnitTests
    ${CHM_SRCS}
    ${TEST_SRCS} 
    ${MODULE_SRCS}
    ${LIBMAW_SRCS}
)

target_link_libraries(
    runUnitTests
    ${LIB_FILES}
    ${GTEST_LINK}
    ${VTK_LIBRARIES}
)

# You can also omit NAME and COMMAND. The second argument could be some other
  # test executable.
  add_test(runUnitTests runUnitTests)

#  add_custom_target(check COMMAND ${CMAKE_CTEST_COMMAND}
#                  DEPENDS ${EXECUTABLE_OUTPUT_PATH}/runUnitTests)
 add_custom_target(check COMMAND ${EXECUTABLE_OUTPUT_PATH}/runUnitTests DEPENDS runUnitTests WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})
 message(STATUS "Copying files to test dir ${EXECUTABLE_OUTPUT_PATH} from ${CMAKE_SOURCE_DIR}/tests/files/")

add_custom_command(TARGET runUnitTests PRE_BUILD
                   COMMAND ${CMAKE_COMMAND} -E copy_directory
                       ${CMAKE_SOURCE_DIR}/tests/files $<TARGET_FILE_DIR:runUnitTests>)

