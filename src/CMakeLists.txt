#meteoio specific
SET(_versionString "2.5.0")
SET(_VERSION "-D_VERSION=\\\"${_versionString}\\\"")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${_VERSION}")
INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/meteoio/CMakeLists.txt)

SET(METEOIO_SRCS "")
FOREACH(f ${meteoio_sources})
#	LIST(APPEND METEOIO_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/meteoio/${f}")
	LIST(APPEND METEOIO_SRCS "meteoio/${f}")
ENDFOREACH(f)


#snowpack specific
INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/modules/snowpack/CMakeLists.txt)
SET(SNOWPACK_SRCS "")
FOREACH(f ${snowpacklib_sources})
#	LIST(APPEND SNOWPACK_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/modules/snowpack/${f}")
LIST(APPEND SNOWPACK_SRCS "modules/snowpack/${f}")
ENDFOREACH(f)

set(FILTER_SRCS
		filters/macdonald_undercatch.cpp

		CACHE INTERNAL "" FORCE
		)

set(SNOBAL_SRC
		modules/snobal/snobal.cpp #<-- main module files
		modules/snobal/sno.cpp
	)

set(MODULE_SRCS

		modules/interp_met/const_llra_ta.cpp
		modules/interp_met/Liston_monthly_llra_ta.cpp
		modules/interp_met/Kunkel_monthlyTd_rh.cpp
		modules/interp_met/Liston_wind.cpp
		modules/interp_met/iswr_from_obs.cpp
		modules/interp_met/Dodson_NSA_ta.cpp
		modules/interp_met/Thornton_var_p.cpp
		modules/interp_met/kunkel_rh.cpp
		modules/interp_met/rh_from_obs.cpp
        modules/interp_met/Longwave_from_obs.cpp

		modules/point_mode.cpp

		modules/slope_iswr.cpp
		modules/Marsh_shading_iswr.cpp
		modules/Sicart_ilwr.cpp
		modules/PenmanMonteith_evaporation.cpp
		modules/interp_met/Thornton_p.cpp
		modules/Walcek_cloud.cpp
		modules/Harder_precip_phase.cpp
		modules/Burridge_iswr.cpp
		modules/Iqbal_iswr.cpp
		modules/Richard_albedo.cpp
		modules/snowpack.cpp
		modules/Gray_inf.cpp
		modules/Simple_Canopy.cpp
		modules/PhysConsts.h
		${SNOBAL_SRC}

		modules/threshold_p_phase.cpp
		${SNOWPACK_SRCS}
		CACHE INTERNAL "" FORCE
		)


set(CHM_SRCS
		#main.cpp needs to be added below so we can re use CHM_SRCS in the gtest build
		core.cpp
		module_factory.cpp
		filter_factory.cpp
		global.cpp
		station.cpp
		variable_map.cpp

		ui/ui.cpp
		ui/ncstream.cpp

		mesh/triangulation.cpp

		interpolation/inv_dist.cpp
		interpolation/TPSpline.cpp

		timeseries/timestep.cpp
		timeseries/timeseries.cpp
		timeseries/daily.cpp

		utility/crc_hash_compare.cpp
		utility/regex_tokenizer.cpp
		utility/timer.cpp
		utility/jsonstrip.cpp
		utility/readjson.cpp

		interpolation/interpolation.cpp
		${METEOIO_SRCS}
		CACHE INTERNAL "" FORCE
		)


set (HEADER_FILES
        mesh
		ui
        modules
        modules/interp_met
		modules/snobal
        libmaw
		meteoio
        interpolation
        timeseries
        utility
		filters
        ${CMAKE_CURRENT_SOURCE_DIR}
		CACHE INTERNAL "" FORCE
		)


set (EXT_HEADER_FILES
		${BOOST_INCLUDE_PATH} #old version of boost???
		${Boost_INCLUDE_DIR}  #watch caps!!
        ${TBB_INCLUDE_DIRS}
        ${GSL_INCLUDE_DIRS}
        ${AMRADILLO_INCLUDE_DIR}
        ${SUNDIALS_INCLUDES}
        ${CGAL_INCLUDE_DIRS}
        ${CGAL_3RD_PARTY_INCLUDE_DIRS}
		${CURSES_INCLUDE_DIR}

		CACHE INTERNAL "" FORCE
    )

set( LIB_FILES
	    ${Boost_LIBRARIES}
        ${TBB_LIBRARIES}
        ${GSL_LIBRARIES}
	    ${CMAKE_THREAD_LIBS_INIT}
        ${ARMADILLO_LIBRARIES} 
        ${SUNDIALS_LIBRARIES}
        ${CGAL_3RD_PARTY_LIBRARIES}
        ${VTK_LIBRARIES}
        ${LAPACK_LIBRARIES}
		${CURSES_LIBRARIES}

		CACHE INTERNAL "" FORCE
)


#are we linking against matlab?
if(MATLAB)
    set(LIBMAW_SRCS
        libmaw/graphics.cpp
        libmaw/matlab_engine.cpp )

    set( LIB_FILES
            ${LIB_FILES}
            ${MATLAB_ENG_LIBRARY}
            ${MATLAB_MX_LIBRARY})

    set(EXT_HEADER_FILES
        ${EXT_HEADER_FILES}
        ${MATLAB_INCLUDE_DIR})

endif()



include_directories(${HEADER_FILES} ${EXT_HEADER_FILES})

link_directories(${TBB_LIBRARY_DIRS}) # as per https://code.google.com/p/findtbb/wiki/HowTo
add_executable(
	CHM
	main.cpp # this needs to be here so we can reuse CHM_SRCS in the gtest build. but since it links it's own main, we cannot have this one.
	${CHM_SRCS}
	${FILTER_SRCS}
    ${MODULE_SRCS}
	${LIBMAW_SRCS}
 )

target_link_libraries(
	CHM
	${LIB_FILES}
)

if (test)

	add_subdirectory (tests/googletest/googletest)
	include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

	set(TEST_SRCS
			#   test_station.cpp
			tests/test_interpolation.cpp
			tests/test_timeseries.cpp
			tests/test_core.cpp
			#    test_mesh.cpp
			tests/test_regexptokenizer.cpp
			#    test_daily.cpp
			)


	set( GTEST_LINK
			gtest
			gtest_main

			)

	add_executable(
			runUnitTests
			${CHM_SRCS}
			${FILTER_SRCS}
			${MODULE_SRCS}
			${LIBMAW_SRCS}
			${TEST_SRCS}
	)

	target_link_libraries(
			runUnitTests
			${LIB_FILES}
			${GTEST_LINK}
	)

	# You can also omit NAME and COMMAND. The second argument could be some other
	# test executable.
	add_test(runUnitTests runUnitTests)

	add_custom_target(check COMMAND ${EXECUTABLE_OUTPUT_PATH}/runUnitTests DEPENDS runUnitTests WORKING_DIRECTORY ${EXECUTABLE_OUTPUT_PATH})


	set(COPY_TO_PATH ${EXECUTABLE_OUTPUT_PATH})

	add_custom_command(TARGET runUnitTests PRE_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy_directory
			${CMAKE_SOURCE_DIR}/src/tests/files
			$<TARGET_FILE_DIR:runUnitTests>
			COMMENT "Copying files to $<TARGET_FILE_DIR:runUnitTests> from ${CMAKE_SOURCE_DIR}/src/tests/files/")


endif()
