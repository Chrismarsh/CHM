
#snowpack specific
SET(_versionString "2.8.0")
SET(_VERSION "-D_VERSION=\\\"${_versionString}\\\"")
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${_VERSION}")

INCLUDE(${CMAKE_CURRENT_SOURCE_DIR}/modules/snowpack/CMakeLists.txt)
SET(SNOWPACK_SRCS "")
FOREACH(f ${snowpacklib_sources})
#	LIST(APPEND SNOWPACK_SRCS "${CMAKE_CURRENT_SOURCE_DIR}/modules/snowpack/${f}")
LIST(APPEND SNOWPACK_SRCS "modules/snowpack/${f}")
ENDFOREACH(f)

set(FILTER_SRCS
		filters/macdonald_undercatch.cpp
		filters/scale_wind_speed.cpp
		filters/goodison_undercatch.cpp
		filters/debias_lw.cpp
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
		modules/interp_met/iswr_from_nwp.cpp
		modules/interp_met/Dodson_NSA_ta.cpp
		modules/interp_met/Thornton_var_p.cpp
		modules/interp_met/kunkel_rh.cpp
		modules/interp_met/rh_from_obs.cpp
		modules/interp_met/Longwave_from_obs.cpp
		modules/interp_met/Longwave_from_obs.cpp
		modules/interp_met/Dist_tlapse.cpp
		modules/interp_met/rh_no_lapse.cpp
		modules/interp_met/t_no_lapse.cpp
		modules/interp_met/p_no_lapse.cpp
		modules/interp_met/p_lapse.cpp
		modules/interp_met/lw_no_lapse.cpp
		modules/interp_met/uniform_wind.cpp
		modules/interp_met/MS_wind.cpp
		modules/interp_met/WindNinja.cpp
		modules/interp_met/Winstral_parameters.cpp
		modules/interp_met/t_monthly_lapse.cpp
		modules/interp_met/Cullen_monthly_llra_ta.cpp

		modules/point_mode.cpp

        modules/iswr.cpp
		modules/fetchr.cpp
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
		modules/scale_wind_vert.cpp
		modules/sub_grid.cpp
    	modules/snow_slide.cpp
		${SNOBAL_SRC}
		modules/threshold_p_phase.cpp
		${SNOWPACK_SRCS}
		modules/solar.cpp
		modules/fast_shadow.cpp
        modules/deform_mesh.cpp
		modules/crop_rotation.cpp
        modules/PBSM3D.cpp

		CACHE INTERNAL "" FORCE )


set(CHM_SRCS
		#main.cpp needs to be added below so we can re use CHM_SRCS in the gtest build
		core.cpp
		global.cpp
		station.cpp
		metdata.cpp

		physics/Atmosphere.cpp

		mesh/triangulation.cpp

		interpolation/inv_dist.cpp
		interpolation/TPSpline.cpp
		interpolation/nearest.cpp

		timeseries/timestep.cpp
		timeseries/timeseries.cpp
		timeseries/daily.cpp
		timeseries/netcdf.cpp
		timeseries/variablestorage.cpp

		utility/regex_tokenizer.cpp
		utility/timer.cpp
		utility/jsonstrip.cpp
		utility/readjson.cpp

		interpolation/interpolation.cpp
        math/coordinates.cpp

		CACHE INTERNAL "" FORCE)


set (HEADER_FILES
        mesh
        modules
		physics
        modules/interp_met
		modules/snobal
        libmaw
        interpolation
        timeseries
        utility
		filters
        ${CMAKE_CURRENT_SOURCE_DIR}
		${CMAKE_BINARY_DIR}/src # for clion generated files / out of source builds
		CACHE INTERNAL "" FORCE
		)


set(EXT_TARGETS
		netcdf-cxx::netcdf-cxx
		Boost::Boost
		tbb::tbb
		gsl::gsl
		Threads::Threads
		armadillo::armadillo
		vtk::vtk
#		Iconv::Iconv # This is for unicode  support in gdal and is not conssitently brought in when building gdal
		gdal::gdal #proj is brought it via this
		viennacl::viennacl
		cgal::cgal
		hdf5::hdf5
		sparsehash::sparsehash
		meteoio::meteoio

		${ARPACK_LIBRARY}
		eigen3::eigen3
		CACHE INTERNAL "" FORCE)

#handle the optionals
if(OpenMP_FOUND)
	list( APPEND EXT_TARGETS	OpenMP::OpenMP_CXX)
endif()

if(USE_OCL AND OPENCL_FOUND)
	list( APPEND EXT_TARGETS OpenCL::OpenCL)
endif()

if(gperftools_FOUND)
	list( APPEND EXT_TARGETS gperftools::gperftools)
endif()

if(MPI_FOUND AND USE_MPI)
	list( APPEND EXT_TARGETS MPI::MPI_CXX)
endif()

#are we linking against matlab?
if(MATLAB)
    set(LIBMAW_SRCS
        libmaw/graphics.cpp
        libmaw/matlab_engine.cpp
		CACHE INTERNAL "" FORCE )

    list( APPEND EXT_TARGETS
            ${MATLAB_ENG_LIBRARY}
            ${MATLAB_MX_LIBRARY})

    list( APPEND HEADER_FILES
        ${MATLAB_INCLUDE_DIR})

endif()

add_executable(
	CHM
	main.cpp # this needs to be here so we can reuse CHM_SRCS in the gtest build. but since it links it's own main, we cannot have this one.
	${CHM_SRCS}
	${FILTER_SRCS}
    ${MODULE_SRCS}
	${LIBMAW_SRCS}
	)

target_include_directories(CHM PRIVATE ${MPI_CXX_INCLUDE_PATH} ${HEADER_FILES} )
target_compile_options(CHM PRIVATE ${MPI_CXX_COMPILE_FLAGS})
target_compile_features(CHM PRIVATE cxx_std_14)
target_link_libraries(
		CHM
		${EXT_TARGETS}
)


if(STATIC_ANLAYSIS)
	include( ${CMAKE_SOURCE_DIR}/CMake/PVS-Studio.cmake)
	pvs_studio_add_target(TARGET analyze ALL
			OUTPUT FORMAT errorfile
			ANALYZE CHM
			LICENSE "/home/chris/Downloads/PVS-Studio.lic"
			LOG CHM.plog
			)
endif()

if (BUILD_TESTS)

	add_subdirectory (tests/googletest/googletest)
	include_directories(${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

	set(TEST_SRCS
			#   test_station.cpp
			tests/test_interpolation.cpp
			tests/test_timeseries.cpp
			tests/test_core.cpp
			tests/test_variablestorage.cpp
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

	target_include_directories(runUnitTests PRIVATE ${MPI_CXX_INCLUDE_PATH} ${HEADER_FILES})
    target_compile_options(runUnitTests PRIVATE ${MPI_CXX_COMPILE_FLAGS})

	target_link_libraries(
			runUnitTests
			${EXT_TARGETS}
			${GTEST_LINK}
			${MPI_CXX_LIBRARIES}
			${MPI_CXX_LINK_FLAGS}
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
