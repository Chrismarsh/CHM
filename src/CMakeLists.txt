

set(FILTER_SRCS
		filters/macdonald_undercatch.cpp
		filters/scale_wind_speed.cpp
		filters/goodison_undercatch.cpp
		filters/debias_lw.cpp
		CACHE INTERNAL "" FORCE
		)



set(MODULE_SRCS

		modules/interp_met/const_llra_ta.cpp
		modules/interp_met/Liston_monthly_llra_ta.cpp
		modules/interp_met/Kunkel_monthlyTd_rh.cpp
		modules/interp_met/Liston_wind.cpp
		modules/interp_met/iswr_from_obs.cpp
		modules/interp_met/iswr_from_nwp.cpp
		modules/interp_met/Dodson_NSA_ta.cpp
		modules/interp_met/p_from_obs.cpp
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
		modules/threshold_p_phase.cpp
		modules/solar.cpp
		modules/fast_shadow.cpp
        modules/deform_mesh.cpp
		modules/crop_rotation.cpp
        modules/PBSM3D.cpp
		modules/snobal.cpp
		modules/fsm.cpp



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
		math
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
		${ARPACK_LIBRARY}
		Threads::Threads
		${GDAL}::${GDAL}
		${GSL}::${GSL}
		${PROJ}::${PROJ}
		${CGAL}::${CGAL}
		${HDF5}::${HDF5}

		${METEOIO}::${METEOIO}
		${TRILINOS}::${TRILINOS}
		${SPARSEHASH}::${SPARSEHASH}

		${FUNC}::${FUNC}
		CACHE INTERNAL "" FORCE)


if(BUILD_WITH_CONAN)
	list( APPEND EXT_TARGETS
			netcdf-cxx::netcdf-cxx
			TBB::tbb
			Boost::Boost
			Eigen3::Eigen3
			armadillo::armadillo
			vtk::vtk

			)
else()
	list( APPEND EXT_TARGETS
			NetCDF::NetCDF
			TBB::tbb
			ARMADILLO::ARMADILLO
			Eigen3::Eigen

			VTK::CommonCore
			VTK::CommonDataModel
			VTK::CommonTransforms
			VTK::FiltersCore
			VTK::FiltersGeneral
			VTK::FiltersGeometry
			VTK::FiltersSources
			VTK::IOXML

			Boost::headers
			Boost::system
			Boost::filesystem
			Boost::date_time
			Boost::log
			Boost::thread
			Boost::regex
			Boost::iostreams
			Boost::program_options
			)
endif()

#handle the optionals
if(OpenMP_FOUND)
	list( APPEND EXT_TARGETS	OpenMP::OpenMP_CXX)
endif()

if(USE_OCL AND OPENCL_FOUND)
	list( APPEND EXT_TARGETS OpenCL::OpenCL)
endif()

if(${GPERFTOOLS}_FOUND)
	list( APPEND EXT_TARGETS ${GPERFTOOLS}::${GPERFTOOLS})
endif()

if(USE_JEMALLOC)
	list( APPEND EXT_TARGETS Jemalloc::Jemalloc)
endif()

if(MPI_FOUND AND USE_MPI)
	list( APPEND EXT_TARGETS MPI::MPI_CXX)
	if(NOT BUILD_WITH_CONAN)
		list( APPEND EXT_TARGETS Boost::mpi)
	endif()
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

# when we run make install all the RPATHs will be stripped out and set to these
# n.b. $ORIGIN should have its dollar sign escaped with a backslash to have it end up with proper
# syntax in the final executable.
# https://gitlab.kitware.com/cmake/community/-/wikis/doc/cmake/RPATH-handling
# !!!  these must be done before targets are defined in add_executable !!!!
if (APPLE)
	set(CMAKE_INSTALL_RPATH "@executable_path/../lib")
else()
	set(CMAKE_INSTALL_RPATH "\$ORIGIN/../lib")
endif()

add_library(CHMmath OBJECT
  math/LinearAlgebra.cpp
  )

target_link_libraries(CHMmath
  ${EXT_TARGETS}
  )

target_include_directories(CHMmath PUBLIC
  ${HEADER_FILES}
)

set_target_properties(CHMmath
		PROPERTIES
		COMPILE_FLAGS ${CHM_BUILD_FLAGS})


add_executable(
		CHM
		main.cpp # this needs to be here so we can reuse CHM_SRCS in the gtest build. but since it links it's own main, we cannot have this one.
		${CHM_SRCS}
		${FILTER_SRCS}
		${MODULE_SRCS}
		${LIBMAW_SRCS}
)


# target_include_directories(CHM PRIVATE ${MPI_CXX_INCLUDE_PATH} ${HEADER_FILES} ${THIRD_PARTY_INCL_PATHS} )
target_include_directories(CHM PRIVATE ${HEADER_FILES} )
if(MPI_FOUND AND USE_MPI)
  target_include_directories(CHM PRIVATE ${MPI_CXX_INCLUDE_PATH} )
  target_compile_options(CHM PRIVATE ${MPI_CXX_COMPILE_FLAGS})
endif()


target_compile_features(CHM PRIVATE cxx_std_14)
target_link_libraries(
		CHM
		CHMmath
		${EXT_TARGETS}
		${THIRD_PARTY_TARGETS}
)
set_target_properties(CHM
		PROPERTIES
		RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
		COMPILE_FLAGS ${CHM_BUILD_FLAGS}
		UNITY_BUILD OFF
		UNITY_BUILD_BATCH_SIZE 2
		)


# the patches suggested by conan https://docs.conan.io/en/latest/howtos/manage_shared_libraries/rpaths.html
# don't work for the gdal target (despite working elsewhere) for some reason. So, the strategy here is:
# 1) The gdal libraries are patched at the conan package stage to fix the rpath issue
# 2) Patch the CHM binary to ensure that we are absolutely looking for @rpath/dylibs
if(APPLE)
	add_custom_command(TARGET CHM POST_BUILD
			#patch CHM to look for @rpath/ for all dylibs that aren't looking for @rpath
			COMMAND bash -c "otool -l ${CMAKE_BINARY_DIR}/bin/CHM | grep name | grep -v segname |  grep -v sectname | grep -v @rpath  | awk '{print $2}' | grep -v '^/' | while read x; do install_name_tool -change $x @rpath/`echo $x | grep -Eo '[a-zA-Z0-9_\.-]+\.dylib'` ${CMAKE_BINARY_DIR}/bin/CHM; done"

			#patch all our lib to look in @rpath
			COMMAND bash -c "ls ${CMAKE_BINARY_DIR}/lib | while read f; do otool -l ${CMAKE_BINARY_DIR}/lib/$f | grep name | grep -v segname |  grep -v sectname | grep -v @rpath  | awk '{print $2}' | grep -v '^/' | while read x; do install_name_tool -change $x @rpath/`echo $x | grep -Eo '[a-zA-Z0-9_\.-]+\.dylib'` ${CMAKE_BINARY_DIR}/lib/$f; done; done"

#			COMMAND bash -c "ls ${CMAKE_BINARY_DIR}/lib | while read f; do otool -l ${CMAKE_BINARY_DIR}/lib/$f | grep name | grep -v segname |  grep -v sectname | grep -v @rpath  | awk '{print $2}' | grep '/Users' | while read x; do install_name_tool -change $x `echo $x | grep -Eo '[a-zA-Z0-9_\.-]+\.dylib'` ${CMAKE_BINARY_DIR}/lib/$f; done; done"

			# These are likely no longer needed
				#patch netcdf-cxx to look for an @rpath
				#COMMAND bash -c "nccxx=`otool -l ${CMAKE_BINARY_DIR}/bin/CHM |grep libnetcdf-cxx | awk '{print $2}' | grep -Eo '[a-zA-Z0-9_+\.-]+\.dylib'`; nc=`otool -l ${CMAKE_BINARY_DIR}/bin/CHM |grep libnetcdf | grep -v cxx | awk '{print $2}' | grep -Eo '[a-zA-Z0-9_\.-]+\.dylib'`; install_name_tool -change $nc @rpath/$nc ${CMAKE_BINARY_DIR}/lib/$nccxx"

			#strip all LC_RPATH from our non-install binary
			COMMAND bash -c "otool -l ${CMAKE_BINARY_DIR}/bin/CHM | grep LC_RPATH -A2 | grep path | awk '{print $2}' | while read x; do install_name_tool -delete_rpath $x ${CMAKE_BINARY_DIR}/bin/CHM; done"

			#strip all LC_RPATH from our libs
			COMMAND bash -c "ls ${CMAKE_BINARY_DIR}/lib | while read f; do otool -l ${CMAKE_BINARY_DIR}/lib/$f | grep LC_RPATH -A2 | grep path | awk '{print $2}' | while read x; do install_name_tool -delete_rpath $x ${CMAKE_BINARY_DIR}/lib/$f; done; done"

			#update LC_RPATH for our non-install binary
			COMMAND bash -c "install_name_tool -add_rpath @executable_path/../lib ${CMAKE_BINARY_DIR}/bin/CHM"

#			COMMAND bash -c "ls ${CMAKE_BINARY_DIR}/lib | while read f; do install_name_tool -add_rpath @executable_path ${CMAKE_BINARY_DIR}/lib/$f; done"

			VERBATIM
			)
else()
	# General note for this whole section: it's important to note that RUNPATH has become the default for ld, superceeding RPATH
	# The search paths are RPATH -> LD_LIBRARY_PATH -> RUNPATH
	# However, RPATH resolves transient libraries (eg., netcdf-cxx requires netcdf-c) whereas RUNPATH does NOT
	# https://stackoverflow.com/a/55283740/410074
	# https://bugs.debian.org/cgi-bin/bugreport.cgi?bug=847298

	# So the approach here will be to use the old style dtags -Wl,--disable-new-dtags
	# ie., use RPATH instead of RUNPATH
	# and update the RPATH to be like our install path. Note RPATH is not used if RUNPATH is present.
	target_link_options(CHM
			PUBLIC "LINKER:--disable-new-dtags" )

	add_custom_command(TARGET CHM POST_BUILD
#			# Some .so will have runpaths, some do not. Just remove them all. The RUNPATHS override the binarie's rpath (I don't know why)
			COMMAND bash -c "cd ${CMAKE_BINARY_DIR}/lib ; ls | while read x; do patchelf --remove-rpath $x; done"
#
#			# RPATH correctly works, so unlikely to be needed but kept just in case.
#			# This removes all runpaths , and then here we will patch each .so to have the correct RUNPATH
			COMMAND bash -c "cd ${CMAKE_BINARY_DIR}/lib ; ls | while read x; do patchelf --set-rpath '\$ORIGIN/../lib' $x; done"
#
#			#see if there are any hard coded paths in your binary and if so replace with just the so name
#			COMMAND  bash -c "readelf ${CMAKE_BINARY_DIR}/bin/CHM -d | grep NEEDED |awk '{print $5}' | grep '\\[/'  | awk -F'[\\[\\]]' '{print $2}' | while read x; do patchelf --replace-needed $x `echo $x | grep -Eo '[a-zA-Z0-9_+\\.-]+so.*?'` ${CMAKE_BINARY_DIR}/bin/CHM; done"
#
#			VERBATIM
#			)

endif()

### Targets for the preprocessing tools
add_executable(
		partition
		preprocessing/partition/main.cpp # this needs to be here so we can reuse CHM_SRCS in the gtest build. but since it links it's own main, we cannot have this one.
		${CHM_SRCS}
)
target_compile_features(partition PRIVATE cxx_std_14)


target_link_libraries(
		partition
		CHMmath
		${EXT_TARGETS}
		${THIRD_PARTY_TARGETS}
)
set_target_properties(partition
		PROPERTIES
		RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin"
		COMPILE_FLAGS ${CHM_BUILD_FLAGS}
		)

#make install will correctly set the rpath for us to find the lib/ dir with the so/dylibs we need
install(TARGETS CHM RUNTIME)
install(TARGETS partition RUNTIME)

#if(BUILD_WITH_CONAN)
#	install(DIRECTORY ${CMAKE_BINARY_DIR}/lib/
#			DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
#endif()

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
	message(STATUS "Tests enabled. Run with make check")
	add_subdirectory (tests/googletest/googletest)

	set(TEST_SRCS
			tests/test_station.cpp
			tests/test_interpolation.cpp
			tests/test_timeseries.cpp
			tests/test_core.cpp
			tests/test_variablestorage.cpp
			tests/test_metdata.cpp
			tests/test_netcdf.cpp
			#    test_mesh.cpp
			tests/test_regexptokenizer.cpp
			#    test_daily.cpp
            tests/test_triangulation.cpp)


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

	target_include_directories(runUnitTests PRIVATE ${HEADER_FILES} )

	if(MPI_FOUND AND USE_MPI)
	  target_include_directories(runUnitTests PRIVATE ${MPI_CXX_INCLUDE_PATH} )
	  target_compile_options(runUnitTests PRIVATE ${MPI_CXX_COMPILE_FLAGS})
	endif()

	target_link_libraries(
			runUnitTests
			CHMmath
			${EXT_TARGETS}
			${GTEST_LINK}
			${THIRD_PARTY_TARGETS}
	)

	set_target_properties(runUnitTests
			PROPERTIES
			RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests"
			)

	#add the `make check` target that runs the unit tests
	set(TEST_DIR ${CMAKE_BINARY_DIR}/tests)
	add_custom_target(check COMMAND ${TEST_DIR}/runUnitTests
						DEPENDS runUnitTests
						WORKING_DIRECTORY ${TEST_DIR})

	if (NOT APPLE)
		target_link_options(runUnitTests
				PUBLIC "LINKER:--disable-new-dtags" )
	endif()
	#copy the data files from test_data to the dir that holds the test binary
	add_custom_command(TARGET runUnitTests PRE_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy_directory
			${CMAKE_SOURCE_DIR}/test_data
			$<TARGET_FILE_DIR:runUnitTests>
			COMMENT "Copying files to $<TARGET_FILE_DIR:runUnitTests> from ${CMAKE_SOURCE_DIR}/test_data/")



endif()
