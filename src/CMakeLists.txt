
set(FILTER_SRCS
		filters/macdonald_undercatch.cpp

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
		${SNOBAL_SRC}

		)


set(CHM_SRCS
		main.cpp
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
		interpolation/interpolation.cpp)

set (HEADER_FILES
        mesh
		ui
        modules
        modules/interp_met
		modules/snobal
        libmaw
        interpolation
        timeseries
        utility
		filters
        ${CMAKE_CURRENT_SOURCE_DIR}
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
        ${METEOIO_INCLUDE_DIR}
		${CURSES_INCLUDE_DIR}

    )

#set(CURSES_LIBRARIES "${CURSES_LIBRARIES} /usr/lib64/libpanel.so")
set( LIB_FILES
	    ${Boost_LIBRARIES} #old var
        ${TBB_LIBRARIES}
        ${GSL_LIBRARIES}
	    ${CMAKE_THREAD_LIBS_INIT}
        ${ARMADILLO_LIBRARIES} 
        ${SUNDIALS_LIBRARIES}
        ${CGAL_3RD_PARTY_LIBRARIES}
        ${VTK_LIBRARIES}
        ${METEOIO_LIBRARY}
        ${LAPACK_LIBRARIES}
		${CURSES_LIBRARIES}
)


#are we linking against matlab?
if(MATLAB)
    set(LIBMAW_SRCS
        libmaw/graphics.cpp
        libmaw/matlab_engine.cpp modules/snobal/snomacros.h)

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
	${CHM_SRCS}
	${FILTER_SRCS}
    ${MODULE_SRCS}
	${LIBMAW_SRCS}
 )

target_link_libraries(
	CHM
	${LIB_FILES}
)
