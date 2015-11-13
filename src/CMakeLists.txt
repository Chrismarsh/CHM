

set(MODULE_SRCS

		modules/interp_met/const_llra_ta.cpp
		modules/interp_met/Liston_monthly_llra_ta.cpp
		modules/interp_met/Liston_monthly_llra_rh.cpp
		modules/interp_met/Liston_wind.cpp

		modules/slope_iswr.cpp
		modules/Marsh_shading_iswr.cpp
		modules/Sicart_ilwr.cpp
		modules/PenmanMonteith_evaporation.cpp
		modules/interp_met/precip.cpp
		modules/Walcek_cloud.cpp
		modules/Harder_precip_phase.cpp
		modules/Burridge_iswr.cpp modules/Iqbal_iswr.cpp modules/Iqbal_iswr.h modules/interp_met/iswr_from_obs.cpp modules/interp_met/iswr_from_obs.h modules/interp_met/Dodson_NSA_ta.cpp modules/interp_met/Dodson_NSA_ta.h)


set(CHM_SRCS
		main.cpp
		core.cpp
		module_factory.cpp
		global.cpp
		station.cpp
		variable_map.cpp

		mesh/triangulation.cpp

		interpolation/inv_dist.cpp
		interpolation/TPSpline.cpp

		timeseries/timestep.cpp
		timeseries/timeseries.cpp
		timeseries/daily.cpp

		utility/crc_hash_compare.cpp
		utility/regex_tokenizer.cpp
		utility/timer.cpp
		)

set (HEADER_FILES
        mesh
        modules
        modules/interp_met
        libmaw
        interpolation
        timeseries
        utility
        ${CMAKE_CURRENT_SOURCE_DIR}
    )

set (EXT_HEADER_FILES
	${BOOST_INCLUDE_PATH}
        ${TBB_INCLUDE_DIRS}
        ${GSL_INCLUDE_DIRS}
        ${AMRADILLO_INCLUDE_DIR}
        ${SUNDIALS_INCLUDES}
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
        ${SUNDIALS_LIBRARIES}
        ${CGAL_LIBRARIES}
        ${CGAL_3RD_PARTY_LIBRARIES}
        ${VTK_LIBRARIES}
        ${METEOIO_LIBRARY}
        ${LAPACK_LIBRARIES}

)


#are we linking against matlab?
if(MATLAB)
    set(LIBMAW_SRCS
        libmaw/graphics.cpp
        libmaw/matlab_engine.cpp)

    set( LIB_FILES
            ${LIB_FILES}
            ${MATLAB_ENG_LIBRARY}
            ${MATLAB_MX_LIBRARY})

    set(EXT_HEADER_FILES
        ${EXT_HEADER_FILES}
        ${MATLAB_INCLUDE_DIR})

endif()

include_directories(${HEADER_FILES} ${EXT_HEADER_FILES})

add_executable(
	CHM 
	${CHM_SRCS} 
        ${MODULE_SRCS}
        ${LIBMAW_SRCS}
 )

target_link_libraries(
	CHM
	${LIB_FILES}
)
