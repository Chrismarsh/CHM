
#http://www.cmake.org/cmake/help/v2.8.10/cmake.html#variable:CMAKE_LANG_COMPILER_ID
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -xhost -axsse4.2 -opt-prefetch -std=c++11")


set(JSON_SPIRIT_SRCS
	json_spirit/json_spirit_reader.cpp 
	json_spirit/json_spirit_value.cpp 
	json_spirit/json_spirit_writer.cpp 
 )

set(MODULE_SRCS

        modules/interp_met/const_llra_ta.cpp
        modules/interp_met/Liston_monthly_llra_ta.cpp
        modules/interp_met/Liston_monthly_llra_rh.cpp
		modules/interp_met/Liston_wind.cpp

        modules/solar.cpp
        modules/Marsh_shading_iswr.cpp

        modules/Sicart_ilwr.cpp
        modules/PenmanMonteith_evaporation.cpp

        modules/precip.cpp

        modules/Walcek_atm_trans.cpp
)


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
	json_spirit
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
        ${EIGEN3_INCLUDE_DIR}
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
        ${EIGEN3_LIBRARIES}
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
	${JSON_SPIRIT_SRCS}
        ${LIBMAW_SRCS}
 )

target_link_libraries(
	CHM
	${LIB_FILES}
)
