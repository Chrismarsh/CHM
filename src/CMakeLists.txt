
#http://www.cmake.org/cmake/help/v2.8.10/cmake.html#variable:CMAKE_LANG_COMPILER_ID
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -xhost -axsse4.2 -opt-prefetch -std=c++11")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -std=c++11") 


if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
  set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -stdlib=libc++")
endif()


#sym link /opt/matlab-7sp1 to the real matlab root will make this work properly
find_package(Matlab REQUIRED)

find_package(Armadillo 3.0 REQUIRED)

find_package(Boost
             1.54.0
             COMPONENTS
             system
             filesystem
             date_time
             log
             thread
             regex
             REQUIRED)

#need this for calling the logger from multiple modules
ADD_DEFINITIONS(-DBOOST_LOG_DYN_LINK)

find_package (Threads REQUIRED)

find_package(TBB REQUIRED)

find_package(GSL REQUIRED)


set(JSON_SPIRIT_SRCS
	json_spirit/json_spirit_reader.cpp 
	json_spirit/json_spirit_value.cpp 
	json_spirit/json_spirit_writer.cpp 
 )

set(MODULE_SRCS
        modules/TestModule.cpp
        modules/solar.cpp
)

set(LIBMAW_SRCS
    libmaw/graphics.cpp
    libmaw/matlab_engine.cpp
)
set(CHM_SRCS 
	main.cpp
	core.cpp
        module_factory.cpp
        
        station.cpp
        timer.cpp

        mesh/mesh.cpp
        mesh/bounding_rect.cpp
        mesh/triangle.cpp
        mesh/triangulation.cpp

        interpolation/interp_2d.cpp
        interpolation/inv_dist.cpp
        interpolation/interp_t_air.cpp
        interpolation/interp_rh.cpp

        timeseries/timestep.cpp
        timeseries/timeseries.cpp

        utility/crc_hash_compare.cpp
        utility/regex_tokenizer.cpp
)
set (HEADER_FILES
	json_spirit
        mesh
        modules
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
        ${MATLAB_INCLUDE_DIR}
    )
set( LIB_FILES
	${Boost_LIBRARIES}
        ${TBB_LIBRARIES}
        ${GSL_LIBRARIES}
	${CMAKE_THREAD_LIBS_INIT}
        ${ARMADILLO_LIBRARIES} 
	${MATLAB_ENG_LIBRARY} 
	${MATLAB_MX_LIBRARY}  
)

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
