
#http://www.cmake.org/cmake/help/v2.8.10/cmake.html#variable:CMAKE_LANG_COMPILER_ID
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -xhost -axsse4.2 -opt-prefetch -std=c++11")

find_package(CGAL REQUIRED)
if(CGAL_FOUND)
    message(STATUS "Found CGAL")
endif()
message(STATUS "${CGAL_USE_FILE}")
include(${CGAL_USE_FILE}) #as per https://www.cgal.org/releases.html release 4.2



set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O3 -std=c++11 ") 


#C#####
#GAL requires strict rounding
if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -fp-model strict") #-frounding-math needed for cgal
endif()

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" OR "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -frounding-math")
endif()
########

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++ -Wno-c++11-extensions")
  set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -stdlib=libc++")
endif()



find_package(OpenMP)

if(OPENMP_FOUND)
    if(NOT "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang") #Clang doesn't have mainline support for openmp 
#        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

#for matlab to run you need csh installed
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



message(STATUS "CHM CXX_FLAGS: ${CMAKE_CXX_FLAGS}")

set(JSON_SPIRIT_SRCS
	json_spirit/json_spirit_reader.cpp 
	json_spirit/json_spirit_value.cpp 
	json_spirit/json_spirit_writer.cpp 
 )

set(MODULE_SRCS
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
        global.cpp
        station.cpp
        timer.cpp
        variable_map.cpp
        
#        mesh/mesh.cpp
#        mesh/bounding_rect.cpp
#        mesh/triangle.cpp
        mesh/triangulation.cpp

        interpolation/interp_2d.cpp
        interpolation/inv_dist.cpp
        interpolation/interp_t_air.cpp
        interpolation/interp_rh.cpp
        interpolation/spline.cpp

        timeseries/timestep.cpp
        timeseries/timeseries.cpp
        timeseries/daily.cpp

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
        ${CGAL_INCLUDE_DIRS}
${CGAL_3RD_PARTY_INCLUDE_DIRS}
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
