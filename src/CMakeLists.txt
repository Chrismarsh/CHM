
#http://www.cmake.org/cmake/help/v2.8.10/cmake.html#variable:CMAKE_LANG_COMPILER_ID
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g -O0 -std=c++11")



if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
  set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -stdlib=libc++")
endif()

#find_package(CGAL REQUIRED)

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
)
 
set(CHM_SOURCES 
	main.cpp
	mesh.cpp
	core.cpp
        module_factory.cpp
        forcing_data.cpp
        regex_tokenizer.cpp
        station.cpp
)
set (HEADER_FILES
	json_spirit
        modules
	${BOOST_INCLUDE_PATH}
        ${TBB_INCLUDE_DIRS}
        ${GSL_INCLUDE_DIRS}
)
set( LIB_FILES
	${Boost_LIBRARIES}
        ${TBB_LIBRARIES}
        ${GSL_LIBRARIES}
	${CMAKE_THREAD_LIBS_INIT}
)

include_directories(${HEADER_FILES})

add_executable(
	CHM 
	${CHM_SOURCES} 
        ${MODULE_SRCS}
	${JSON_SPIRIT_SRCS} 
)

target_link_libraries(
	CHM
	${LIB_FILES}
)