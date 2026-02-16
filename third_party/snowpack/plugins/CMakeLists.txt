INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/")
LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/tools/cmake/")

CONFIGURE_FILE( ${PROJECT_SOURCE_DIR}/snowpack/plugins/SnowpackIO.cmake.cc ${PROJECT_SOURCE_DIR}/snowpack/plugins/SnowpackIO.cc )

SET(plugins_sources ${plugins_sources} plugins/SnowpackIO.cc)

#for now, we always compile with AsciiIO and SmetIO
SET(plugins_sources ${plugins_sources} plugins/AsciiIO.cc)
SET(plugins_sources ${plugins_sources} plugins/SmetIO.cc)

IF(PLUGIN_IMISIO)
	SET(plugins_sources ${plugins_sources} plugins/ImisDBIO.cc)
	FIND_PACKAGE(OCCI REQUIRED)
	INCLUDE_DIRECTORIES(SYSTEM ${OCCI_INCLUDE_DIR})
	LINK_DIRECTORIES(${OCCI_LIBRARIES})
	SET(extra_libs ${extra_libs} ${OCCI_LIBRARIES})
ENDIF(PLUGIN_IMISIO)

IF(PLUGIN_CAAMLIO)
	SET(plugins_sources ${plugins_sources} plugins/CaaMLIO.cc)
	SET(plugins_sources ${plugins_sources} plugins/pugixml/pugixml.cpp)
	set_source_files_properties(plugins/pugixml/pugixml.cpp PROPERTIES COMPILE_FLAGS "-w")
ENDIF(PLUGIN_CAAMLIO)

