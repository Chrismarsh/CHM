#INCLUDE_DIRECTORIES("${PROJECT_SOURCE_DIR}/")
#LIST(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/tools/cmake/")

#clear this from when we setup meteoio as they reuse this variable
SET(plugins_sources  plugins/SnowpackIO.cc)
#for now, we always compile with AsciiIO and SmetIO
SET(plugins_sources ${plugins_sources} plugins/AsciiIO.cc)
SET(plugins_sources ${plugins_sources} plugins/SmetIO.cc)

IF(PLUGIN_IMISIO)
	SET(plugins_sources ${plugins_sources} plugins/ImisDBIO.cc)
	FIND_PACKAGE(OCCI REQUIRED)
	INCLUDE_DIRECTORIES(SYSTEM ${OCCI_INCLUDE_DIR})
	LINK_DIRECTORIES(${OCCI_LIBRARIES}) #this, so that if the users types in paths, it works
ENDIF(PLUGIN_IMISIO)

IF(PLUGIN_CAAMLIO)
	SET(plugins_sources ${plugins_sources} plugins/CaaMLIO.cc)
	FIND_PACKAGE(LibXml2 REQUIRED)
	INCLUDE_DIRECTORIES(SYSTEM ${LIBXML2_INCLUDE_DIR})
	LINK_DIRECTORIES(${LIBXML2_LIBRARIES})
ENDIF(PLUGIN_CAAMLIO)
