include("LibFindMacros")

# Finally the library itself
GET_FILENAME_COMPONENT(SRC_DIR ${CMAKE_SOURCE_DIR} PATH) #ie goes up one level
STRING(REPLACE " " "\\ " SRC_DIR ${SRC_DIR})

IF(WIN32)
	GET_FILENAME_COMPONENT(LIBSNOWPACK_ROOT1 "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Microsoft\\Windows\\CurrentVersion\\Uninstall\\Snowpack;UninstallString]" PATH CACHE INTERNAL)
	GET_FILENAME_COMPONENT(LIBSNOWPACK_ROOT2 "[HKEY_LOCAL_MACHINE\\SOFTWARE\\Wow6432Node\\Microsoft\\Windows\\CurrentVersion\\Uninstall\\Snowpack;UninstallString]" PATH CACHE INTERNAL)
	GET_FILENAME_COMPONENT(LIBSNOWPACK_ROOT3 "[HKEY_LOCAL_MACHINE\\SOFTWARE\\WSL Institute for Snow and Avalanche Research\\Snowpack]" ABSOLUTE CACHE INTERNAL)
	SET(SEARCH_PATH
		ENV LIB
		${LIBSNOWPACK_ROOT1}/bin
		${LIBSNOWPACK_ROOT2}/bin
		${LIBSNOWPACK_ROOT3}/bin
		${SRC_DIR}/snowpack/bin
		"C:/Program Files/Snowpack/bin" )

	IF(MSVC)
		FIND_LIBRARY(LIBSNOWPACK_LIBRARY
			NAMES libsnowpack.lib
			PATHS ${SEARCH_PATH}
			DOC "Location of the libsnowpack, like c:/Program Files/Snowpack-2.0.0/lib/libsnowpack.lib"
			)
	ELSE(MSVC)
		FIND_LIBRARY(LIBSNOWPACK_LIBRARY
			NAMES libsnowpack.dll.a
			PATHS ${SEARCH_PATH}
			DOC "Location of the libsnowpack, like c:/Program Files/Snowpack-2.0.0/lib/libsnowpack.dll.a"
			)
	ENDIF(MSVC)
ELSE(WIN32)
	IF(APPLE)
		FIND_LIBRARY(LIBSNOWPACK_LIBRARY
		NAMES snowpack
		PATHS
			"/Applications/Snowpack/lib"
			ENV LD_LIBRARY_PATH
			ENV DYLD_FALLBACK_LIBRARY_PATH
			"~/usr/lib"
			"/usr/local/lib"
			"/usr/lib"
			"/opt/lib"
				"/Users/chris/Documents/PhD/code/snowpack-3.3.0/lib"
			${SRC_DIR}/snowpack/lib
		DOC "Location of the libsnowpack, like /usr/lib/libsnowpack.dylib"
		)
	ELSE(APPLE)
		FIND_LIBRARY(LIBSNOWPACK_LIBRARY
		NAMES snowpack
		PATHS
			ENV LD_LIBRARY_PATH
			"~/usr/lib"
			"/usr/local/lib"
			"/usr/lib"
			"/opt/lib"
				"/home/chris/Documents/PhD/code/snowpack-3.3.0/lib"
			${SRC_DIR}/snowpack/lib
				DOC "Location of the libsnowpack, like /usr/lib"
		)
	ENDIF(APPLE)
ENDIF(WIN32)

#build LIBSNOWPACK_ROOT so we can provide a hint for searching for the header file
GET_FILENAME_COMPONENT(snowpack_libs_root ${LIBSNOWPACK_LIBRARY} PATH)
IF(${CMAKE_VERSION} VERSION_GREATER "2.8.11")
	GET_FILENAME_COMPONENT(LIBSNOWPACK_ROOT ${snowpack_libs_root} DIRECTORY)
ELSE(${CMAKE_VERSION} VERSION_GREATER "2.8.11")
	SET(LIBSNOWPACK_ROOT "${snowpack_libs_root}/../")
	STRING(REPLACE  " " "\\ " LIBSNOWPACK_ROOT ${LIBSNOWPACK_ROOT})
ENDIF(${CMAKE_VERSION} VERSION_GREATER "2.8.11")

# locate main header file
FIND_PATH(LIBSNOWPACK_INCLUDE_DIR
  NAMES snowpack/libsnowpack.h
  #HINTS ${LIBSNOWPACK_ROOT}/include
  PATHS
	"${LIBSNOWPACK_ROOT}/include"
	"${LIBSNOWPACK_ROOT}"
	"~/usr/include"
	"/usr/local/include"
	"/usr/include"
	"/opt/include"
		"/Users/chris/Documents/PhD/code/snowpack-3.3.0/"
		"/home/chris/Documents/PhD/code/snowpack-3.3.0/"
		DOC "Location of the libsnowpack headers, like /usr/include"
)

# Set the include dir variables and the libraries and let libfind_process do the rest.
# NOTE: Singular variables for this library, plural for libraries this this lib depends on.
SET(LIBSNOWPACK_PROCESS_INCLUDES LIBSNOWPACK_INCLUDE_DIR)
SET(LIBSNOWPACK_PROCESS_LIBS LIBSNOWPACK_LIBRARY)
libfind_process(LIBSNOWPACK)
