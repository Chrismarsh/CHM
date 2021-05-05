# - Find NetCDF
# Find the native NetCDF includes and library
#
#  NetCDF_INCLUDE_DIR  - user modifiable choice of where netcdf headers are
#  NetCDF_LIBRARY      - user modifiable choice of where netcdf libraries are
#
# Your package can require certain interfaces to be FOUND by setting these
#
#  NetCDF_CXX         - require the C++ interface and link the C++ library
#  NetCDF_F77         - require the F77 interface and link the fortran library
#  NetCDF_F90         - require the F90 interface and link the fortran library
#
# Or equivalently by calling FindNetCDF with a COMPONENTS argument containing one or
# more of "CXX;F77;F90".
#
# When interfaces are requested the user has access to interface specific hints:
#
#  NetCDF_${LANG}_INCLUDE_DIR - where to search for interface header files
#  NetCDF_${LANG}_LIBRARY     - where to search for interface libraries
#
# This module returns these variables for the rest of the project to use.
#
#  NetCDF_FOUND          - True if NetCDF found including required interfaces (see below)
#  NetCDF_LIBRARIES      - All netcdf related libraries.
#  NetCDF_INCLUDE_DIRS   - All directories to include.
#  NetCDF_HAS_INTERFACES - Whether requested interfaces were found or not.
#  NetCDF_${LANG}_INCLUDE_DIRS/NetCDF_${LANG}_LIBRARIES - C/C++/F70/F90 only interface
#
# Normal usage would be:
#  set (NetCDF_F90 "YES")
#  find_package (NetCDF REQUIRED)
#  target_link_libraries (uses_everthing ${NetCDF_LIBRARIES})
#  target_link_libraries (only_uses_f90 ${NetCDF_F90_LIBRARIES})

IF( DEFINED ENV{NetCDF_DIR} )
	SET( NetCDF_DIR "$ENV{NetCDF_DIR}" )
ENDIF()

# first look for the base c library
find_path(NetCDF_C_INCLUDE_DIR
		netcdf.h
    	PATHS "${NetCDF_DIR}/include")

find_library (NetCDF_C_LIBRARY
		NAMES netcdf
	    PATHS "${NetCDF_DIR}/lib")


find_package_handle_standard_args(NetCDF DEFAULT_MSG
		NetCDF_C_LIBRARY NetCDF_C_INCLUDE_DIR)

#start finding requested language components
set (NetCDF_libs "${NetCDF_C_LIBRARY}")
set (NetCDF_includes "${NetCDF_C_INCLUDE_DIR}")

get_filename_component (NetCDF_lib_dirs "${NetCDF_C_LIBRARY}" PATH)

if(NOT NetCDF_FIND_COMPONENTS)
	set(NetCDF_LANGUAGE_BINDINGS "C")
else()
	foreach(_component IN LISTS NetCDF_FIND_COMPONENTS)
		set(NetCDF_${_component} "YES")
	endforeach()
endif()


macro (NetCDF_check_interface lang header libs)
    if (NetCDF_${lang})
		#search starting from user modifiable cache var
		find_path (NetCDF_${lang}_INCLUDE_DIR
			NAMES ${header}

			HINTS "${NetCDF_C_INCLUDE_DIR}"
			HINTS "${NetCDF_${lang}_ROOT}/include"
			)

		# The C++ libraries can have different names depending on the platform. If we just search in a random order
		# we might pick up a system one that doesn't correspond to where we found the header. Here, we enforce that
		# the lib is ../ from the header files we found previously
		find_library (NetCDF_${lang}_LIBRARY 
				NAMES ${libs}

				NO_DEFAULT_PATH
				HINTS "${NetCDF_${lang}_INCLUDE_DIR}/../lib"
				HINTS "${NetCDF_${lang}_INCLUDE_DIR}/../lib64"
				PATHS "${NetCDF_${lang}_ROOT}/lib"
				PATHS "${NetCDF_${lang}_ROOT}/lib64"
			)


		find_package_handle_standard_args(NetCDF DEFAULT_MSG
				NetCDF_${lang}_INCLUDE_DIR NetCDF_${lang}_LIBRARY)

		mark_as_advanced (NetCDF_${lang}_INCLUDE_DIR NetCDF_${lang}_LIBRARY)

		#export to internal varS that rest of project can use directly
		set (NetCDF_${lang}_LIBRARIES ${NetCDF_${lang}_LIBRARY})
		set (NetCDF_${lang}_INCLUDE_DIRS ${NetCDF_${lang}_INCLUDE_DIR})

		if (NetCDF_${lang}_INCLUDE_DIR AND NetCDF_${lang}_LIBRARY)
			list (APPEND NetCDF_libs ${NetCDF_${lang}_LIBRARY})
			list (APPEND NetCDF_includes ${NetCDF_${lang}_INCLUDE_DIR})
		endif ()

    endif ()
endmacro ()

list (FIND NetCDF_FIND_COMPONENTS "CXX" _nextcomp)
if (_nextcomp GREATER -1)
    set (NetCDF_CXX 1)
endif ()
list (FIND NetCDF_FIND_COMPONENTS "F77" _nextcomp)
if (_nextcomp GREATER -1)
    set (NetCDF_F77 1)
endif ()
list (FIND NetCDF_FIND_COMPONENTS "F90" _nextcomp)
if (_nextcomp GREATER -1)
    set (NetCDF_F90 1)
endif ()

# try a few naming schemes, depends on the platform
list(APPEND lib_names "netcdf_c++4" "netcdf_c++" "netcdf-cxx4")
NetCDF_check_interface (CXX netcdf "${lib_names}")

NetCDF_check_interface (F77 netcdf.inc  netcdff)
NetCDF_check_interface (F90 netcdf.mod  netcdff)

#export accumulated results to internal varS that rest of project can depend on

set (NetCDF_LIBRARIES ${NetCDF_libs})
set (NetCDF_INCLUDE_DIRS ${NetCDF_includes})

find_package_handle_standard_args(NetCDF DEFAULT_MSG
		NetCDF_LIBRARIES NetCDF_INCLUDE_DIRS)

if(NetCDF_FOUND)
	mark_as_advanced (NetCDF_C_LIBRARY NetCDF_C_INCLUDE_DIR NetCDF_DIR)

	add_library(NetCDF::NetCDF INTERFACE IMPORTED)

	set_target_properties(NetCDF::NetCDF PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${NetCDF_INCLUDE_DIRS}")
	set_property(TARGET NetCDF::NetCDF PROPERTY INTERFACE_LINK_LIBRARIES "${NetCDF_LIBRARIES}")

	message(STATUS "NetCDF incl for all components -- ${NetCDF_INCLUDE_DIRS}")
	message(STATUS "NetCDF lib for all components -- ${NetCDF_LIBRARIES}")
endif()



