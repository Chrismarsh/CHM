# - Find NetCDF
# Find the native NetCDF includes and library
#
#  NETCDF_INCLUDE_DIR  - user modifiable choice of where netcdf headers are
#  NETCDF_LIBRARY      - user modifiable choice of where netcdf libraries are
#
# Your package can require certain interfaces to be FOUND by setting these
#
#  NETCDF_CXX         - require the C++ interface and link the C++ library
#  NETCDF_F77         - require the F77 interface and link the fortran library
#  NETCDF_F90         - require the F90 interface and link the fortran library
#
# Or equivalently by calling FindNetCDF with a COMPONENTS argument containing one or
# more of "CXX;F77;F90".
#
# When interfaces are requested the user has access to interface specific hints:
#
#  NETCDF_${LANG}_INCLUDE_DIR - where to search for interface header files
#  NETCDF_${LANG}_LIBRARY     - where to search for interface libraries
#
# This module returns these variables for the rest of the project to use.
#
#  NETCDF_FOUND          - True if NetCDF found including required interfaces (see below)
#  NETCDF_LIBRARIES      - All netcdf related libraries.
#  NETCDF_INCLUDE_DIRS   - All directories to include.
#  NETCDF_HAS_INTERFACES - Whether requested interfaces were found or not.
#  NETCDF_${LANG}_INCLUDE_DIRS/NETCDF_${LANG}_LIBRARIES - C/C++/F70/F90 only interface
#
# Normal usage would be:
#  set (NETCDF_F90 "YES")
#  find_package (NetCDF REQUIRED)
#  target_link_libraries (uses_everthing ${NETCDF_LIBRARIES})
#  target_link_libraries (only_uses_f90 ${NETCDF_F90_LIBRARIES})

IF( DEFINED ENV{NetCDF_DIR} )
	SET( NetCDF_DIR "$ENV{NetCDF_DIR}" )
ENDIF()

find_path(NETCDF_INCLUDE_DIR
		include/netcdf.h
    	HINTS ${NETCDF_DIR})

find_library (NETCDF_LIBRARY
		NAMES netcdf
	    HINTS "${NETCDF_INCLUDE_DIR}/../lib")

#start finding requested language components
set (NetCDF_libs "")
set (NetCDF_includes "${NETCDF_INCLUDE_DIR}")

get_filename_component (NetCDF_lib_dirs "${NETCDF_LIBRARY}" PATH)
set (NETCDF_HAS_INTERFACES "YES") # will be set to NO if we're missing any interfaces

if(NOT NetCDF_FIND_COMPONENTS)
	set(NetCDF_LANGUAGE_BINDINGS "C")
else()
	foreach(_component IN LISTS NetCDF_FIND_COMPONENTS)
		set(NETCDF_${_component} "YES")
	endforeach()
endif()


macro (NetCDF_check_interface lang header libs)
    if (NETCDF_${lang})
		#search starting from user modifiable cache var
		find_path (NETCDF_${lang}_INCLUDE_DIR
			NAMES ${header}
			HINTS "${NETCDF_INCLUDE_DIR}"
			HINTS "${NETCDF_${lang}_ROOT}/include"
			)

		find_library (NETCDF_${lang}_LIBRARY NAMES ${libs}
			HINTS "${NetCDF_lib_dirs}"
			HINTS "${NETCDF_${lang}_ROOT}/lib"
			HINTS "${NETCDF_${lang}_ROOT}/lib64"
			)

		mark_as_advanced (NETCDF_${lang}_INCLUDE_DIR NETCDF_${lang}_LIBRARY)

		#export to internal varS that rest of project can use directly
		set (NETCDF_${lang}_LIBRARIES ${NETCDF_${lang}_LIBRARY})
		set (NETCDF_${lang}_INCLUDE_DIRS ${NETCDF_${lang}_INCLUDE_DIR})

		if (NETCDF_${lang}_INCLUDE_DIR AND NETCDF_${lang}_LIBRARY)
			list (APPEND NetCDF_libs ${NETCDF_${lang}_LIBRARY})
			list (APPEND NetCDF_includes ${NETCDF_${lang}_INCLUDE_DIR})
		else ()
			set (NETCDF_HAS_INTERFACES "NO")
		   # message (STATUS "Failed to find NetCDF interface for ${lang}")
		endif ()

		if(NOT NETCDF_HAS_INTERFACES)
			message (STATUS "Failed to find NetCDF interface for ${lang}")
		endif()
    endif ()
endmacro ()

list (FIND NetCDF_FIND_COMPONENTS "CXX" _nextcomp)
if (_nextcomp GREATER -1)
    set (NETCDF_CXX 1)
endif ()
list (FIND NetCDF_FIND_COMPONENTS "F77" _nextcomp)
if (_nextcomp GREATER -1)
    set (NETCDF_F77 1)
endif ()
list (FIND NetCDF_FIND_COMPONENTS "F90" _nextcomp)
if (_nextcomp GREATER -1)
    set (NETCDF_F90 1)
endif ()

# try a few naming schemes, depends on the platform
list(APPEND lib_names "netcdf_c++4" "netcdf_c++" "netcdf-cxx4")
NetCDF_check_interface (CXX netcdf "${lib_names}")

NetCDF_check_interface (F77 netcdf.inc  netcdff)
NetCDF_check_interface (F90 netcdf.mod  netcdff)

#export accumulated results to internal varS that rest of project can depend on
list (APPEND NetCDF_libs "${NETCDF_C_LIBRARIES}")
set (NETCDF_LIBRARIES ${NetCDF_libs})
set (NETCDF_INCLUDE_DIRS ${NetCDF_includes})

find_package_handle_standard_args(NetCDF DEFAULT_MSG
		NETCDF_LIBRARIES NETCDF_INCLUDE_DIRS NETCDF_HAS_INTERFACES)

if(NetCDF_FOUND)
	set (NETCDF_C_INCLUDE_DIRS ${NETCDF_INCLUDE_DIR})
	set (NETCDF_C_LIBRARIES ${NETCDF_LIBRARY})

	mark_as_advanced (NETCDF_LIBRARY)
	mark_as_advanced (NETCDF_INCLUDE_DIR)
	mark_as_advanced(NetCDF_DIR)

	add_library(NetCDF::NetCDF INTERFACE IMPORTED)

	set_target_properties(NetCDF::NetCDF PROPERTIES INTERFACE_INCLUDE_DIRECTORIES "${NETCDF_INCLUDE_DIRS}")
	set_property(TARGET NetCDF::NetCDF PROPERTY INTERFACE_LINK_LIBRARIES "${NETCDF_LIBRARIES}")

endif()



