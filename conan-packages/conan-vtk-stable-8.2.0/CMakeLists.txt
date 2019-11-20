cmake_minimum_required(VERSION 2.8.11)
project(cmake_wrapper)

# Instead of the following two lines include some of the
# Conan logic directly:
# include(${CMAKE_BINARY_DIR}/../conanbuildinfo.cmake)
# conan_basic_setup()

# conan_set_fpic():
if(DEFINED CONAN_CMAKE_POSITION_INDEPENDENT_CODE)
    message(STATUS "Conan: Adjusting fPIC flag (${CONAN_CMAKE_POSITION_INDEPENDENT_CODE})")
    set(CMAKE_POSITION_INDEPENDENT_CODE ${CONAN_CMAKE_POSITION_INDEPENDENT_CODE})
endif()

# conan_set_vs_runtime():
if(CONAN_LINK_RUNTIME)
    foreach(flag CMAKE_C_FLAGS_RELEASE CMAKE_CXX_FLAGS_RELEASE
                 CMAKE_C_FLAGS_RELWITHDEBINFO CMAKE_CXX_FLAGS_RELWITHDEBINFO
                 CMAKE_C_FLAGS_MINSIZEREL CMAKE_CXX_FLAGS_MINSIZEREL)
        if(DEFINED ${flag})
            string(REPLACE "/MD" ${CONAN_LINK_RUNTIME} ${flag} "${${flag}}")
        endif()
    endforeach()
    foreach(flag CMAKE_C_FLAGS_DEBUG CMAKE_CXX_FLAGS_DEBUG)
        if(DEFINED ${flag})
            string(REPLACE "/MDd" ${CONAN_LINK_RUNTIME} ${flag} "${${flag}}")
        endif()
    endforeach()
endif()

# conan_set_std():
set(ignorevar "${CONAN_STD_CXX_FLAG}${CONAN_CMAKE_CXX_STANDARD}${CONAN_CMAKE_CXX_EXTENSIONS}")
if (CMAKE_VERSION VERSION_LESS "3.1" OR
    (CMAKE_VERSION VERSION_LESS "3.12" AND ("${CONAN_CMAKE_CXX_STANDARD}" STREQUAL "20" OR "${CONAN_CMAKE_CXX_STANDARD}" STREQUAL "gnu20")))
if(CONAN_STD_CXX_FLAG)
    message(STATUS "Conan setting CXX_FLAGS flags: ${CONAN_STD_CXX_FLAG}")
    set(CMAKE_CXX_FLAGS "${CONAN_STD_CXX_FLAG} ${CMAKE_CXX_FLAGS}")
endif()
else()
    if(CONAN_CMAKE_CXX_STANDARD)
        message(STATUS "Conan setting CPP STANDARD: ${CONAN_CMAKE_CXX_STANDARD} WITH EXTENSIONS ${CONAN_CMAKE_CXX_EXTENSIONS}")
        set(CMAKE_CXX_STANDARD ${CONAN_CMAKE_CXX_STANDARD})
        set(CMAKE_CXX_EXTENSIONS ${CONAN_CMAKE_CXX_EXTENSIONS})
    endif()
endif()

# conan_set_libcxx():
if(DEFINED CONAN_LIBCXX)
    message(STATUS "Conan: C++ stdlib: ${CONAN_LIBCXX}")
    if(CONAN_COMPILER STREQUAL "clang" OR CONAN_COMPILER STREQUAL "apple-clang")
        if(CONAN_LIBCXX STREQUAL "libstdc++" OR CONAN_LIBCXX STREQUAL "libstdc++11" )
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libstdc++")
        elseif(CONAN_LIBCXX STREQUAL "libc++")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
        endif()
    endif()
    if(CONAN_COMPILER STREQUAL "sun-cc")
        if(CONAN_LIBCXX STREQUAL "libCstd")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -library=Cstd")
        elseif(CONAN_LIBCXX STREQUAL "libstdcxx")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -library=stdcxx4")
        elseif(CONAN_LIBCXX STREQUAL "libstlport")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -library=stlport4")
        elseif(CONAN_LIBCXX STREQUAL "libstdc++")
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -library=stdcpp")
        endif()
    endif()
    if(CONAN_LIBCXX STREQUAL "libstdc++11")
        add_definitions(-D_GLIBCXX_USE_CXX11_ABI=1)
    elseif(CONAN_LIBCXX STREQUAL "libstdc++")
        add_definitions(-D_GLIBCXX_USE_CXX11_ABI=0)
    endif()
endif()

if(MSVC)
    # Remove pdb-related flags
    string(REGEX REPLACE "/Z[iI7]" ""
        CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
    # Embed debug info into .lib-files
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_RELEASE} /Z7")
endif()

add_subdirectory(source_subfolder)
