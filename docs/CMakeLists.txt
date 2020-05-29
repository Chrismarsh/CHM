# add a target to generate API documentation with Doxygen
# generate this via 'make doc'
find_package(Doxygen)
if(DOXYGEN_FOUND)


    file(GLOB_RECURSE CHM_HEADERS ${CMAKE_SOURCE_DIR}/src/*.h)

    set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
    set(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/html/index.html)
    set(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in)
    set(DOXYFILE_OUT ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile)

    set( DOXYGEN_INPUT_DIR
#            ${CMAKE_SOURCE_DIR}/src/filters
#            ${CMAKE_SOURCE_DIR}/src/interpolation
#            ${CMAKE_SOURCE_DIR}/src/math
#            ${CMAKE_SOURCE_DIR}/src/mesh

            ${CMAKE_SOURCE_DIR}/src/modules
#            ${CMAKE_SOURCE_DIR}/src/modules/FSM
#            ${CMAKE_SOURCE_DIR}/src/modules/interp_met
#            ${CMAKE_SOURCE_DIR}/src/modules/snobal
#            ${CMAKE_SOURCE_DIR}/src/modules/snowpack
#
#            ${CMAKE_SOURCE_DIR}/src/physics
#            ${CMAKE_SOURCE_DIR}/src/timeseries
#            ${CMAKE_SOURCE_DIR}/src/utility
#            ${CMAKE_SOURCE_DIR}/src

    )

    string (REPLACE ";" " " DOXYGEN_INPUT_DIR "${DOXYGEN_INPUT_DIR}")

    configure_file(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)


#    file(MAKE_DIRECTORY ${DOXYGEN_OUTPUT_DIR})
#    add_custom_command(OUTPUT ${DOXYGEN_INDEX_FILE}
#            DEPENDS ${CHM_HEADERS}
#            COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
#            MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
#            COMMENT "Generating docs")
#
#    add_custom_target(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})

    find_program(SPHINX_EXECUTABLE
            NAMES sphinx-build
            DOC "Path to sphinx-build executable")

    set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
    set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)

    add_custom_target(Sphinx ALL
            COMMAND
            ${SPHINX_EXECUTABLE} -b html
            -Dbreathe_projects.CHM=${DOXYGEN_OUTPUT_DIR}/xml
#            -Dexhale_args.exhaleDoxygenStdin="INPUT = ${DOXYGEN_INPUT_DIR}"

            ${SPHINX_SOURCE} ${SPHINX_BUILD}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
            DEPENDS
            ${DOXYFILE_OUT}
#            MAIN_DEPENDENCY ${SPHINX_SOURCE}/conf.py
            COMMENT "Generating documentation with Sphinx")

endif(DOXYGEN_FOUND)