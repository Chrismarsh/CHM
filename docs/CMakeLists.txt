# add a target to generate API documentation with Doxygen
# generate this via 'make doc'
find_package(Doxygen)
if(DOXYGEN_FOUND)

    set(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doxygen)

    find_program(SPHINX_EXECUTABLE
            NAMES sphinx-build
            DOC "Path to sphinx-build executable")

    set(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR})
    set(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)

    add_custom_target(Sphinx ALL
            COMMAND
            ${SPHINX_EXECUTABLE} -b html
            -Dbreathe_projects.CHM=${DOXYGEN_OUTPUT_DIR}/xml

            ${SPHINX_SOURCE} ${SPHINX_BUILD}
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}

            COMMENT "Generating documentation with Sphinx")

endif(DOXYGEN_FOUND)