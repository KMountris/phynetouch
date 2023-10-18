#-----------------------------------------
# Find Doxygen and Sphinx installation.
FIND_PACKAGE(Doxygen REQUIRED)
FIND_PACKAGE(Sphinx REQUIRED)

# Find all the public headers
GET_TARGET_PROPERTY(${PROJECT_NAME}_PUBLIC_HEADER_DIR ${PROJECT_NAME} INTERFACE_INCLUDE_DIRECTORIES)
FILE(GLOB_RECURSE ${PROJECT_NAME}_PUBLIC_HEADERS ${${PROJECT_NAME}_PUBLIC_HEADER_DIR}/*.hpp)

# Set up folder alias for Doxygen.
SET(DOXYGEN_INPUT_DIR ${PROJECT_SOURCE_DIR}/phynetouch/docs)
SET(DOXYGEN_OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}/doxygen)
SET(DOXYGEN_INDEX_FILE ${DOXYGEN_OUTPUT_DIR}/html/index.html)
SET(DOXYFILE_IN ${CMAKE_CURRENT_SOURCE_DIR}/doxyfile.in)
SET(DOXYFILE_OUT ${CMAKE_CURRENT_BINARY_DIR}/doxygen/doxyfile)

# Set up folder alias for Sphinx.
SET(SPHINX_SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/source)
SET(SPHINX_BUILD ${CMAKE_CURRENT_BINARY_DIR}/sphinx)
SET(SPHINX_CONF_IN ${SPHINX_SOURCE}/conf.py.in)
SET(SPHINX_INDEXRST_IN ${SPHINX_SOURCE}/index.rst.in)
SET(SPHINX_CONF_OUT ${SPHINX_BUILD}/conf.py)
SET(SPHINX_INDEXRST_OUT ${SPHINX_BUILD}/index.rst)
SET(SPHINX_INDEX_FILE ${SPHINX_BUILD}/index.html)

# Replace variables inside @@ with the current values.
CONFIGURE_FILE(${DOXYFILE_IN} ${DOXYFILE_OUT} @ONLY)
CONFIGURE_FILE(${SPHINX_CONF_IN} ${SPHINX_CONF_OUT} @ONLY)
CONFIGURE_FILE(${SPHINX_INDEXRST_IN} ${SPHINX_INDEXRST_OUT} @ONLY)

# Make custom command and target for Doxygen.
ADD_CUSTOM_COMMAND(OUTPUT ${DOXYGEN_INDEX_FILE}
                   DEPENDS ${${PROJECT_NAME}_PUBLIC_HEADERS}
                   COMMAND ${DOXYGEN_EXECUTABLE} ${DOXYFILE_OUT}
                   MAIN_DEPENDENCY ${DOXYFILE_OUT} ${DOXYFILE_IN}
                   COMMENT "Generating the ${PROJECT_NAME} Doxygen documentation..."
)
ADD_CUSTOM_TARGET(Doxygen ALL DEPENDS ${DOXYGEN_INDEX_FILE})

# Make custom command and target for Sphinx.
ADD_CUSTOM_COMMAND(OUTPUT ${SPHINX_INDEX_FILE}
                   COMMAND
                     ${SPHINX_EXECUTABLE} -b html
                     # Tell Breathe where to find the Doxygen output
                     -Dbreathe_projects.${PROJECT_NAME}=${DOXYGEN_OUTPUT_DIR}/xml
                   ${SPHINX_SOURCE} ${SPHINX_BUILD}
                   WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                   DEPENDS
                   # Other docs files you want to track should go here (or in some variable)
                   ${SPHINX_INDEXRST_OUT}
                   ${DOXYGEN_INDEX_FILE}
                   MAIN_DEPENDENCY ${SPHINX_CONF_OUT}
                   COMMENT "Generating the ${PROJECT_NAME} Sphinx documentation..."
)
ADD_CUSTOM_TARGET(Sphinx ALL DEPENDS ${SPHINX_INDEX_FILE})
