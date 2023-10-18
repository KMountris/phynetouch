#----------------------------------------------
set(LIB_NAME ${LC_PROJECT_NAME}Physics)
set(LIB_FOLDER physics)

# Module header files.
set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/bioheat.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bioheat.tpp
    ${CMAKE_CURRENT_SOURCE_DIR}/contact_handle.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/contact_handle.tpp
    ${CMAKE_CURRENT_SOURCE_DIR}/deformation.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/deformation.tpp
    ${CMAKE_CURRENT_SOURCE_DIR}/electrical.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/electrical.tpp
)

#-------- Build INTERFACE library --------
add_library(${LIB_NAME} INTERFACE)
target_link_libraries(${LIB_NAME} INTERFACE Eigen3::Eigen)

target_include_directories(${LIB_NAME} INTERFACE
    $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>/${PROJECT_NAME}
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)
set_target_properties(${LIB_NAME} PROPERTIES
    PUBLIC_HEADER "${HEADERS}"
)
#--------------------------------------

install(TARGETS
        ${LIB_NAME}
    EXPORT
        ${LIB_NAME}Targets
    RUNTIME
        DESTINATION ${INSTALL_BINDIR}
        COMPONENT ${PROJECT_NAME}_Runtime
    LIBRARY
        DESTINATION ${INSTALL_LIBDIR}
        COMPONENT ${PROJECT_NAME}_Runtime
    ARCHIVE
        DESTINATION ${INSTALL_LIBDIR}
        COMPONENT ${PROJECT_NAME}_Development
    PUBLIC_HEADER
        DESTINATION ${INSTALL_INCLUDEDIR}/${PROJECT_NAME}/${LIB_FOLDER}
        COMPONENT ${PROJECT_NAME}_Development
)

install(EXPORT ${LIB_NAME}Targets
    NAMESPACE "${PROJECT_NAME}::"
    DESTINATION ${INSTALL_CMAKEDIR}
    COMPONENT ${PROJECT_NAME}_Development
)
