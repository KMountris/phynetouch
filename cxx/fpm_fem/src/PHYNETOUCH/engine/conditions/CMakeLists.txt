#----------------------------------------------
set(LIB_NAME ${LC_PROJECT_NAME}Conditions)
set(LIB_FOLDER conditions)

# Module header files.
set(HEADERS
    ${CMAKE_CURRENT_SOURCE_DIR}/body_load_bc.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/body_load_bc.tpp
    ${CMAKE_CURRENT_SOURCE_DIR}/bound_conds.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/dirichlet_bc.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/dirichlet_bc.tpp
    ${CMAKE_CURRENT_SOURCE_DIR}/initial_bc.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/load_curve_basic.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/load_curve_factory.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/load_curve_smooth.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/load_curve_step.hpp
    ${CMAKE_CURRENT_SOURCE_DIR}/prescribed_bc.hpp
)

# Module header files.
set(SOURCES
${CMAKE_CURRENT_SOURCE_DIR}/load_curve_smooth.cpp
${CMAKE_CURRENT_SOURCE_DIR}/load_curve_step.cpp
${CMAKE_CURRENT_SOURCE_DIR}/prescribed_bc.cpp
)

#-------- Build library --------
add_library(${LIB_NAME} "")
target_sources(${LIB_NAME} PRIVATE ${SOURCES})
target_link_libraries(${LIB_NAME} PUBLIC Eigen3::Eigen IMP::IMP)

include(GenerateExportHeader)
generate_export_header(${LIB_NAME}
    BASE_NAME "${LIB_NAME}"
    EXPORT_MACRO_NAME "${LIB_NAME}_EXPORT"
    EXPORT_FILE_NAME "${CMAKE_BINARY_DIR}/src/${PROJECT_NAME}/${LIB_FOLDER}/${LIB_NAME}Export.h"
    DEPRECATED_MACRO_NAME "${LIB_NAME}_DEPRECATED"
    NO_EXPORT_MACRO_NAME "${LIB_NAME}_NO_EXPORT"
    STATIC_DEFINE "${LIB_NAME}_STATIC_DEFINE"
    NO_DEPRECATED_MACRO_NAME "${LIB_NAME}_NO_DEPRECATED"
    DEFINE_NO_DEPRECATED
)

target_include_directories(${LIB_NAME} PUBLIC
    $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>/${PROJECT_NAME}
    $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/src>
)

set_target_properties(${LIB_NAME}
    PROPERTIES
        LINKER_LANGUAGE CXX
        POSITION_INDEPENDENT_CODE 1
        DEBUG_POSTFIX "_d"
        PUBLIC_HEADER "${HEADERS};${CMAKE_BINARY_DIR}/src/${PROJECT_NAME}/${LIB_FOLDER}/${LIB_NAME}Export.h"
        $<$<BUILD_SHARED_LIBS>:CXX_VISIBILITY_PRESET hidden>
        $<$<BUILD_SHARED_LIBS>:VISIBILITY_INLINES_HIDDEN 1>
        $<$<BUILD_SHARED_LIBS>:SOVERSION ${PROJECT_VERSION}>
        $<$<BUILD_SHARED_LIBS>:OUTPUT_NAME ${LIB_NAME}>
        $<$<BUILD_SHARED_LIBS>:MACOSX_RPATH ON>
        $<$<NOT:$<BUILD_SHARED_LIBS>>:ARCHIVE_OUTPUT_NAME ${LIB_NAME}>
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