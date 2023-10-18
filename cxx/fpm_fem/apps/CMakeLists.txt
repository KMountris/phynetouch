####
# Build PhyNeTouch apps

# Options for building PhyNeTouch apps
option(BUILD_APPS_WITH_EXTERNAL_JSON "Use an external build of the NLOHMANN JSON library" OFF)
option(BUILD_APPS_WITH_EXTERNAL_TERMCOLOR "Use an external build of the TERMCOLOR library" OFF)

add_subdirectory(${CMAKE_CURRENT_SOURCE_DIR}/tools)

# Third party libraries for building PhyNeTouch apps
add_subdirectory(${CMAKE_SOURCE_DIR}/thirdparty/nlohmann ${CMAKE_BINARY_DIR}/thirdparty/nlohmann)
add_subdirectory(${CMAKE_SOURCE_DIR}/thirdparty/termcolor ${CMAKE_BINARY_DIR}/thirdparty/termcolor)


# Build PhyNeTouchSim app
add_executable(PhyNeTouchSim phynetouch_sim.cpp)
target_compile_definitions(PhyNeTouchSim PRIVATE PHYNETOUCH_VERSION="${PROJECT_VERSION}")
target_link_libraries(PhyNeTouchSim PRIVATE ${PROJECT_NAME} Tools nlohmann_json::nlohmann_json termcolor::termcolor pthread IMP::IMP)
target_include_directories(PhyNeTouchSim PUBLIC
    $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/${INSTALL_INCLUDEDIR}>/${PROJECT_NAME}
    $<INSTALL_INTERFACE:${INSTALL_INCLUDEDIR}>/${PROJECT_NAME}
)

# Install apps in targets
install(TARGETS
    PhyNeTouchSim
    RUNTIME
        DESTINATION ${INSTALL_BINDIR}
        COMPONENT ${PROJECT_NAME}_Runtime
)