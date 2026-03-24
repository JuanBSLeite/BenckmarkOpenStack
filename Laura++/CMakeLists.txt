# Top level CMakeLists.txt for Laura++

# Enforce an out-of-source build.
# Should be the first action in the top level CMakeLists.txt
if(${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_BINARY_DIR})
    message(STATUS "Laura++ requires an out-of-source build.")
    message(STATUS "Please remove these files from ${CMAKE_BINARY_DIR} first:")
    message(STATUS "    CMakeCache.txt")
    message(STATUS "    CMakeFiles")
    message(STATUS "Once these files are removed, create a separate directory")
    message(STATUS "and run CMake from there, pointing it to:")
    message(STATUS "    ${CMAKE_SOURCE_DIR}")
    message(FATAL_ERROR "in-source build detected")
endif()

# Also require a minimum version of CMake
cmake_minimum_required(VERSION 3.13.0)

# Project setup
project(Laura++ VERSION 3.7.0
    DESCRIPTION "A maximum likelihood fitting package for performing Dalitz-plot analysis."
    HOMEPAGE_URL "https://laura.hepforge.org"
    )

# Option to enable/disable compilation of the RooFit task class
option(LAURA_BUILD_ROOFIT_TASK "Enable/disable compilation of RooFit task class" OFF)

# Option to enable/disable building of the Doxygen documentation
option(LAURA_BUILD_DOCS "Enable/disable building of Doxygen documentation" OFF)

# Options to enable/disable compilation of the example/test executables
option(LAURA_BUILD_EXAMPLES "Enable/disable compilation of example executables" OFF)
option(LAURA_BUILD_TESTS "Enable/disable compilation of test executables" OFF)

message(STATUS "Laura++: Optional compilation of RooFit task class  LAURA_BUILD_ROOFIT_TASK  ${LAURA_BUILD_ROOFIT_TASK}")
message(STATUS "Laura++: Optional building of Doxygen documentation LAURA_BUILD_DOCS         ${LAURA_BUILD_DOCS}")
message(STATUS "Laura++: Optional building of example executables   LAURA_BUILD_EXAMPLES     ${LAURA_BUILD_EXAMPLES}")
message(STATUS "Laura++: Optional building of test executables      LAURA_BUILD_TESTS        ${LAURA_BUILD_TESTS}")

# Prepend this project's custom module path(s) to CMAKE_MODULE_PATH
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})

# Include needed modules and perform any custom setup
# Install paths
include(GNUInstallDirs)
# Compilation/linking flags and related settings
include(LauraCompilerFlags)
# External dependencies
include(LauraExternalDependencies)

# Build and install the library
add_subdirectory(src)

# Build the ROOT dictionary and install the headers
add_subdirectory(inc)

# Create the documentation
if(LAURA_BUILD_DOCS)
    add_subdirectory(doc)
endif()

# Build and install the executables in the examples directory
if(LAURA_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(LAURA_BUILD_TESTS)
    add_subdirectory(test)
endif()

# Generate CMake config files, which can be used by other projects
include(CMakePackageConfigHelpers)

set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR})
set(LIB_INSTALL_DIR ${CMAKE_INSTALL_LIBDIR})
set(DATA_INSTALL_DIR ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME})

configure_package_config_file(cmake/Templates/LauraConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/LauraConfig.cmake
    INSTALL_DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/cmake
    PATH_VARS INCLUDE_INSTALL_DIR LIB_INSTALL_DIR DATA_INSTALL_DIR
    )

write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/LauraConfigVersion.cmake
    COMPATIBILITY SameMajorVersion
    )

install(FILES ${CMAKE_CURRENT_BINARY_DIR}/LauraConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/LauraConfigVersion.cmake
    DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/cmake
    )

install(
    EXPORT "LauraTargets"
    NAMESPACE "Laura::"
    DESTINATION ${CMAKE_INSTALL_DATADIR}/${CMAKE_PROJECT_NAME}/cmake
    )
