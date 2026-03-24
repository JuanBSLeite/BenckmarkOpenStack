# Use glob to find the sources
file(GLOB LAURA_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/*.cc)
if (NOT LAURA_BUILD_ROOFIT_TASK)
    list(REMOVE_ITEM LAURA_SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/LauRooFitTask.cc)
endif()

# Build the shared library
add_library(Laura++ SHARED ${LAURA_SOURCES})
set_target_properties(Laura++ PROPERTIES OUTPUT_NAME Laura++)
set_target_properties(Laura++ PROPERTIES VERSION ${CMAKE_PROJECT_VERSION} SOVERSION ${CMAKE_PROJECT_VERSION_MAJOR})
set_target_properties(Laura++ PROPERTIES LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/${CMAKE_INSTALL_LIBDIR})
target_include_directories(Laura++ PUBLIC $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/inc> $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/${CMAKE_PROJECT_NAME}>)
target_link_libraries(Laura++ PUBLIC ROOT::Core ROOT::RIO ROOT::Hist ROOT::Matrix ROOT::Physics ROOT::Minuit ROOT::EG ROOT::Tree)
if (LAURA_BUILD_ROOFIT_TASK)
    target_link_libraries(Laura++ PUBLIC ROOT::RooFit ROOT::RooFitCore)
endif()

# Install the libraries
install(
    TARGETS Laura++
    EXPORT "LauraTargets"
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )

