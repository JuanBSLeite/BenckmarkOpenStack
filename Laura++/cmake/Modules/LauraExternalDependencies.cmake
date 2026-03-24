# Laura++ needs ROOT, we just require that it is a version that exports CMake targets properly
# We may depend on other things in the future, e.g. JSON or YAML for easier configuration of particle properties etc.

if(DEFINED ENV{ROOTSYS})
    list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})
endif()

if(LAURA_BUILD_ROOFIT_TASK)
    find_package(ROOT 6.18 REQUIRED COMPONENTS EG RooFitCore RooFit)
else()
    find_package(ROOT 6.18 REQUIRED COMPONENTS EG)
endif()

message(STATUS "Laura++: Using ROOT installation from: ${ROOT_DIR}")
