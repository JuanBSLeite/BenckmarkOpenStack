
# Set the build type (if not already specified)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Laura++: Setting build type to 'Release' as none was specified")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build, options are: Release, MinSizeRel, Debug, RelWithDebInfo" FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
elseif(CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    message(STATUS "Laura++: Build type '${CMAKE_BUILD_TYPE}'")
endif()

# Set the warning/optimise/debug flags for each build type
if( ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" OR ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang" )
    message(STATUS "Laura++: Customising warning/optimise/debug flags for each build type")

    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsigned-char -Wall -Wextra -Wshadow -Woverloaded-virtual")

    if( ${CMAKE_CXX_COMPILER_ID} STREQUAL "GNU" )
        set(CMAKE_CXX_FLAGS_DEBUG          "-Og -g3")
        set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
        set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")
        set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g3")
    elseif( ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang" )
        set(CMAKE_CXX_FLAGS_DEBUG          "-O0 -g")
        set(CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
        set(CMAKE_CXX_FLAGS_RELEASE        "-O3 -DNDEBUG")
        set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")
    endif()
else()
    message(WARNING "Laura++: No customisation of warning/optimise/debug flags implemented for compiler: ${CMAKE_CXX_COMPILER_ID}")
endif()

# Make sure our project's include directories always come first
set(CMAKE_INCLUDE_DIRECTORIES_PROJECT_BEFORE ON)

# Control verbosity of the build
set(CMAKE_VERBOSE_MAKEFILE OFF CACHE BOOL "Control verbosity of generated Makefiles")

# Control verbosity of the install - suppresses 'up-to-date' messages
set(CMAKE_INSTALL_MESSAGE LAZY)

# C++ standard settings
set(CMAKE_CXX_EXTENSIONS OFF)
set(LAURA_CXX_STANDARD 17)
if(DEFINED ENV{CMAKE_CXX_STANDARD})
    if($ENV{CMAKE_CXX_STANDARD} GREATER ${LAURA_CXX_STANDARD})
        set(LAURA_CXX_STANDARD $ENV{CMAKE_CXX_STANDARD})
    elseif($ENV{CMAKE_CXX_STANDARD} LESS ${LAURA_CXX_STANDARD})
        message(WARNING "Laura++: Version of C++ standard requested via CMAKE_CXX_STANDARD environment variable, C++$ENV{CMAKE_CXX_STANDARD}, is below required minimum of C++${LAURA_CXX_STANDARD}, not changing it")
    endif()
endif()
set(CMAKE_CXX_STANDARD ${LAURA_CXX_STANDARD} CACHE STRING "C++ standard")
set(CMAKE_CXX_STANDARD_REQUIRED ON)
message(STATUS "Laura++: Using C++${CMAKE_CXX_STANDARD} standard")

# Special linker flags for MacOSX
if (APPLE)
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -single_module -undefined dynamic_lookup")
endif()

# RPATH handling
set(CMAKE_MACOSX_RPATH TRUE)
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

