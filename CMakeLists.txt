cmake_minimum_required(VERSION 3.10)

project(spz
  DESCRIPTION "A 3D Gaussians format"
  LANGUAGES C CXX
  VERSION 1.1.0)

include(GNUInstallDirs)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)


# zlib is required to build the project
find_package(ZLIB REQUIRED)

set(spz_sources
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cc/load-spz.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cc/splat-c-types.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cc/splat-types.cc"
)

set(spz_headers
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cc/load-spz.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cc/splat-c-types.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cc/splat-types.h"
)

# create the library and configure it
add_library(spz ${spz_sources})
add_library(spz::spz ALIAS spz)

target_link_libraries(spz PRIVATE ZLIB::ZLIB)
if(ANDROID)
  target_link_libraries(spz PRIVATE log)
endif()

target_include_directories(spz
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/cc>
  INTERFACE $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

set_target_properties(spz PROPERTIES
  PUBLIC_HEADER "${spz_headers}"
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON
)

# Installation

## Install spzConfig.cmake spz::spz target can be found when calling find_package(spz)
include(CMakePackageConfigHelpers)
configure_package_config_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/cmake/spzConfig.cmake.in"
  "${CMAKE_BINARY_DIR}/cmake/spzConfig.cmake"
  INSTALL_DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/spz")
write_basic_package_version_file(
  "${CMAKE_BINARY_DIR}/cmake/spzConfigVersion.cmake"
  VERSION "${spz_VERSION}"
  COMPATIBILITY SameMajorVersion)
install(
  FILES
    "${CMAKE_BINARY_DIR}/cmake/spzConfig.cmake"
    "${CMAKE_BINARY_DIR}/cmake/spzConfigVersion.cmake"
  DESTINATION
    "${CMAKE_INSTALL_LIBDIR}/cmake/spz"
)

install(TARGETS spz
  EXPORT spzTargets
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT spzTargets
  NAMESPACE spz::
  DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/spz"
)

# --- Command Line Tools Option ---
option(BUILD_TOOLS "Build command line tools" ON)

if(BUILD_TOOLS)
  # CLI Tools
  add_executable(ply_to_spz cli_tools/src/ply_to_spz.cpp)
  target_link_libraries(ply_to_spz PRIVATE spz)
  target_include_directories(ply_to_spz PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
  set_target_properties(ply_to_spz PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
  )

  add_executable(spz_to_ply cli_tools/src/spz_to_ply.cpp)
  target_link_libraries(spz_to_ply PRIVATE spz)
  target_include_directories(spz_to_ply PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
  set_target_properties(spz_to_ply PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
  )

  add_executable(spz_info cli_tools/src/spz_info.cpp)
  target_link_libraries(spz_info PRIVATE spz)
  target_include_directories(spz_info PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
  set_target_properties(spz_info PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
  )

  install(TARGETS ply_to_spz spz_to_ply spz_info
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )
endif()

# --- Python Bindings Option ---
option(BUILD_PYTHON_BINDINGS "Build Python bindings using nanobind" OFF)

# --- Python Bindings (nanobind) ---
if(BUILD_PYTHON_BINDINGS)
  # Find Python with the Development component, which is required by nanobind.
  find_package(Python 3.8
    REQUIRED COMPONENTS Interpreter Development.Module
    OPTIONAL_COMPONENTS Development.SABIModule)

  # Detect the installed nanobind package and import it into CMake
  execute_process(
    COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
    OUTPUT_STRIP_TRAILING_WHITESPACE OUTPUT_VARIABLE nanobind_DIR)
  find_package(nanobind CONFIG REQUIRED)


  nanobind_add_module(
    # Name of the extension
    spz_py

    # Target the stable ABI for Python 3.12+, which reduces
    # the number of binary wheels that must be built. This
    # does nothing on older Python versions
    STABLE_ABI

    # Build libnanobind statically and merge it into the
    # extension (which itself remains a shared library)
    NB_STATIC

    # Source code
    src/python/spz/spz.cc
  )

  # Rename the output file to "spz" (no prefix/lib, correct extension)
  # This is necessary for the Python module to be recognized correctly
  # by Python's import system.
  set_target_properties(spz_py PROPERTIES
      OUTPUT_NAME spz
  )

  # Add the project root to the include path for the python module
  target_include_directories(spz_py PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

  # Link the python module against the spz library
  target_link_libraries(spz_py PRIVATE spz::spz)

  # Install directive for scikit-build-core
  install(TARGETS spz_py LIBRARY DESTINATION spz)
endif()
