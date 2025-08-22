cmake_minimum_required(VERSION 3.10)

project(spz
  DESCRIPTION "A 3D Gaussians format"
  LANGUAGES C CXX
  VERSION 1.1.0)

include(GNUInstallDirs)

option(SPZ_BUILD_PYTHON_BINDINGS "Build Python bindings using pybind11" OFF)

if(CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
  set(SPZ_BUILD_PYTHON_BINDINGS OFF CACHE BOOL "Build Python bindings using pybind11" FORCE)
  set(SPZ_BUILD_WASM ON)
  set(SPZ_USE_EMSCRIPTEN_ZLIB ON)
endif()

set(SPZ_SHOULD_INSTALL TRUE)

# zlib is required to build the project
if(NOT SPZ_USE_EMSCRIPTEN_ZLIB)
  find_package(ZLIB QUIET)

  if (NOT ZLIB_FOUND)
    message(STATUS "pybind11 not found, using FetchContent to download it.")
    include(FetchContent)
    FetchContent_Declare(
      zlib
      URL https://zlib.net/zlib-1.3.1.tar.gz
      URL_HASH SHA256=9a93b2b7dfdac77ceba5a558a580e74667dd6fede4585b91eefb60f03b72df23
    )
    FetchContent_MakeAvailable(zlib)
    set(zlib_linkage_target zlib)

    # Skipping install of SPZ target because zlib was fetched from the source and non-exportable
    set(SPZ_SHOULD_INSTALL FALSE)
  else()
    set(zlib_linkage_target ZLIB::ZLIB)
  endif()
endif()

if(SPZ_BUILD_PYTHON_BINDINGS)
  find_package(Python3 COMPONENTS Interpreter Development.Module REQUIRED)

  # Try to find pybind11, otherwise fetch it
  find_package(pybind11 CONFIG QUIET)
  if(NOT pybind11_FOUND)
    message(STATUS "pybind11 not found, using FetchContent to download it.")
    include(FetchContent)
    FetchContent_Declare(
      pybind11
      GIT_REPOSITORY https://github.com/pybind/pybind11.git
      GIT_TAG        v2.13.1
    )
    FetchContent_MakeAvailable(pybind11)
  else()
    message(STATUS "Using system-installed pybind11")
  endif()
endif()

set(spz_sources
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cc/load-spz.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cc/splat-c-types.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cc/splat-types.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cc/splat-extensions.cc"
)

set(spz_headers
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cc/load-spz.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cc/splat-c-types.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cc/splat-types.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cc/splat-extensions.h"
)

# create the library and configure it
add_library(spz ${spz_sources})
add_library(spz::spz ALIAS spz)

if(NOT SPZ_USE_EMSCRIPTEN_ZLIB)
  target_link_libraries(spz PRIVATE ${zlib_linkage_target})
else()
  target_compile_options(spz PRIVATE --use-port=zlib)
  target_link_options(spz PRIVATE --use-port=zlib)
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

if(SPZ_BUILD_WASM)
  file(MAKE_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/dist")

  set(spz_binding_sources
    "${CMAKE_CURRENT_SOURCE_DIR}/src/emscripten/spz-bindings.cc"
  )

  add_executable(spz-wasm ${spz_binding_sources})
  target_link_libraries(spz-wasm PRIVATE spz)

  target_link_options(spz-wasm PRIVATE
    "-O3"
    "-sUSE_ZLIB=1"
    "--bind"
    "-sMODULARIZE=1"
    "-sEXPORT_NAME=createSpzModule"
    "-sALLOW_MEMORY_GROWTH=1"
    "-sENVIRONMENT=web"
    "-sEXPORT_ES6=1"
    "-sEXPORTED_RUNTIME_METHODS=[\"ccall\",\"cwrap\"]"
  )

  set_target_properties(spz-wasm PROPERTIES
    OUTPUT_NAME "spz"
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/dist"
  )

  # Copy the TypeScript declaration file
  add_custom_command(
    OUTPUT "${CMAKE_CURRENT_SOURCE_DIR}/dist/spz.d.ts"
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
            "${CMAKE_CURRENT_SOURCE_DIR}/src/emscripten/spz.d.ts"
            "${CMAKE_CURRENT_SOURCE_DIR}/dist/spz.d.ts"
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/src/emscripten/spz.d.ts"
  )

  add_custom_target(copy_spz_dts ALL
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/dist/spz.d.ts"
  )
else()
  set_target_properties(spz PROPERTIES POSITION_INDEPENDENT_CODE ON)
endif()

if(SPZ_BUILD_PYTHON_BINDINGS)
  file(GLOB_RECURSE SPZ_PYBIND_SOURCES
    "${CMAKE_CURRENT_SOURCE_DIR}/src/pybindings/*.cpp"
  )

  pybind11_add_module(spz_bindings MODULE ${SPZ_PYBIND_SOURCES})

  target_link_libraries(spz_bindings
    PRIVATE spz
    PRIVATE pybind11::module
  )

  set_target_properties(spz_bindings PROPERTIES
    CXX_STANDARD 17
    CXX_STANDARD_REQUIRED ON
    POSITION_INDEPENDENT_CODE ON
  )

  if(WIN32)
    # Windows (.pyd)
    set_target_properties(spz_bindings PROPERTIES SUFFIX ".pyd")
  else()
    # Unix (.so)
    set_target_properties(spz_bindings PROPERTIES SUFFIX ".so")
  endif()

  if(APPLE)
    set_target_properties(spz_bindings PROPERTIES
      XCODE_ATTRIBUTE_PRODUCT_BUNDLE_IDENTIFIER "com.nianticlabs.spz.spz_bindings"
    )
  endif()
endif()

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

if (SPZ_SHOULD_INSTALL)
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
endif()

if(SPZ_BUILD_PYTHON_BINDINGS)
  install(TARGETS spz_bindings
    LIBRARY DESTINATION "${Python3_SITEARCH}"
    RUNTIME DESTINATION "${Python3_SITEARCH}"  # For Windows .pyd
  )
endif()
