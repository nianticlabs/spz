cmake_minimum_required(VERSION 3.10)

project(spz
  DESCRIPTION "A 3D Gaussians format"
  LANGUAGES C CXX
  VERSION 1.1.0)

include(GNUInstallDirs)

if(CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
  set(SPZ_BUILD_WASM ON)
  set(SPZ_USE_EMSCRIPTEN_ZLIB ON)
endif()

# zlib is required to build the project
if(NOT SPZ_USE_EMSCRIPTEN_ZLIB)
  find_package(ZLIB REQUIRED)
endif()

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

if(NOT SPZ_USE_EMSCRIPTEN_ZLIB)
  target_link_libraries(spz PRIVATE ZLIB::ZLIB)
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
