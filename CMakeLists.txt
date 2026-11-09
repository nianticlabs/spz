cmake_minimum_required(VERSION 3.14...3.31)

project(spz
  DESCRIPTION "A 3D Gaussians format"
  LANGUAGES C CXX
  VERSION 1.1.0)

include(GNUInstallDirs)

# --- Python Bindings Option ---
option(SPZ_BUILD_PYTHON_BINDINGS "Build Python bindings using nanobind" OFF)

# --- Command Line Tools Option ---
option(SPZ_BUILD_TOOLS "Build command line tools" ON)

# --- Extensions Option ---
option(SPZ_BUILD_EXTENSIONS "Build SPZ extensions support" OFF)

# --- WASM output directory ---
set(SPZ_WASM_OUTPUT_DIR "${PROJECT_BINARY_DIR}/dist" CACHE PATH
  "Directory that receives spz.js and its TypeScript declarations (WASM builds only)")

# --- Emscripten guard ---
# The wasm presets set this so that forgetting `emcmake` fails loudly instead of
# silently producing a native build in build/wasm*.
option(SPZ_REQUIRE_EMSCRIPTEN "Fail configuration unless the Emscripten toolchain is active" OFF)
if(SPZ_REQUIRE_EMSCRIPTEN AND NOT CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
  message(FATAL_ERROR
    "SPZ_REQUIRE_EMSCRIPTEN is ON but the Emscripten toolchain is not active. "
    "Configure with `emcmake cmake --preset <wasm-preset>` or use scripts/build.sh. "
    "If this build directory was first configured without emcmake, add --fresh.")
endif()

if(CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
  set(SPZ_BUILD_PYTHON_BINDINGS OFF CACHE BOOL "Build Python bindings using nanobind" FORCE)
  set(SPZ_BUILD_TOOLS OFF CACHE BOOL "Build command line tools" FORCE)
  set(SPZ_BUILD_WASM ON)
  set(SPZ_USE_EMSCRIPTEN_ZLIB ON)
endif()

set(SPZ_SHOULD_INSTALL TRUE)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)


# zlib is required to build the project
if(NOT SPZ_USE_EMSCRIPTEN_ZLIB)
  find_package(ZLIB QUIET)

  if (NOT ZLIB_FOUND)
    message(STATUS "ZLIB not found, using FetchContent to download it.")
    include(FetchContent)
    FetchContent_Declare(
      zlib
      URL https://github.com/madler/zlib/releases/download/v1.3.2/zlib-1.3.2.tar.gz
      URL_HASH SHA256=bb329a0a2cd0274d05519d61c667c062e06990d72e125ee2dfa8de64f0119d16
    )
    FetchContent_MakeAvailable(zlib)
    set(zlib_linkage_target zlib)

    # Skipping install of SPZ target because zlib was fetched from the source and non-exportable
    set(SPZ_SHOULD_INSTALL FALSE)
  else()
    set(zlib_linkage_target ZLIB::ZLIB)
  endif()
endif()

# zstd is required for ZSTD compression/decompression
find_package(zstd QUIET)
if(NOT zstd_FOUND)
  message(STATUS "zstd not found, using FetchContent to download it.")
  include(FetchContent)
  FetchContent_Declare(
    zstd
    URL https://github.com/facebook/zstd/releases/download/v1.5.6/zstd-1.5.6.tar.gz
    URL_HASH SHA256=8c29e06cf42aacc1eafc4077ae2ec6c6fcb96a626157e0593d5e82a34fd403c1
    SOURCE_SUBDIR build/cmake
  )
  set(ZSTD_BUILD_SHARED OFF CACHE BOOL "" FORCE)
  set(ZSTD_BUILD_PROGRAMS OFF CACHE BOOL "" FORCE)
  set(ZSTD_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  if(CMAKE_SYSTEM_NAME STREQUAL "Emscripten")
    # The WASM build is single-threaded (see runParallel in zstd-batch.cc), and
    # zstd's find_package(Threads) fails under newer Emscripten + CMake 4.x.
    set(ZSTD_MULTITHREAD_SUPPORT OFF CACHE BOOL "" FORCE)
  endif()
  FetchContent_MakeAvailable(zstd)
  add_library(zstd::zstd ALIAS libzstd_static)
  set(SPZ_SHOULD_INSTALL FALSE)
else()
  if(TARGET zstd::libzstd_static)
    add_library(zstd::zstd ALIAS zstd::libzstd_static)
  elseif(NOT TARGET zstd::zstd)
    add_library(zstd::zstd ALIAS zstd::libzstd_shared)
  endif()
endif()

set(spz_sources
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cc/load-spz.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cc/splat-c-types.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cc/splat-types.cc"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cc/zstd-batch.cc"
)

set(spz_headers
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cc/load-spz.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cc/splat-c-types.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/src/cc/splat-types.h"
)

if(SPZ_BUILD_EXTENSIONS)
  list(APPEND spz_sources
    "${CMAKE_CURRENT_SOURCE_DIR}/extensions/cc/splat-extensions.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/extensions/cc/safe-orbit-camera-adobe.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/extensions/cc/coordinate-system-adobe.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/extensions/cc/georeference-niantic.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/extensions/cc/custom-attribute-niantic.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/extensions/cc/extension-frames.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/extensions/cc/quantized-planes.cc"
  )
  list(APPEND spz_headers
    "${CMAKE_CURRENT_SOURCE_DIR}/extensions/cc/splat-extensions.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/extensions/cc/safe-orbit-camera-adobe.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/extensions/cc/coordinate-system-adobe.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/extensions/cc/georeference-niantic.h"
    "${CMAKE_CURRENT_SOURCE_DIR}/extensions/cc/custom-attribute-niantic.h"
  )
endif()

# create the library and configure it
add_library(spz ${spz_sources})
add_library(spz::spz ALIAS spz)

if(NOT SPZ_USE_EMSCRIPTEN_ZLIB)
  target_link_libraries(spz PRIVATE ${zlib_linkage_target} zstd::zstd)
else()
  target_compile_options(spz PUBLIC --use-port=zlib)
  target_link_options(spz PUBLIC --use-port=zlib)
  target_link_libraries(spz PRIVATE zstd::zstd)
endif()

if(ANDROID)
  target_link_libraries(spz PRIVATE log)
endif()

target_include_directories(spz
  PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/src/cc>
  INTERFACE $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

if(SPZ_BUILD_EXTENSIONS)
  target_include_directories(spz
    PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/extensions/cc>
  )
  target_compile_definitions(spz PUBLIC SPZ_BUILD_EXTENSIONS)
endif()

if (CMAKE_CXX_STANDARD)
    set(SPZ_CXX ${CMAKE_CXX_STANDARD})
else()
    set(SPZ_CXX 17)
endif()

set_target_properties(spz PROPERTIES
  PUBLIC_HEADER "${spz_headers}"
  CXX_STANDARD ${SPZ_CXX}
  CXX_STANDARD_REQUIRED ON
)

if(SPZ_BUILD_WASM)
  file(MAKE_DIRECTORY "${SPZ_WASM_OUTPUT_DIR}")

  set(spz_binding_sources
    "${CMAKE_CURRENT_SOURCE_DIR}/src/emscripten/spz-bindings.cc"
    "${CMAKE_CURRENT_SOURCE_DIR}/src/emscripten/chunked-unpack.cc"
  )

  if(SPZ_BUILD_EXTENSIONS)
    list(APPEND spz_binding_sources
      "${CMAKE_CURRENT_SOURCE_DIR}/extensions/emscripten/splat-extensions.cc"
    )
  endif()

  add_executable(spz-wasm ${spz_binding_sources})
  target_link_libraries(spz-wasm PRIVATE spz)

  if(SPZ_BUILD_EXTENSIONS)
    target_include_directories(spz-wasm PRIVATE
      $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    )
  endif()

  target_link_options(spz-wasm PRIVATE
    "$<$<NOT:$<CONFIG:Debug>>:-O3>"
    "--bind"
    "-sMODULARIZE=1"
    "-sEXPORT_NAME=createSpzModule"
    "-sALLOW_MEMORY_GROWTH=1"
    "-sMAXIMUM_MEMORY=4294967296"
    "-sENVIRONMENT=web"
    "-sEXPORT_ES6=1"
    "-sSINGLE_FILE=1"
    "-sEXPORTED_RUNTIME_METHODS=[\"ccall\",\"cwrap\",\"HEAPU8\"]"
    "-sEXPORTED_FUNCTIONS=[\"_malloc\",\"_free\"]"
  )

  set_target_properties(spz-wasm PROPERTIES
    OUTPUT_NAME "spz"
    RUNTIME_OUTPUT_DIRECTORY "${SPZ_WASM_OUTPUT_DIR}"
  )

  # Generate the TypeScript declaration file conditionally based on SPZ_BUILD_EXTENSIONS
  if(SPZ_BUILD_EXTENSIONS)
    # Read the extensions TypeScript definitions
    file(READ "${CMAKE_CURRENT_SOURCE_DIR}/extensions/emscripten/splat-extensions.d.ts.in" SPZ_EXTENSIONS_TYPESCRIPT_DEFINITIONS)
    file(READ "${CMAKE_CURRENT_SOURCE_DIR}/extensions/emscripten/splat-extensions-modules.d.ts.in" SPZ_EXTENSIONS_MODULE_INTERFACE)
    set(SPZ_EXTENSIONS_GAUSSIAN_CLOUD_FIELD "  extensions: SpzExtensionBase[];")
  else()
    set(SPZ_EXTENSIONS_TYPESCRIPT_DEFINITIONS "")
    set(SPZ_EXTENSIONS_GAUSSIAN_CLOUD_FIELD "")
    set(SPZ_EXTENSIONS_MODULE_INTERFACE "")
  endif()

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/src/emscripten/spz.d.ts.in"
    "${SPZ_WASM_OUTPUT_DIR}/spz.d.ts"
    @ONLY
  )

  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/src/emscripten/spz-helpers.js"
    "${SPZ_WASM_OUTPUT_DIR}/spz-helpers.js"
    COPYONLY
  )
  configure_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/src/emscripten/spz-helpers.d.ts"
    "${SPZ_WASM_OUTPUT_DIR}/spz-helpers.d.ts"
    COPYONLY
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

if (SPZ_SHOULD_INSTALL)
  install(TARGETS spz
    EXPORT spzTargets
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    PUBLIC_HEADER DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}/spz"
  )

  install(EXPORT spzTargets
    NAMESPACE spz::
    DESTINATION "${CMAKE_INSTALL_LIBDIR}/cmake/spz"
  )
endif()

if(SPZ_BUILD_TOOLS)
  # CLI Tools
  add_executable(ply_to_spz cli_tools/src/ply_to_spz.cpp)
  target_link_libraries(ply_to_spz PRIVATE spz)
  target_include_directories(ply_to_spz PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
  set_target_properties(ply_to_spz PROPERTIES
    CXX_STANDARD ${SPZ_CXX}
    CXX_STANDARD_REQUIRED ON
  )

  add_executable(spz_to_ply cli_tools/src/spz_to_ply.cpp)
  target_link_libraries(spz_to_ply PRIVATE spz)
  target_include_directories(spz_to_ply PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
  set_target_properties(spz_to_ply PROPERTIES
    CXX_STANDARD ${SPZ_CXX}
    CXX_STANDARD_REQUIRED ON
  )

  add_executable(spz_info cli_tools/src/spz_info.cpp)
  target_link_libraries(spz_info PRIVATE spz)
  target_include_directories(spz_info PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src)
  set_target_properties(spz_info PROPERTIES
    CXX_STANDARD ${SPZ_CXX}
    CXX_STANDARD_REQUIRED ON
  )

  install(TARGETS ply_to_spz spz_to_ply spz_info
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
  )

  # Round-trip smoke tests through the CLI tools. Only when spz is the top-level
  # project so that consumers using add_subdirectory don't inherit our tests.
  if(BUILD_TESTING AND CMAKE_SOURCE_DIR STREQUAL PROJECT_SOURCE_DIR)
    enable_testing()
    set(spz_sample "${CMAKE_CURRENT_SOURCE_DIR}/samples/hornedlizard.spz")
    set(spz_sample_point_count 786233)
    set(spz_test_dir "${CMAKE_CURRENT_BINARY_DIR}/test-output")
    file(MAKE_DIRECTORY "${spz_test_dir}")

    add_test(NAME spz_info.sample
      COMMAND spz_info "${spz_sample}")
    add_test(NAME spz_to_ply.sample
      COMMAND spz_to_ply "${spz_sample}" "${spz_test_dir}/roundtrip.ply")
    add_test(NAME ply_to_spz.roundtrip
      COMMAND ply_to_spz "${spz_test_dir}/roundtrip.ply" "${spz_test_dir}/roundtrip.spz")
    add_test(NAME spz_info.roundtrip
      COMMAND spz_info "${spz_test_dir}/roundtrip.spz")

    # The CLI tools exit 0 when they cannot open a file, so assert on the point
    # count rather than trusting the exit code.
    set_tests_properties(spz_info.sample spz_info.roundtrip PROPERTIES
      PASS_REGULAR_EXPRESSION "Number of points: ${spz_sample_point_count}\n")
    set_tests_properties(spz_to_ply.sample PROPERTIES FIXTURES_SETUP spz_ply)
    set_tests_properties(ply_to_spz.roundtrip PROPERTIES
      FIXTURES_REQUIRED spz_ply FIXTURES_SETUP spz_roundtrip)
    set_tests_properties(spz_info.roundtrip PROPERTIES FIXTURES_REQUIRED spz_roundtrip)
  endif()
endif()

# --- Python Bindings (nanobind) ---
if(SPZ_BUILD_PYTHON_BINDINGS)
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
  # The ARCHIVE_OUTPUT_NAME is also set to avoid name clashes on Windows.
  set_target_properties(spz_py PROPERTIES
      OUTPUT_NAME spz
      ARCHIVE_OUTPUT_NAME spz_py
  )

  # Add the project root to the include path for the python module
  target_include_directories(spz_py PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

  # Add extensions include directory and source file if extensions are enabled
  if(SPZ_BUILD_EXTENSIONS)
    target_include_directories(spz_py PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/extensions/cc)
    target_sources(spz_py PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/extensions/python/splat-extensions.cc)
  endif()

  # Link the python module against the spz library
  target_link_libraries(spz_py PRIVATE spz::spz)

  # Install directive for scikit-build-core
  install(TARGETS spz_py LIBRARY DESTINATION spz)
endif()
