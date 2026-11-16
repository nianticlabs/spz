#!/usr/bin/env bash
# Configure and build one CMake preset. To run tests, use scripts/test.sh.
#
#   scripts/build.sh                 # default preset (dev)
#   scripts/build.sh wasm-ext        # WASM with extensions, via emcmake
#   scripts/build.sh dev -DBUILD_SHARED_LIBS=ON --fresh
#   scripts/build.sh --list          # show available presets
#
# Extra arguments after the preset name are passed to the configure step.
# Output goes to build/<preset>/ (WASM artifacts: build/<preset>/dist/).
set -euo pipefail

readonly DEFAULT_PRESET="dev"
readonly MIN_CMAKE_MAJOR=3
readonly MIN_CMAKE_MINOR=21
readonly WASM_PRESET_PREFIX="wasm"
readonly ANDROID_PRESET_PREFIX="android"
readonly BUILD_ROOT="build"

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

die() {
  echo "error: $*" >&2
  exit 1
}

usage() {
  sed -n '2,10s/^# \{0,1\}//p' "${BASH_SOURCE[0]}"
}

require_command() {
  command -v "$1" >/dev/null 2>&1 || die "'$1' not found on PATH. $2"
}

check_cmake_version() {
  require_command cmake "Install CMake ${MIN_CMAKE_MAJOR}.${MIN_CMAKE_MINOR}+ (e.g. 'brew install cmake' or 'pip install cmake')."
  local version major minor
  version="$(cmake --version | head -n1 | awk '{print $3}')"
  major="${version%%.*}"
  minor="$(echo "$version" | cut -d. -f2)"
  if (( major < MIN_CMAKE_MAJOR || (major == MIN_CMAKE_MAJOR && minor < MIN_CMAKE_MINOR) )); then
    die "CMake ${version} is too old for presets; need ${MIN_CMAKE_MAJOR}.${MIN_CMAKE_MINOR}+."
  fi
}

# Prefer Ninja for fresh build trees; an existing tree keeps its generator because
# CMake refuses to switch generators in place.
choose_generator() {
  local preset="$1"
  if [[ -z "${CMAKE_GENERATOR:-}" && ! -f "${BUILD_ROOT}/${preset}/CMakeCache.txt" ]] \
      && command -v ninja >/dev/null 2>&1; then
    export CMAKE_GENERATOR="Ninja"
  fi
}

main() {
  if [[ "${1:-}" == "-h" || "${1:-}" == "--help" ]]; then
    usage
    exit 0
  fi

  check_cmake_version

  if [[ "${1:-}" == "--list" ]]; then
    cmake --list-presets=all
    exit 0
  fi

  local preset="${1:-$DEFAULT_PRESET}"
  if (( $# > 0 )); then shift; fi

  choose_generator "$preset"

  case "$preset" in
    "${WASM_PRESET_PREFIX}"*)
      require_command emcmake "Install and activate the Emscripten SDK: https://emscripten.org/docs/getting_started/downloads.html"
      emcmake cmake --preset "$preset" "$@"
      ;;
    "${ANDROID_PRESET_PREFIX}"*)
      [[ -n "${ANDROID_NDK_HOME:-}" ]] || die "ANDROID_NDK_HOME is not set; point it at your NDK (r27+)."
      [[ -d "$ANDROID_NDK_HOME" ]] || die "ANDROID_NDK_HOME='${ANDROID_NDK_HOME}' is not a directory."
      cmake --preset "$preset" "$@"
      ;;
    *)
      cmake --preset "$preset" "$@"
      ;;
  esac

  cmake --build --preset "$preset"

  echo "Built preset '${preset}' into ${BUILD_ROOT}/${preset}/"
}

main "$@"
