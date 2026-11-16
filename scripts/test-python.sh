#!/usr/bin/env bash
# Build the Python bindings into a local virtualenv and run the pytest suite.
#
#   scripts/test-python.sh                 # bindings without extensions
#   scripts/test-python.sh --ext           # bindings with SPZ_BUILD_EXTENSIONS=ON
#   scripts/test-python.sh --ext -- -k io  # arguments after -- go to pytest
#
# Uses $PYTHON (default: python3) to create .venv/ on first run.
set -euo pipefail

readonly VENV_DIR=".venv"
readonly TEST_DIR="tests/python"
readonly BUILD_DIR_PLAIN="build/py/{wheel_tag}"
readonly BUILD_DIR_EXT="build/py-ext/{wheel_tag}"

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

die() {
  echo "error: $*" >&2
  exit 1
}

usage() {
  sed -n '2,8s/^# \{0,1\}//p' "${BASH_SOURCE[0]}"
}

main() {
  local extensions="OFF"
  local build_dir="$BUILD_DIR_PLAIN"
  local pytest_args=()

  while (( $# > 0 )); do
    case "$1" in
      --ext) extensions="ON"; build_dir="$BUILD_DIR_EXT" ;;
      -h|--help) usage; exit 0 ;;
      --) shift; pytest_args=("$@"); break ;;
      *) die "unknown argument '$1' (see --help)" ;;
    esac
    shift
  done

  local python="${PYTHON:-python3}"
  command -v "$python" >/dev/null 2>&1 || die "'$python' not found; set PYTHON to a Python 3.8+ interpreter."

  if [[ ! -x "${VENV_DIR}/bin/python" ]]; then
    "$python" -m venv "$VENV_DIR"
  fi
  local venv_python="${VENV_DIR}/bin/python"

  # Separate build dirs keep the extension and non-extension builds from
  # reconfiguring each other on every switch.
  SKBUILD_BUILD_DIR="$build_dir" \
  SKBUILD_CMAKE_DEFINE="SPZ_BUILD_EXTENSIONS=${extensions}" \
    "$venv_python" -m pip install --quiet ".[tests]"

  "$venv_python" -m pytest "$TEST_DIR" "${pytest_args[@]+"${pytest_args[@]}"}"
}

main "$@"
