#!/usr/bin/env bash
# Build one CMake preset (incrementally) and run its CTest suite.
#
#   scripts/test.sh                  # default preset (dev)
#   scripts/test.sh dev-ext
#   scripts/test.sh dev -R roundtrip # extra arguments go to ctest
#
# Only presets with a matching test preset can be tested; see --list.
# For the Python bindings, use scripts/test-python.sh.
set -euo pipefail

readonly DEFAULT_PRESET="dev"

repo_root="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$repo_root"

die() {
  echo "error: $*" >&2
  exit 1
}

usage() {
  sed -n '2,9s/^# \{0,1\}//p' "${BASH_SOURCE[0]}"
}

has_test_preset() {
  cmake --list-presets=test 2>/dev/null | grep -q "^  \"$1\""
}

main() {
  case "${1:-}" in
    -h|--help) usage; exit 0 ;;
    --list) cmake --list-presets=test; exit 0 ;;
  esac

  local preset="${1:-$DEFAULT_PRESET}"
  if (( $# > 0 )); then shift; fi

  command -v cmake >/dev/null 2>&1 || die "'cmake' not found on PATH."
  has_test_preset "$preset" \
    || die "preset '${preset}' has no tests. Testable presets: $(cmake --list-presets=test | grep '^  "' | tr -d ' "' | paste -sd, -)"

  # build.sh owns the remaining prerequisite checks and configuration.
  "${repo_root}/scripts/build.sh" "$preset"

  ctest --preset "$preset" "$@"
}

main "$@"
