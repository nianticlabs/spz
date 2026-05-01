#include <cctype>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "cc/load-spz.h"

static std::string toUpper(std::string s) {
  for (char &c : s) {
    c = static_cast<char>(std::toupper(static_cast<unsigned char>(c)));
  }
  return s;
}

static spz::CoordinateSystem parseCoordinateSystem(const char *name) {
  const std::string u = toUpper(name);
  if (u == "UNSPECIFIED")
    return spz::CoordinateSystem::UNSPECIFIED;
  if (u == "LDB")
    return spz::CoordinateSystem::LDB;
  if (u == "RDB")
    return spz::CoordinateSystem::RDB;
  if (u == "LUB")
    return spz::CoordinateSystem::LUB;
  if (u == "RUB")
    return spz::CoordinateSystem::RUB;
  if (u == "LDF")
    return spz::CoordinateSystem::LDF;
  if (u == "RDF")
    return spz::CoordinateSystem::RDF;
  if (u == "LUF")
    return spz::CoordinateSystem::LUF;
  if (u == "RUF")
    return spz::CoordinateSystem::RUF;
  if (u == "LBD")
    return spz::CoordinateSystem::LBD;
  if (u == "RBD")
    return spz::CoordinateSystem::RBD;
  if (u == "LBU")
    return spz::CoordinateSystem::LBU;
  if (u == "RBU")
    return spz::CoordinateSystem::RBU;
  if (u == "LFD")
    return spz::CoordinateSystem::LFD;
  if (u == "RFD")
    return spz::CoordinateSystem::RFD;
  if (u == "LFU")
    return spz::CoordinateSystem::LFU;
  if (u == "RFU")
    return spz::CoordinateSystem::RFU;
  throw std::invalid_argument(std::string("unknown coordinate system: ") + name);
}

static void printUsage(const char *prog) {
  std::cerr << "Usage: " << prog << " [--from <COORD>] <input.ply> <output.spz>\n"
            << "  --from COORD   source frame of the PLY data for packing (default: UNSPECIFIED).\n"
            << "                 Examples: RDF (typical 3DGS PLY), RUB, GLB, LUF, RUF, ...\n";
}

int main(int argc, char *argv[]) {
  const char *fromArg = nullptr;
  std::vector<const char *> positional;
  for (int i = 1; i < argc; ++i) {
    if (std::strcmp(argv[i], "--from") == 0) {
      if (i + 1 >= argc) {
        printUsage(argv[0]);
        std::cerr << "Error: --from requires a value\n";
        return 1;
      }
      fromArg = argv[++i];
    } else if (argv[i][0] == '-') {
      printUsage(argv[0]);
      std::cerr << "Error: unknown option " << argv[i] << "\n";
      return 1;
    } else {
      positional.push_back(argv[i]);
    }
  }

  if (positional.size() != 2) {
    printUsage(argv[0]);
    return 1;
  }

  try {
    spz::UnpackOptions unpack_options;
    spz::GaussianCloud splat = spz::loadSplatFromPly(positional[0], unpack_options);

    spz::PackOptions pack_options;
    if (fromArg != nullptr) {
      pack_options.from = parseCoordinateSystem(fromArg);
    }

    spz::saveSpz(splat, pack_options, positional[1]);

    return 0;
  } catch (const std::exception &e) {
    std::cerr << "Error: " << e.what() << std::endl;
    return 1;
  }
}
