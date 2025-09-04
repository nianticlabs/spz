#include <iostream>
#include "cc/load-spz.h"

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: spz_to_ply <input.spz> <output.ply>" << std::endl;
        return 1;
    }

    try {
        // Load SPZ File
        spz::UnpackOptions unpack_options;
        spz::GaussianCloud splat = spz::loadSpz(argv[1], unpack_options);

        // Save PLY File
        spz::PackOptions pack_options;
        spz::saveSplatToPly(splat, pack_options, argv[2]);

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
} 