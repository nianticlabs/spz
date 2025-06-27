#include <iostream>
#include "cc/load-spz.h"

int main(int argc, char *argv[]) {
    if (argc < 3) {
        std::cerr << "Usage: ply_to_spz <input.ply> <output.spz>" << std::endl;
        return 1;
    }

    try {
        // Load PLY File
        spz::UnpackOptions unpack_options;
        spz::GaussianCloud splat = spz::loadSplatFromPly(argv[1], unpack_options);

        // Save SPZ File
        spz::PackOptions pack_options;
        spz::saveSpz(splat, pack_options, argv[2]);

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
} 