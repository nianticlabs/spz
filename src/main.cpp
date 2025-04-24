#include "load-spz.h"
#include <iostream>
#include <filesystem>
#include <string>

namespace fs = std::filesystem;

int main(int argc, char* argv[]) {
    bool antialiased = false;
    bool fast_rot_quantization = false;
    int arg_offset = 1;

    // Parse flags
    while (arg_offset < argc && argv[arg_offset][0] == '-') {
        std::string flag = argv[arg_offset];
        if (flag == "-a") {
            antialiased = true;
        }
        else if (flag == "-f") {
            fast_rot_quantization = true;
        }
        arg_offset++;
    }

    // Check remaining arguments
    if (argc - arg_offset != 2) {
        std::cerr << "Usage: spz_convert [-a] [-f] <input_file> <output_file>\n";
        std::cerr << "Options:\n";
        std::cerr << "  -a    Enable antialiasing\n";
        std::cerr << "  -f    Enable fast rotation quantization\n";
        std::cerr << "Supported conversions:\n";
        std::cerr << "  .ply -> .spz\n";
        std::cerr << "  .spz -> .ply\n";
        return 1;
    }

    std::string input_file = argv[arg_offset];
    std::string output_file = argv[arg_offset + 1];

    // Determine the input file extension
    fs::path input_path(input_file);
    if (!fs::exists(input_path)) {
        std::cerr << "Input file does not exist: " << input_file << "\n";
        return 1;
    }

    std::string input_extension = input_path.extension().string();
    // Convert extension to lowercase for case-insensitive comparison
    std::transform(input_extension.begin(), input_extension.end(), input_extension.begin(),
                   [](unsigned char c) { return std::tolower(c); });

    if (input_extension == ".ply") {
        // Convert PLY to SPZ
        spz::GaussianCloud cloud = spz::loadSplatFromPly(input_file, spz::UnpackOptions());
        if (cloud.numPoints == 0) {
            std::cerr << "Failed to load GaussianCloud from " << input_file << "\n";
            return 1;
        }
        cloud.antialiased = antialiased;

        spz::PackOptions po;
        po.fast_rot_quantization = fast_rot_quantization;
        bool success = spz::saveSpz(cloud, po, output_file);
        if (!success) {
            std::cerr << "Failed to save GaussianCloud to " << output_file << "\n";
            return 1;
        }

        std::cout << "Successfully converted " << input_file << " to " << output_file << "\n";
    }
    else if (input_extension == ".spz") {
        // Convert SPZ to PLY
        spz::GaussianCloud cloud = spz::loadSpz(input_file, spz::UnpackOptions());
        if (cloud.numPoints == 0) {
            std::cerr << "Failed to load GaussianCloud from " << input_file << "\n";
            return 1;
        }

        bool success = spz::saveSplatToPly(cloud, spz::PackOptions(), output_file);
        if (!success) {
            std::cerr << "Failed to save GaussianCloud to " << output_file << "\n";
            return 1;
        }

        std::cout << "Successfully converted " << input_file << " to " << output_file << "\n";
    }
    else {
        std::cerr << "Unsupported input file extension: " << input_extension << "\n";
        std::cerr << "Supported extensions are .ply and .spz\n";
        return 1;
    }

    return 0;
} 