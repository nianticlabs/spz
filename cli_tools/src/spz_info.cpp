#include <iostream>
#include "cc/load-spz.h"

void printCloudInfo(const spz::GaussianCloud& cloud) {
    std::cout << "Number of points: " << cloud.positions.size() / 3 << std::endl;
    
    if (!cloud.positions.empty()) {
        float min_x = cloud.positions[0];
        float max_x = cloud.positions[0];
        float min_y = cloud.positions[1];
        float max_y = cloud.positions[1];
        float min_z = cloud.positions[2];
        float max_z = cloud.positions[2];

        for (size_t i = 0; i < cloud.positions.size(); i += 3) {
            min_x = std::min(min_x, cloud.positions[i]);
            max_x = std::max(max_x, cloud.positions[i]);
            min_y = std::min(min_y, cloud.positions[i + 1]);
            max_y = std::max(max_y, cloud.positions[i + 1]);
            min_z = std::min(min_z, cloud.positions[i + 2]);
            max_z = std::max(max_z, cloud.positions[i + 2]);
        }

        std::cout << "Bounding box:" << std::endl;
        std::cout << "  X: " << min_x << " to " << max_x << std::endl;
        std::cout << "  Y: " << min_y << " to " << max_y << std::endl;
        std::cout << "  Z: " << min_z << " to " << max_z << std::endl;
    }
}

int main(int argc, char *argv[]) {
    if (argc < 2) {
        std::cerr << "Usage: spz_info <input.spz>" << std::endl;
        return 1;
    }

    try {
        // Load SPZ File
        spz::UnpackOptions unpack_options;
        spz::GaussianCloud cloud = spz::loadSpz(argv[1], unpack_options);
        printCloudInfo(cloud);
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
} 