#include <iostream>
#include <string>
#include "cc/load-spz.h"
#ifdef SPZ_BUILD_EXTENSIONS
#include "temporal-niantic.h"
#endif

namespace {
constexpr char kTimeFlag[] = "--time";
constexpr int kArgsWithTime = 5;
}  // namespace

int main(int argc, char *argv[]) {
    const bool hasTime = argc == kArgsWithTime && std::string(argv[3]) == kTimeFlag;
    if (argc < 3 || (argc > 3 && !hasTime)) {
        std::cerr << "Usage: spz_to_ply <input.spz> <output.ply> [--time <seconds>]\n"
                     "  --time exports an animated clip's frame at that clip time; without it the\n"
                     "  poster frame is exported." << std::endl;
        return 1;
    }

    try {
        // Load SPZ File
        spz::UnpackOptions unpack_options;
#ifdef SPZ_BUILD_EXTENSIONS
        if (hasTime) unpack_options.temporal = spz::TemporalLoad::Include;
#endif
        spz::GaussianCloud splat = spz::loadSpz(argv[1], unpack_options);

        if (hasTime) {
#ifdef SPZ_BUILD_EXTENSIONS
            const double t = std::stod(argv[4]);
            const std::optional<spz::SpzTemporalModel> model = spz::bindTemporalModel(splat);
            if (!model) {
                std::cerr << "Error: " << argv[1] << " isn't a valid animated clip" << std::endl;
                return 1;
            }
            splat = spz::evaluateTemporal(splat, *model, t);
#else
            std::cerr << "Error: --time needs a build with SPZ_BUILD_EXTENSIONS" << std::endl;
            return 1;
#endif
        }

        // Save PLY File
        spz::PackOptions pack_options;
        spz::saveSplatToPly(splat, pack_options, argv[2]);

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
} 