#include <algorithm>
#include <iostream>
#include <string>
#include "cc/load-spz.h"
#ifdef SPZ_BUILD_EXTENSIONS
#include "temporal-niantic.h"

namespace {
std::string encodingName(const spz::SpzTemporalChannel& channel) {
    switch (channel.encoding) {
        case spz::SpzChannelEncoding::Float32: return "float32";
        case spz::SpzChannelEncoding::Range: return std::to_string(channel.bits) + "-bit range";
        case spz::SpzChannelEncoding::Step: return "step";
    }
    return "unknown encoding";
}

const char* loopName(spz::SpzTemporalLoop loop) {
    switch (loop) {
        case spz::SpzTemporalLoop::Once: return "once";
        case spz::SpzTemporalLoop::Loop: return "loop";
        case spz::SpzTemporalLoop::PingPong: return "ping-pong";
    }
    return "unknown";
}

void printClipInfo(const spz::GaussianCloud& cloud) {
    const std::optional<spz::SpzTemporalModel> model = spz::bindTemporalModel(cloud);
    if (!model) return;
    const spz::SpzExtensionTemporalNiantic& clip = *model->clip;
    std::cout << "Animated clip: " << clip.clipStart << "s to " << clip.clipEnd
              << "s, poster frame at " << clip.posterTime << "s, " << loopName(clip.loop);
    if (clip.nominalFps > 0.0f) std::cout << ", " << clip.nominalFps << " fps";
    std::cout << std::endl;
    for (const spz::SpzTemporalChannel& channel : clip.channels) {
        std::cout << "  channel target " << static_cast<int>(channel.target) << " type "
                  << static_cast<int>(channel.type) << ": dim " << channel.dim << ", "
                  << encodingName(channel);
        const std::vector<float> maxError = spz::channelMaxError(channel);
        if (!maxError.empty()) {
            // A skin channel's first half are driver indices, which round back exactly.
            const bool skin = channel.type == spz::SpzTemporalChannelType::Skin;
            const auto first = maxError.begin() + (skin ? maxError.size() / 2 : 0);
            std::cout << ", max error " << *std::max_element(first, maxError.end())
                      << (skin ? " (weights; indices exact)" : "");
        }
        if (channel.pointRange) {
            std::cout << ", points " << channel.pointRange->begin << "-" << channel.pointRange->end;
        }
        std::cout << std::endl;
    }
    if (clip.driverCount > 0) {
        std::cout << "  " << clip.driverCount << " drivers, " << clip.driverKeyframeCount()
                  << " keyframes";
        if (!clip.driverNames.empty()) {
            std::cout << ":";
            for (const std::string& name : clip.driverNames) std::cout << " " << name;
        }
        std::cout << std::endl;
    }
}
}  // namespace
#endif

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
#ifdef SPZ_BUILD_EXTENSIONS
        unpack_options.temporal = spz::TemporalLoad::Include;
#endif
        spz::GaussianCloud cloud = spz::loadSpz(argv[1], unpack_options);
        printCloudInfo(cloud);
#ifdef SPZ_BUILD_EXTENSIONS
        printClipInfo(cloud);
#endif
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
} 