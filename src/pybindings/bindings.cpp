/*
MIT License

Copyright (c) 2025 Adobe Inc.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>  // For STL containers like std::vector
#include <pybind11/numpy.h>

#include <sstream>  // For std::ostringstream

#include "load-spz.h"
#include "splat-c-types.h"
#include "splat-extensions.h"
#include "splat-types.h"

namespace py = pybind11;
using namespace spz;

py::bytes saveSpzToBytes(const GaussianCloud& g, const PackOptions& options) {
  std::vector<uint8_t> buffer;
  if (!saveSpz(g, options, &buffer)) {
    throw std::runtime_error("saveSpz failed");
  }
  return py::bytes(reinterpret_cast<const char*>(buffer.data()), buffer.size());
}

// Overload with default options
py::bytes saveSpzToBytesDefault(const GaussianCloud& g) {
  PackOptions options;
  return saveSpzToBytes(g, options);
}

// Wrapper functions with default options
bool saveSpzDefault(const GaussianCloud& g, const std::string& filename) {
  PackOptions options;
  return saveSpz(g, options, filename);
}

bool saveSpzVectorDefault(const GaussianCloud& g, std::vector<uint8_t>* output) {
  PackOptions options;
  return saveSpz(g, options, output);
}

GaussianCloud loadSpzDefault(const std::string& filename) {
  UnpackOptions options;
  return loadSpz(filename, options);
}

GaussianCloud loadSpzVectorDefault(const std::vector<uint8_t>& data) {
  UnpackOptions options;
  return loadSpz(data, options);
}

void saveSplatToPlyDefault(const GaussianCloud& gaussians, const std::string& filename) {
  PackOptions options;
  saveSplatToPly(gaussians, options, filename);
}

GaussianCloud loadSplatFromPlyDefault(const std::string& filename) {
  UnpackOptions options;
  return loadSplatFromPly(filename, options);
}

py::bytes serializePackedGaussiansToBytes(const PackedGaussians& packed) {
  std::ostringstream oss;
  serializePackedGaussians(packed, &oss);
  std::string str = oss.str();
  return py::bytes(str);
}

py::bytes compressGzippedFromBytes(const py::bytes& data) {
  std::string str = data;
  std::vector<uint8_t> output;
  bool success = compressGzipped(reinterpret_cast<const uint8_t*>(str.data()), str.size(), &output);
  if (!success) {
    throw std::runtime_error("compressGzipped failed");
  }
  return py::bytes(reinterpret_cast<const char*>(output.data()), output.size());
}

PYBIND11_MODULE(spz_bindings, m) {
  m.doc() = "Python bindings for spz library";

  // Bind CoordinateSystem enum
  py::enum_<CoordinateSystem>(m, "CoordinateSystem")
      .value("UNSPECIFIED", CoordinateSystem::UNSPECIFIED)
      .value("LDB", CoordinateSystem::LDB)  // Left Down Back
      .value("RDB", CoordinateSystem::RDB)  // Right Down Back
      .value("LUB", CoordinateSystem::LUB)  // Left Up Back
      .value("RUB", CoordinateSystem::RUB)  // Right Up Back, Three.js coordinate system
      .value("LDF", CoordinateSystem::LDF)  // Left Down Front
      .value("RDF", CoordinateSystem::RDF)  // Right Down Front, PLY coordinate system
      .value("LUF", CoordinateSystem::LUF)  // Left Up Front, GLB coordinate system
      .value("RUF", CoordinateSystem::RUF)  // Right Up Front, Unity coordinate system
      .export_values();

  py::enum_<SpzExtensionType>(m, "SpzExtensionType", py::arithmetic())
      .value("SPZ_ADOBE_sh_quantization", SpzExtensionType::SPZ_ADOBE_sh_quantization)
      .value("SPZ_ADOBE_safe_orbit_camera", SpzExtensionType::SPZ_ADOBE_safe_orbit_camera)
      .export_values();

  py::class_<SpzExtensionBase, std::shared_ptr<SpzExtensionBase>>(m, "SpzExtensionBase");

  py::class_<SpzExtensionSHQuantizationAdobe, SpzExtensionBase, std::shared_ptr<SpzExtensionSHQuantizationAdobe>>(m, "SpzExtensionSHQuantizationAdobe")
      .def(py::init<>())
      .def_readwrite("sh1Bits", &SpzExtensionSHQuantizationAdobe::sh1Bits)
      .def_readwrite("shRestBits", &SpzExtensionSHQuantizationAdobe::shRestBits)
      .def_readwrite("shMin", &SpzExtensionSHQuantizationAdobe::shMin)
      .def_readwrite("shMax", &SpzExtensionSHQuantizationAdobe::shMax)
      .def_static("type", &SpzExtensionSHQuantizationAdobe::type);

  py::class_<SpzExtensionSafeOrbitCameraAdobe, SpzExtensionBase, std::shared_ptr<SpzExtensionSafeOrbitCameraAdobe>>(m, "SpzExtensionSafeOrbitCameraAdobe")
      .def(py::init<>())
      .def_readwrite("safeOrbitElevationMin", &SpzExtensionSafeOrbitCameraAdobe::safeOrbitElevationMin)
      .def_readwrite("safeOrbitElevationMax", &SpzExtensionSafeOrbitCameraAdobe::safeOrbitElevationMax)
      .def_readwrite("safeOrbitRadiusMin", &SpzExtensionSafeOrbitCameraAdobe::safeOrbitRadiusMin)
      .def_static("type", &SpzExtensionSafeOrbitCameraAdobe::type);

  // Bind CoordinateConverter
  py::class_<CoordinateConverter>(m, "CoordinateConverter")
      .def(py::init<>())
      .def_readwrite("flipP", &CoordinateConverter::flipP)
      .def_readwrite("flipQ", &CoordinateConverter::flipQ)
      .def_readwrite("flipSh", &CoordinateConverter::flipSh);

  // Bind coordinateConverter function
  m.def("coordinateConverter", &coordinateConverter, py::arg("from"), py::arg("to"));

  // Bind PackOptions
  py::class_<PackOptions>(m, "PackOptions")
      .def(py::init<>())
      .def_readwrite("version", &PackOptions::version)
      .def_readwrite("from", &PackOptions::from)
      .def_readwrite("sh1Bits", &PackOptions::sh1Bits)
      .def_readwrite("shRestBits", &PackOptions::shRestBits)
      .def_readwrite("enableSHMinMaxScaling", &PackOptions::enableSHMinMaxScaling);

  // Bind UnpackOptions
  py::class_<UnpackOptions>(m, "UnpackOptions").def(py::init<>()).def_readwrite("to", &UnpackOptions::to);

  // Bind UnpackedGaussian
  py::class_<UnpackedGaussian>(m, "UnpackedGaussian")
      .def(py::init<>())
      .def_readwrite("position", &UnpackedGaussian::position)
      .def_readwrite("rotation", &UnpackedGaussian::rotation)
      .def_readwrite("scale", &UnpackedGaussian::scale)
      .def_readwrite("color", &UnpackedGaussian::color)
      .def_readwrite("alpha", &UnpackedGaussian::alpha)
      .def_readwrite("shR", &UnpackedGaussian::shR)
      .def_readwrite("shG", &UnpackedGaussian::shG)
      .def_readwrite("shB", &UnpackedGaussian::shB);

  // Bind PackedGaussian
  py::class_<PackedGaussian>(m, "PackedGaussian")
      .def(py::init<>())
      .def("unpack",
           &PackedGaussian::unpack,
           py::arg("usesFloat16"),
           py::arg("fractionalBits"),
           py::arg("c"),
           py::arg("shMin") = -1.0f,
           py::arg("shMax") = 1.0f);

  // Bind PackedGaussians
  py::class_<PackedGaussians>(m, "PackedGaussians")
      .def(py::init<>())
      .def_readwrite("version", &PackedGaussians::version)
      .def_readwrite("numPoints", &PackedGaussians::numPoints)
      .def_readwrite("shDegree", &PackedGaussians::shDegree)
      .def_readwrite("fractionalBits", &PackedGaussians::fractionalBits)
      .def_readwrite("antialiased", &PackedGaussians::antialiased)
      .def_readwrite("extensions", &PackedGaussians::extensions)
      .def_readwrite("positions", &PackedGaussians::positions)
      .def_readwrite("scales", &PackedGaussians::scales)
      .def_readwrite("rotations", &PackedGaussians::rotations)
      .def_readwrite("alphas", &PackedGaussians::alphas)
      .def_readwrite("colors", &PackedGaussians::colors)
      .def_readwrite("sh", &PackedGaussians::sh)
      .def("usesFloat16", &PackedGaussians::usesFloat16)
      .def("at", &PackedGaussians::at)
      .def("unpack", &PackedGaussians::unpack, py::arg("i"), py::arg("c"));

  // Bind SpzFloatBuffer
  py::class_<SpzFloatBuffer>(m, "SpzFloatBuffer").def(py::init<>()).def_readwrite("count", &SpzFloatBuffer::count).def_readwrite("data", &SpzFloatBuffer::data);

  py::class_<SpzExtensionNode>(m, "SpzExtensionNode")
      .def_property_readonly("type", [](SpzExtensionNode& self) { return self.type; })
      .def_property_readonly("next", [](SpzExtensionNode& self) -> SpzExtensionNode* { return self.next; })
      .def_property_readonly("dataAddr", [](SpzExtensionNode& self) { return reinterpret_cast<std::uintptr_t>(self.data); });

  // Bind GaussianCloudData
  py::class_<GaussianCloudData>(m, "GaussianCloudData")
      .def(py::init<>())
      .def_readwrite("numPoints", &GaussianCloudData::numPoints)
      .def_readwrite("shDegree", &GaussianCloudData::shDegree)
      .def_readwrite("antialiased", &GaussianCloudData::antialiased)
      .def_readwrite("positions", &GaussianCloudData::positions)
      .def_readwrite("scales", &GaussianCloudData::scales)
      .def_readwrite("rotations", &GaussianCloudData::rotations)
      .def_readwrite("alphas", &GaussianCloudData::alphas)
      .def_readwrite("colors", &GaussianCloudData::colors)
      .def_readwrite("sh", &GaussianCloudData::sh)
      .def_property_readonly("extensions", [](GaussianCloudData& self) -> SpzExtensionNode* { return self.extensions; });

  // Bind GaussianCloud
  py::class_<GaussianCloud>(m, "GaussianCloud")
      .def(py::init<>())
      .def_readwrite("numPoints", &GaussianCloud::numPoints)
      .def_readwrite("shDegree", &GaussianCloud::shDegree)
      .def_readwrite("antialiased", &GaussianCloud::antialiased)
      .def_readwrite("extensions", &GaussianCloud::extensions)
      .def_readwrite("positions", &GaussianCloud::positions)
      .def_readwrite("scales", &GaussianCloud::scales)
      .def_readwrite("rotations", &GaussianCloud::rotations)
      .def_readwrite("alphas", &GaussianCloud::alphas)
      .def_readwrite("colors", &GaussianCloud::colors)
      .def_readwrite("sh", &GaussianCloud::sh)
      .def("data", &GaussianCloud::data)
      .def("convertCoordinates", &GaussianCloud::convertCoordinates, py::arg("from"), py::arg("to"))
      .def("rotate180DegAboutX", &GaussianCloud::rotate180DegAboutX)
      .def("medianVolume", &GaussianCloud::medianVolume);

  // Bind saveSpz with PackOptions (existing)
  m.def("saveSpz",
        py::overload_cast<const GaussianCloud&, const PackOptions&, std::vector<uint8_t>*>(&saveSpz),
        py::arg("g"),
        py::arg("options"),
        py::arg("output"));
  m.def("saveSpz",
        py::overload_cast<const GaussianCloud&, const PackOptions&, const std::string&>(&saveSpz),
        py::arg("g"),
        py::arg("options"),
        py::arg("filename"));

  // Bind saveSpz with default options (new)
  m.def("saveSpz", &saveSpzDefault, py::arg("g"), py::arg("filename"));
  m.def("saveSpz", &saveSpzVectorDefault, py::arg("g"), py::arg("output"));

  // Bind loadSpz with UnpackOptions (existing)
  m.def("loadSpz", py::overload_cast<const std::vector<uint8_t>&, const UnpackOptions&>(&loadSpz), py::arg("data"), py::arg("options"));
  m.def("loadSpz", py::overload_cast<const std::string&, const UnpackOptions&>(&loadSpz), py::arg("filename"), py::arg("options"));

  // Bind loadSpz with default options (new)
  m.def("loadSpz", &loadSpzDefault, py::arg("filename"));
  m.def("loadSpz", &loadSpzVectorDefault, py::arg("data"));

  // Bind loadSpzPacked
  m.def("loadSpzPacked", py::overload_cast<const std::string&>(&loadSpzPacked), py::arg("filename"));
  m.def("loadSpzPacked", py::overload_cast<const uint8_t*, int32_t>(&loadSpzPacked), py::arg("data"), py::arg("size"));
  m.def("loadSpzPacked", py::overload_cast<const std::vector<uint8_t>&>(&loadSpzPacked), py::arg("data"));

  // Bind saveSplatToPly and loadSplatFromPly with options (existing)
  m.def("saveSplatToPly", &saveSplatToPly, py::arg("gaussians"), py::arg("options"), py::arg("filename"));
  m.def("loadSplatFromPly", &loadSplatFromPly, py::arg("filename"), py::arg("options"));

  // Bind saveSplatToPly and loadSplatFromPly with default options (new)
  m.def("saveSplatToPly", &saveSplatToPlyDefault, py::arg("gaussians"), py::arg("filename"));
  m.def("loadSplatFromPly", &loadSplatFromPlyDefault, py::arg("filename"));

  // Bind saveSpzToBytes with options (existing)
  m.def("saveSpzToBytes", &saveSpzToBytes, py::arg("gaussians"), py::arg("options"));

  // Bind saveSpzToBytes with default options (new)
  m.def("saveSpzToBytes", &saveSpzToBytesDefault, py::arg("gaussians"));

  // Bind serializePackedGaussians
  m.def("serializePackedGaussians", &serializePackedGaussiansToBytes, py::arg("packed"));

  // Bind compressGzipped
  m.def("compressGzipped", &compressGzippedFromBytes, py::arg("data"));

  // Fast factory to build GaussianCloud directly from NumPy arrays without Python boxing
  m.def(
      "createGaussianCloudFromNumpy",
      [](py::array_t<float, py::array::c_style | py::array::forcecast> xyz,
         py::array_t<float, py::array::c_style | py::array::forcecast> opacities,
         py::array_t<float, py::array::c_style | py::array::forcecast> scales,
         py::array_t<float, py::array::c_style | py::array::forcecast> rotations,
         py::array_t<float, py::array::c_style | py::array::forcecast> rgb,
         py::object f_rest,
         int sh_degree) {
        // Basic shape validation
        auto xyz_info = xyz.request();
        auto op_info = opacities.request();
        auto sc_info = scales.request();
        auto rot_info = rotations.request();
        auto rgb_info = rgb.request();

        if (xyz_info.ndim != 2 || xyz_info.shape[1] != 3) {
          throw std::invalid_argument("xyz must have shape [N,3]");
        }
        if (!((op_info.ndim == 1) || (op_info.ndim == 2 && op_info.shape[1] == 1))) {
          throw std::invalid_argument("opacities must have shape [N] or [N,1]");
        }
        if (sc_info.ndim != 2 || sc_info.shape[1] != 3) {
          throw std::invalid_argument("scale must have shape [N,3]");
        }
        if (rot_info.ndim != 2 || rot_info.shape[1] != 4) {
          throw std::invalid_argument("rotation must have shape [N,4]");
        }
        if (rgb_info.ndim != 2 || rgb_info.shape[1] != 3) {
          throw std::invalid_argument("rgb must have shape [N,3]");
        }

        const int32_t n = static_cast<int32_t>(xyz_info.shape[0]);
        if (op_info.shape[0] != n || sc_info.shape[0] != n || rot_info.shape[0] != n || rgb_info.shape[0] != n) {
          throw std::invalid_argument("all inputs must have the same length N");
        }

        GaussianCloud g;
        g.numPoints = n;
        g.shDegree = sh_degree;

        // Bulk copy contiguous memory
        g.positions.assign(xyz.data(), xyz.data() + xyz.size());
        g.scales.assign(scales.data(), scales.data() + scales.size());
        g.colors.assign(rgb.data(), rgb.data() + rgb.size());

        // Handle opacities possibly shaped [N,1]
        if (op_info.ndim == 2) {
          // Flatten [N,1] to [N]
          g.alphas.resize(static_cast<size_t>(n));
          const float* src = static_cast<const float*>(op_info.ptr);
          for (int32_t i = 0; i < n; ++i) {
            g.alphas[static_cast<size_t>(i)] = src[static_cast<size_t>(i) * static_cast<size_t>(op_info.shape[1])];
          }
        } else {
          g.alphas.assign(opacities.data(), opacities.data() + opacities.size());
        }

        // Normalize rotations per-row and copy
        g.rotations.resize(static_cast<size_t>(n) * 4);
        const float* r = rotations.data();
        for (int32_t i = 0; i < n; ++i) {
          const float* q = r + static_cast<size_t>(i) * 4;
          float s = std::sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);
          if (s == 0.0f) s = 1.0f;
          g.rotations[static_cast<size_t>(i) * 4 + 0] = q[0] / s;
          g.rotations[static_cast<size_t>(i) * 4 + 1] = q[1] / s;
          g.rotations[static_cast<size_t>(i) * 4 + 2] = q[2] / s;
          g.rotations[static_cast<size_t>(i) * 4 + 3] = q[3] / s;
        }

        if (!f_rest.is_none()) {
          auto f = py::cast<py::array_t<float, py::array::c_style | py::array::forcecast>>(f_rest);
          g.sh.assign(f.data(), f.data() + f.size());
        }

        return g;
      },
      py::arg("xyz"),
      py::arg("opacities"),
      py::arg("scale"),
      py::arg("rotation"),
      py::arg("rgb"),
      py::arg("f_rest") = py::none(),
      py::arg("sh_degree") = 0);

  m.attr("LATEST_SPZ_HEADER_VERSION") = spz::LATEST_SPZ_HEADER_VERSION;
}
