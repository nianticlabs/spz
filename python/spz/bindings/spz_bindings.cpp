#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <string>
#include <stdexcept>

#include "load-spz.h"
#include "splat-types.h"

namespace py = pybind11;

// Helper function to convert string to CoordinateSystem enum
spz::CoordinateSystem stringToCoordinateSystem(const std::string& cs) {
    if (cs == "UNSPECIFIED") return spz::CoordinateSystem::UNSPECIFIED;
    if (cs == "LDB") return spz::CoordinateSystem::LDB;
    if (cs == "RDB") return spz::CoordinateSystem::RDB;
    if (cs == "LUB") return spz::CoordinateSystem::LUB;
    if (cs == "RUB") return spz::CoordinateSystem::RUB;
    if (cs == "LDF") return spz::CoordinateSystem::LDF;
    if (cs == "RDF") return spz::CoordinateSystem::RDF;
    if (cs == "LUF") return spz::CoordinateSystem::LUF;
    if (cs == "RUF") return spz::CoordinateSystem::RUF;
    throw std::invalid_argument("Invalid coordinate system: " + cs);
}

// PLY to SPZ conversion function
bool ply_to_spz_impl(const std::string& ply_path, const std::string& spz_path, const std::string& coordinate_system) {
    // Validate coordinate system first (let invalid_argument exceptions propagate)
    spz::CoordinateSystem coord_sys = stringToCoordinateSystem(coordinate_system);
    
    try {
        // Set up unpacking options (for loading PLY)
        spz::UnpackOptions unpack_options;
        unpack_options.to = coord_sys;
        
        // Load the PLY file
        spz::GaussianCloud gaussians = spz::loadSplatFromPly(ply_path, unpack_options);
        
        if (gaussians.numPoints == 0) {
            return false;  // Failed to load PLY file
        }
        
        // Set up packing options (for saving SPZ)
        spz::PackOptions pack_options;
        pack_options.from = coord_sys;
        
        // Save as SPZ
        bool success = spz::saveSpz(gaussians, pack_options, spz_path);
        
        return success;
    } catch (const std::invalid_argument& e) {
        // Re-throw coordinate system errors
        throw;
    } catch (const std::exception& e) {
        // Other errors return false
        return false;
    }
}

// SPZ to PLY conversion function
bool spz_to_ply_impl(const std::string& spz_path, const std::string& ply_path, const std::string& coordinate_system) {
    // Validate coordinate system first (let invalid_argument exceptions propagate)
    spz::CoordinateSystem coord_sys = stringToCoordinateSystem(coordinate_system);
    
    try {
        // Set up unpacking options (for loading SPZ)
        spz::UnpackOptions unpack_options;
        unpack_options.to = coord_sys;
        
        // Load the SPZ file
        spz::GaussianCloud gaussians = spz::loadSpz(spz_path, unpack_options);
        
        if (gaussians.numPoints == 0) {
            return false;  // Failed to load SPZ file
        }
        
        // Set up packing options (for saving PLY)
        spz::PackOptions pack_options;
        pack_options.from = coord_sys;
        
        // Save as PLY
        bool success = spz::saveSplatToPly(gaussians, pack_options, ply_path);
        
        return success;
    } catch (const std::invalid_argument& e) {
        // Re-throw coordinate system errors
        throw;
    } catch (const std::exception& e) {
        // Other errors return false
        return false;
    }
}

PYBIND11_MODULE(_core, m) {
    m.doc() = "SPZ Python bindings for compressed 3D Gaussian splats";
    
    m.def("_ply_to_spz_impl", &ply_to_spz_impl, 
          "Convert PLY file to SPZ format",
          py::arg("ply_path"), py::arg("spz_path"), py::arg("coordinate_system"));
    
    m.def("_spz_to_ply_impl", &spz_to_ply_impl,
          "Convert SPZ file to PLY format", 
          py::arg("spz_path"), py::arg("ply_path"), py::arg("coordinate_system"));
} 