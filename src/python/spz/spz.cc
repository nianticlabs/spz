// Nanobind Python bindings for the spz library.

#include <nanobind/nanobind.h>
#include <nanobind/stl/string.h>
#include <nanobind/ndarray.h>

#include <vector>
#include <cstring>
#include <stdexcept>
#include <string>

#include "src/cc/load-spz.h"
#include "src/cc/splat-types.h"

namespace nb = nanobind;

// Strict 1‑D, C‑contiguous, CPU, float32 ndarray alias
// This enforces that only float32 numpy arrays are accepted, though nanobind
// will automatically convert compatible numeric types (int, float64, etc.) to float32.
// Non-numeric types (strings, objects, etc.) will be rejected at binding time.
using NdArray1D = nb::ndarray<nb::numpy, float, nb::ndim<1>, nb::c_contig, nb::device::cpu>;

// -----------------------------------------------------------------------------
// Helpers
// -----------------------------------------------------------------------------

// Ensure an array length is a multiple of k; throw ValueError with a helpful message otherwise.
static inline void ensure_multiple(const char *name, size_t size, size_t k) {
    if (k == 0) {
        throw nb::value_error("internal error: divisor cannot be zero");
    }
    if (size % k != 0) {
        std::string msg = std::string(name) + " length must be a multiple of " + std::to_string(k) +
                          ", got " + std::to_string(size);
        throw nb::value_error(msg.c_str());
    }
}

// ──────────────────────────────────────────────────────────────────────────────
// Convert a `const std::vector<float>&` into a *new* NumPy float32 1-D array.
// The array owns its own copy of the data, so Python can out-live the C++ vector.
// 
// This lambda is used as a getter for properties like positions, scales, etc.
// It creates a new numpy array that owns its own copy of the data, ensuring
// memory safety when the original C++ object is destroyed.
// ──────────────────────────────────────────────────────────────────────────────
static auto vector_getter = [](const std::vector<float>& vec) -> NdArray1D {
    // If the source vector is empty, hand back an empty ndarray that owns nothing.
    if (vec.empty())
        return NdArray1D(/*data=*/nullptr,
                         /*shape=*/{0},
                         /*owner capsule=*/nb::handle());

    //----------------------------------------------------------------------
    // 1. Allocate a raw buffer and copy the vector’s bytes into it
    //----------------------------------------------------------------------
    const size_t n = vec.size();          // element count
    float* data    = new float[n];        // plain new[] so we can later delete[]
    std::memcpy(data, vec.data(), n * sizeof(float));

    //----------------------------------------------------------------------
    // 2. Create a nanobind capsule that *knows how to free that buffer*
    //    when Python’s ref-count for the ndarray reaches zero.
    //----------------------------------------------------------------------
    nb::capsule owner(
        data,
        [](void* p) noexcept {            // custom deleter (must *not* throw)
            delete[] static_cast<float*>(p);
        });

    //----------------------------------------------------------------------
    // 3. Build and return a 1-D ndarray that:
    //      • points at `data`
    //      • sees its length as {n}
    //      • holds a reference to `owner`, thereby tying the buffer’s lifetime
    //        to the ndarray’s lifetime on the Python side.
    //----------------------------------------------------------------------
    return NdArray1D(data, {n}, owner);
};

// ──────────────────────────────────────────────────────────────────────────────
// Copy data from a NumPy array back into a std::vector<float>&.
// Expands or shrinks the vector as needed.
// 
// This lambda is used as a setter for properties like positions, scales, etc.
// It copies the data from the numpy array into the C++ vector. The nanobind
// NdArray1D type alias already ensures type safety - only float32 arrays are
// accepted (though nanobind will convert compatible numeric types automatically).
// Non-numeric types are rejected at the binding level.
// ──────────────────────────────────────────────────────────────────────────────
static auto vector_setter = [](std::vector<float>& vec,
                               const NdArray1D&   arr) {
    //----------------------------------------------------------------------
    // Resize the destination vector to exactly match the ndarray length.
    //----------------------------------------------------------------------
    const size_t n = arr.shape(0);
    vec.resize(n);

    //----------------------------------------------------------------------
    // Copy bytes. If n == 0 the memcpy is skipped to avoid UB.
    //----------------------------------------------------------------------
    if (n)
        std::memcpy(vec.data(), arr.data(), n * sizeof(float));
};

// -----------------------------------------------------------------------------
// Python module
// -----------------------------------------------------------------------------

NB_MODULE(spz, m) {
    m.doc() = "Python bindings for the spz library (Gaussian splatting).";

    // -------------------------------------------------------------------------
    // Enums
    // Coordinate system enum with comprehensive documentation
    nb::enum_<spz::CoordinateSystem>(m, "CoordinateSystem", R"doc(
        Coordinate system conventions for 3D Gaussian splats.
        
        These enums define the handedness and axis orientations used in different 3D graphics
        systems and file formats. Each value encodes three axes directions:
        - First letter: Right (R) or Left (L) - defines the X-axis direction
        - Second letter: Up (U) or Down (D) - defines the Y-axis direction  
        - Third letter: Front (F) or Back (B) - defines the Z-axis direction
        
        Common usage patterns:
        - RDF: Standard PLY file format (right-handed, Y-down, Z-forward)
        - RUB: Three.js coordinate system (right-handed, Y-up, Z-backward)
        - LUF: glTF/GLB format (left-handed, Y-up, Z-forward)
        - RUF: Unity game engine (left-handed, Y-up, Z-forward)
        
        Use convert_coordinates() to transform between different coordinate systems.
    )doc")
        .value("UNSPECIFIED", spz::CoordinateSystem::UNSPECIFIED, "Unspecified coordinate system")
        .value("LDB",         spz::CoordinateSystem::LDB, "Left Down Back - left-handed, Y-down, Z-backward")
        .value("RDB",         spz::CoordinateSystem::RDB, "Right Down Back - right-handed, Y-down, Z-backward")
        .value("LUB",         spz::CoordinateSystem::LUB, "Left Up Back - left-handed, Y-up, Z-backward")
        .value("RUB",         spz::CoordinateSystem::RUB, "Right Up Back - right-handed, Y-up, Z-backward (Three.js)")
        .value("LDF",         spz::CoordinateSystem::LDF, "Left Down Front - left-handed, Y-down, Z-forward")
        .value("RDF",         spz::CoordinateSystem::RDF, "Right Down Front - right-handed, Y-down, Z-forward (PLY format)")
        .value("LUF",         spz::CoordinateSystem::LUF, "Left Up Front - left-handed, Y-up, Z-forward (GLB format)")
        .value("RUF",         spz::CoordinateSystem::RUF, "Right Up Front - right-handed, Y-up, Z-forward (Unity)")
        .export_values();

    // -------------------------------------------------------------------------
    // Options structs
    // -------------------------------------------------------------------------
    nb::class_<spz::PackOptions>(m, "PackOptions")
        .def(nb::init<>())
        .def_rw("from_coord", &spz::PackOptions::from,
                "Coordinate system of the input splat");

    nb::class_<spz::UnpackOptions>(m, "UnpackOptions")
        .def(nb::init<>())
        .def_rw("to_coord", &spz::UnpackOptions::to,
                "Desired coordinate system of the output splat");

    // -------------------------------------------------------------------------
    // GaussianCloud - Main data structure for 3D Gaussian splats
    // -------------------------------------------------------------------------
    // This class represents a collection of 3D Gaussian splats, each defined by:
    // - Position (xyz coordinates)
    // - Scale (log-scale size factors for x,y,z axes)  
    // - Rotation (quaternion for orientation)
    // - Alpha (opacity before sigmoid activation)
    // - Color (base RGB values)
    // - Spherical Harmonics (view-dependent color coefficients)
    //
    // The property bindings use custom lambdas (vector_getter/vector_setter) to:
    // 1. Convert between C++ std::vector<float> and Python numpy arrays
    // 2. Ensure memory safety by copying data (not sharing pointers)
    // 3. Validate array types (enforces float32 dtype)
    // 4. Handle empty arrays gracefully
    // -------------------------------------------------------------------------
    nb::class_<spz::GaussianCloud> cloud(m, "GaussianCloud");
    cloud
        .def(nb::init<>(), "Construct an empty GaussianCloud.")
        .def_prop_ro("num_points",
                     [](const spz::GaussianCloud &self) {
                         return static_cast<int32_t>(self.positions.size() / 3);
                     },
                     "Total number of points (gaussians) in this splat (derived from positions).")
        .def("__len__", [](const spz::GaussianCloud &self) {
            return static_cast<int32_t>(self.positions.size() / 3);
        })
        .def("__repr__", [](const spz::GaussianCloud &self) {
            const int32_t n = static_cast<int32_t>(self.positions.size() / 3);
            return nb::str("GaussianCloud(num_points={}, sh_degree={}, antialiased={})")
                .format(n, self.shDegree, self.antialiased);
        })
        .def_prop_rw("sh_degree",
                     [](const spz::GaussianCloud &self) { return self.shDegree; },
                     [](spz::GaussianCloud &self, int32_t deg) {
                         if (deg < 0 || deg > 3) {
                             throw nb::value_error("sh_degree must be in [0, 3]");
                         }
                         self.shDegree = deg;
                     },
                     "Degree of spherical harmonics for this splat (restricted to 0..3).")
        .def_rw("antialiased", &spz::GaussianCloud::antialiased, "Whether the gaussians should be rendered in antialiased mode (mip splatting).")
        // Bind the 'positions' property as a read-write numpy array
        // - Getter: Creates a new numpy array copy of the positions vector
        // - Setter: Copies data from numpy array into the positions vector
        // - nb::rv_policy::move: Optimizes return value handling
        // - Array contains 3D coordinates as [x1,y1,z1,x2,y2,z2,...] 
        // - Accepts numeric arrays (int, float64, etc.) but converts to float32
        .def_prop_rw("positions",
                     [](const spz::GaussianCloud &self){ return vector_getter(self.positions); },
                     [](spz::GaussianCloud &self, const NdArray1D &arr){
                         const size_t n = arr.shape(0);
                         ensure_multiple("positions", n, 3);
                         vector_setter(self.positions, arr);
                         self.numPoints = static_cast<int32_t>(self.positions.size() / 3);
                     },
                     nb::rv_policy::move, R"doc(Gaussian centers as a flat xyz array.

                - dtype: any numeric (auto-converts to float32)
                - shape: (num_points * 3,), must be a multiple of 3
                - semantics: setting positions defines `num_points = len(positions) / 3`
                - layout: [x1, y1, z1, x2, y2, z2, ...])doc")
        // Bind the 'scales' property - log-scale values for Gaussian ellipsoid sizes
        // Array contains 3D scale factors as [sx1,sy1,sz1,sx2,sy2,sz2,...]
        .def_prop_rw("scales",
                     [](const spz::GaussianCloud &self){ return vector_getter(self.scales); },
                     [](spz::GaussianCloud &self, const NdArray1D &arr){
                         const size_t n = arr.shape(0);
                         ensure_multiple("scales", n, 3);
                         vector_setter(self.scales, arr);
                         if (self.numPoints > 0 && self.scales.size() != static_cast<size_t>(self.numPoints) * 3) {
                             throw nb::value_error("scales length must equal num_points * 3");
                         }
                     },
                     nb::rv_policy::move, R"doc(Log-scale ellipsoid radii per axis.

                - dtype: numeric → float32
                - shape: (num_points * 3,), must be a multiple of 3 and equal to num_points * 3 when num_points > 0
                - layout: [sx1, sy1, sz1, sx2, sy2, sz2, ...]
                - note: stored on log scale; actual radius = exp(scale))doc")
        // Bind the 'rotations' property - quaternions for Gaussian orientations
        // Array contains quaternions as [x1,y1,z1,w1,x2,y2,z2,w2,...]
        .def_prop_rw("rotations",
                     [](const spz::GaussianCloud &self){ return vector_getter(self.rotations); },
                     [](spz::GaussianCloud &self, const NdArray1D &arr){
                         const size_t n = arr.shape(0);
                         ensure_multiple("rotations", n, 4);
                         vector_setter(self.rotations, arr);
                         if (self.numPoints > 0 && self.rotations.size() != static_cast<size_t>(self.numPoints) * 4) {
                             throw nb::value_error("rotations length must equal num_points * 4");
                         }
                     },
                     nb::rv_policy::move, R"doc(Orientation quaternions, xyzw per point.

                - dtype: numeric → float32
                - shape: (num_points * 4,), must be a multiple of 4 and equal to num_points * 4 when num_points > 0
                - layout: [x1, y1, z1, w1, x2, y2, z2, w2, ...])doc")
        // Bind the 'alphas' property - opacity values (before sigmoid activation)
        // Array contains one float per Gaussian: [a1,a2,a3,...]
        .def_prop_rw("alphas",
                     [](const spz::GaussianCloud &self){ return vector_getter(self.alphas); },
                     [](spz::GaussianCloud &self, const NdArray1D &arr){
                         const size_t n = arr.shape(0);
                         // multiple of 1 is always true; enforce against num_points when set
                         vector_setter(self.alphas, arr);
                         if (self.numPoints > 0 && self.alphas.size() != static_cast<size_t>(self.numPoints)) {
                             throw nb::value_error("alphas length must equal num_points");
                         }
                     },
                     nb::rv_policy::move, R"doc(Opacity values before sigmoid activation.

                - dtype: numeric → float32
                - shape: (num_points,), must equal num_points when num_points > 0
                - note: effective opacity is sigmoid(alpha))doc")
        // Bind the 'colors' property - base RGB colors (DC component of spherical harmonics)
        // Array contains RGB triplets as [r1,g1,b1,r2,g2,b2,...]
        .def_prop_rw("colors",
                     [](const spz::GaussianCloud &self){ return vector_getter(self.colors); },
                     [](spz::GaussianCloud &self, const NdArray1D &arr){
                         const size_t n = arr.shape(0);
                         ensure_multiple("colors", n, 3);
                         vector_setter(self.colors, arr);
                         if (self.numPoints > 0 && self.colors.size() != static_cast<size_t>(self.numPoints) * 3) {
                             throw nb::value_error("colors length must equal num_points * 3");
                         }
                     },
                     nb::rv_policy::move, R"doc(Base RGB (SH DC component), triplets per point.

                - dtype: numeric → float32
                - shape: (num_points * 3,), must be a multiple of 3 and equal to num_points * 3 when num_points > 0
                - layout: [r1, g1, b1, r2, g2, b2, ...])doc")
        // Bind the 'sh' property - spherical harmonics coefficients for view-dependent colors
        // Array layout is complex - see detailed documentation in the docstring below
        .def_prop_rw("sh",
                     [](const spz::GaussianCloud &self){ return vector_getter(self.sh); },
                     [](spz::GaussianCloud &self, const NdArray1D &arr){
                         const size_t n = arr.shape(0);
                         ensure_multiple("sh", n, 3);
                         // Compute dim per point per channel based on degree: 0, 3, 8, 15
                         int deg = self.shDegree;
                         int shDimPerChannel = (deg == 0) ? 0 : ((deg + 1) * (deg + 1) - 1);
                         if (shDimPerChannel == 0) {
                             if (n != 0) {
                                 throw nb::value_error("sh must be empty when sh_degree == 0");
                             }
                         } else {
                             ensure_multiple("sh", n, static_cast<size_t>(shDimPerChannel) * 3);
                         }
                         vector_setter(self.sh, arr);
                         if (self.numPoints > 0) {
                             size_t expected = static_cast<size_t>(self.numPoints) * static_cast<size_t>(shDimPerChannel) * 3;
                             if (self.sh.size() != expected) {
                                 throw nb::value_error("sh length must equal num_points * ((sh_degree+1)^2 - 1) * 3");
                             }
                         }
                     },
                     nb::rv_policy::move, R"doc(Spherical harmonics coefficients.
                The number of coefficients per point depends on shDegree:
                    0 -> 0
                    1 -> 9   (3 coeffs x 3 channels)
                    2 -> 24  (8 coeffs x 3 channels)
                    3 -> 45  (15 coeffs x 3 channels)
                    The color channel is the inner (fastest varying) axis, and the coefficient is the outer
                    (slower varying) axis, i.e. for degree 1, the order of the 9 values is:
                    sh1n1_r, sh1n1_g, sh1n1_b, sh10_r, sh10_g, sh10_b, sh1p1_r, sh1p1_g, sh1p1_b.
                    Shapes and constraints:
                    - dtype: numeric → float32
                    - if sh_degree == 0: sh must be empty
                    - else: len(sh) must be a multiple of (((sh_degree+1)^2 - 1) * 3)
                    - if num_points > 0: len(sh) must equal num_points * (((sh_degree+1)^2 - 1) * 3)
                    Tip: set sh_degree first, then provide sh with the exact length.)doc")
        // Methods
        .def("convert_coordinates", &spz::GaussianCloud::convertCoordinates,
             nb::arg("from_coord"), nb::arg("to_coord"),
             R"doc(Convert between two coordinate systems in-place.
                  Example:
                  convert_coordinates(from_coord=spz.CoordinateSystem.RDF,
                                      to_coord=spz.CoordinateSystem.RUB))doc")
        .def("rotate_180_deg_about_x", &spz::GaussianCloud::rotate180DegAboutX,
             R"doc(Convenience for RUB ↔ RDF conversion (180° about X).

                 Applies the same internal logic as convert_coordinates(RUB, RDF).)doc")
        .def("median_volume", &spz::GaussianCloud::medianVolume,
             "Return the median Gaussian volume.");

    // -------------------------------------------------------------------------
    // Functions
    // -------------------------------------------------------------------------

    m.def("load_spz", (spz::GaussianCloud (*)(const std::string &, const spz::UnpackOptions &)) &spz::loadSpz,
          nb::arg("filename"), nb::arg("options") = spz::UnpackOptions(),
          "Load a *.spz* file and return a GaussianCloud.");

    m.def("save_spz", (bool (*)(const spz::GaussianCloud &, const spz::PackOptions &, const std::string &)) &spz::saveSpz,
          nb::arg("gaussians"), nb::arg("options"), nb::arg("filename"),
          "Save a GaussianCloud to a *.spz* file.");

    m.def("load_splat_from_ply", (spz::GaussianCloud (*)(const std::string &, const spz::UnpackOptions &)) &spz::loadSplatFromPly,
          nb::arg("filename"), nb::arg("options") = spz::UnpackOptions(),
          "Read GaussianCloud data from a *.ply* file.");

    m.def("save_splat_to_ply", (bool (*)(const spz::GaussianCloud &, const spz::PackOptions &, const std::string &)) &spz::saveSplatToPly,
          nb::arg("gaussians"), nb::arg("options"), nb::arg("filename"),
          "Write GaussianCloud data to a *.ply* file.");
}
