/*
MIT License

Copyright (c) 2026 Niantic Spatial

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

// Internal to the library (not installed): linear quantization of per-point float components to
// 1-16 bit levels, stored as byte planes. For extensions that store quantized per-point data.

#ifndef SPZ_QUANTIZED_PLANES_H_
#define SPZ_QUANTIZED_PLANES_H_

#include <cstddef>
#include <cstdint>
#include <optional>
#include <vector>

namespace spz {

constexpr uint8_t kMaxQuantizedBits = 16;
// bits value meaning lossless float32 instead of quantized levels.
constexpr uint8_t kFloat32Bits = 32;

inline bool isQuantizedBits(uint8_t bits) { return bits >= 1 && bits <= kMaxQuantizedBits; }

// Bytes per level, and so the number of byte planes.
size_t quantizedLevelBytes(uint8_t bits);

// Quantizes values[p * stride] for p < points to levels
// round((v - min) / (max - min) * (2^bits - 1)), and writes them to out as
// quantizedLevelBytes(bits) planes of points bytes, least significant first.
// Preconditions: isQuantizedBits(bits); min <= max, both finite; out holds
// points * quantizedLevelBytes(bits) bytes.
void encodeQuantizedPlanes(const float* values, size_t stride, size_t points, float min, float max,
                           uint8_t bits, uint8_t* out);

// Inverse of encodeQuantizedPlanes, writing values[p * stride]. Returns false if a level exceeds
// 2^bits - 1. Preconditions as for encodeQuantizedPlanes, with in holding the planes.
bool decodeQuantizedPlanes(const uint8_t* in, size_t points, float min, float max, uint8_t bits,
                           float* values, size_t stride);

// ---- Component-major per-point data, as custom attributes and temporal channels store it ----
// Point-major values (dim floats per point) are stored component-major: every point's first
// component, then every point's second, and so on. float32 components (bits 32) take 4 bytes,
// little-endian; quantized components (bits 1-16) are byte planes between per-component min/max.

// Bytes the component-major data of points points takes, or std::nullopt if it can't be addressed
// on this platform (e.g. on wasm32). Precondition: bits is 32 or isQuantizedBits(bits).
std::optional<size_t> componentMajorBytes(size_t points, uint16_t dim, uint8_t bits);

// Each component's min and max over point-major values (zeros when there are no points).
void componentRanges(const std::vector<float>& values, uint16_t dim, std::vector<float>* min,
                     std::vector<float>* max);

// Encodes point-major values. Preconditions: values holds points * dim floats; when quantized,
// min and max hold dim finite bounds with min <= max.
std::vector<uint8_t> encodeComponentMajor(const std::vector<float>& values, size_t points,
                                          uint16_t dim, uint8_t bits, const std::vector<float>& min,
                                          const std::vector<float>& max);

// Decodes bytes into point-major values. Returns false if a quantized level exceeds 2^bits - 1.
// Preconditions: bytes holds componentMajorBytes(points, dim, bits) bytes; when quantized, min and
// max hold dim bounds.
bool decodeComponentMajor(const uint8_t* bytes, size_t points, uint16_t dim, uint8_t bits,
                          const std::vector<float>& min, const std::vector<float>& max,
                          std::vector<float>* values);

// ---- Fixed-step quantization ----
// Each component has a step; a value is stored as the integer level round(v / step), zigzag-coded
// to a u32 and split into 4 byte planes (least significant first), component-major. The error is
// at most step / 2 whatever the value, unlike min/max ranges, which one outlier coarsens for every
// point. Small levels leave the high planes zero, which zstd compresses to almost nothing.

// Bytes per component level in fixed-step data.
constexpr size_t kStepLevelBytes = 4;

// Whether every step is finite and positive.
bool validSteps(const std::vector<float>& steps);

// Encodes point-major values into *out. Returns false if a value is non-finite or its level doesn't
// fit in 32 bits. Preconditions: values holds points * dim floats; steps holds dim valid steps.
bool encodeComponentMajorSteps(const std::vector<float>& values, size_t points, uint16_t dim,
                               const std::vector<float>& steps, std::vector<uint8_t>* out);

// Decodes bytes into point-major values: level * step, computed in double precision. Preconditions:
// bytes holds points * dim * kStepLevelBytes bytes; steps holds dim steps.
void decodeComponentMajorSteps(const uint8_t* bytes, size_t points, uint16_t dim,
                               const std::vector<float>& steps, std::vector<float>* values);

// The value fixed-step coding with this step stores for v (what decoding returns). Precondition:
// the level fits in 32 bits.
float quantizeToStep(float v, float step);

}  // namespace spz

#endif  // SPZ_QUANTIZED_PLANES_H_
