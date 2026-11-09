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

namespace spz {

constexpr uint8_t kMaxQuantizedBits = 16;

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

}  // namespace spz

#endif  // SPZ_QUANTIZED_PLANES_H_
