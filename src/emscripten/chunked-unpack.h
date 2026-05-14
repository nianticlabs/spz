/*
MIT License

Copyright (c) 2026 Niantic Labs
Copyright (c) 2026 Adobe Inc.

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

#pragma once

#include <cstdint>

#include "load-spz.h"
#include "splat-types.h"

namespace spz {

// Unpack `pointCount` points starting at `pointOffset` of `attr` into `out`.
// Must provide at least `pointCount * floatsPerPoint(attr, shDegree)` floats.
//
// Returns false on validation failure (range error, attribute already released,
// null `out`). Not thread-safe on the same PackedGaussians.
bool unpackChunk(PackedGaussians &packed,
                 SplatAttribute attr,
                 int32_t pointOffset,
                 int32_t pointCount,
                 float *out,
                 CoordinateSystem coordTo);

}  // namespace spz
