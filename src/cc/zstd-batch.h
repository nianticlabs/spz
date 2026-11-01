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


// Internal to the library (not installed): compresses and decompresses many independent zstd
// frames in parallel, for the core streams and extension frames alike.

#ifndef SPZ_ZSTD_BATCH_H_
#define SPZ_ZSTD_BATCH_H_

#include <cstddef>
#include <cstdint>
#include <vector>

namespace spz {

// Bytes to compress, owned elsewhere.
struct ByteSpan {
  const uint8_t* data;
  size_t size;
};

// Decompresses src into exactly dstSize bytes at dst.
struct DecompressTask {
  const uint8_t* src;
  size_t srcSize;
  uint8_t* dst;
  size_t dstSize;
};

// Compresses each span into its own zstd frame, in parallel; compressed[i] holds srcs[i]. Returns
// false (after logging) if any compression fails. Precondition: no span is empty.
bool compressAll(const std::vector<ByteSpan>& srcs, std::vector<std::vector<uint8_t>>* compressed);

// Decompresses every task in parallel. Returns one flag per task, nonzero when the task produced
// exactly its dstSize bytes.
std::vector<uint8_t> decompressAll(const std::vector<DecompressTask>& tasks);

}  // namespace spz

#endif  // SPZ_ZSTD_BATCH_H_
