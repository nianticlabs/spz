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


#include "zstd-batch.h"
#include "load-spz.h"

#include <zstd.h>

#include <algorithm>
#include <atomic>
#include <future>
#include <thread>

namespace spz {
namespace {
constexpr int kCompressionLevel = 12;

bool compressZstd(ByteSpan src, std::vector<uint8_t>* out) {
  const size_t bound = ZSTD_compressBound(src.size);
  out->resize(bound);
  ZSTD_CCtx* cctx = ZSTD_createCCtx();
  if (!cctx) return false;
  size_t result = ZSTD_CCtx_setParameter(cctx, ZSTD_c_compressionLevel, kCompressionLevel);
  if (!ZSTD_isError(result)) {
    result = ZSTD_compress2(cctx, out->data(), bound, src.data, src.size);
  }
  ZSTD_freeCCtx(cctx);
  if (ZSTD_isError(result)) return false;
  out->resize(result);
  return true;
}

// Runs task(0) ... task(count - 1), spread over at most one thread per core.
template <typename Task>
void runParallel(size_t count, const Task& task) {
#if defined(__EMSCRIPTEN__)
  // TODO: Add support for parallel execution on WASM.
  for (size_t i = 0; i < count; i++) {
    task(i);
  }
#else
  const size_t workers =
    std::min<size_t>(count, std::max(1u, std::thread::hardware_concurrency()));
  std::atomic<size_t> next{0};
  std::vector<std::future<void>> futures;
  futures.reserve(workers);
  for (size_t w = 0; w < workers; w++) {
    futures.push_back(std::async(std::launch::async, [&next, count, &task]() {
      for (size_t i = next++; i < count; i = next++) {
        task(i);
      }
    }));
  }
  for (auto& f : futures) {
    f.get();
  }
#endif
}
}  // namespace

bool compressAll(const std::vector<ByteSpan>& srcs, std::vector<std::vector<uint8_t>>* compressed) {
  compressed->assign(srcs.size(), {});
  std::vector<uint8_t> succeeded(srcs.size(), 0);
  runParallel(srcs.size(), [&srcs, compressed, &succeeded](size_t i) {
    succeeded[i] = compressZstd(srcs[i], &(*compressed)[i]);
  });
  if (std::find(succeeded.begin(), succeeded.end(), 0) != succeeded.end()) {
    SpzLog("[SPZ ERROR] compressAll: ZSTD compression failed");
    return false;
  }
  return true;
}

std::vector<uint8_t> decompressAll(const std::vector<DecompressTask>& tasks) {
  std::vector<uint8_t> succeeded(tasks.size(), 0);
  runParallel(tasks.size(), [&tasks, &succeeded](size_t i) {
    const DecompressTask& t = tasks[i];
    const size_t ret = ZSTD_decompress(t.dst, t.dstSize, t.src, t.srcSize);
    succeeded[i] = !ZSTD_isError(ret) && ret == t.dstSize;
  });
  return succeeded;
}

}  // namespace spz
