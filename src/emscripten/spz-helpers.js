/*
MIT License

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

// JS-side helpers that wrap the raw wasm streaming APIs with the plumbing
// every consumer ends up writing: heap allocation around an input buffer,
// Promise wrapping around the synchronous chunk callbacks, and error
// propagation. The helpers stay layout-agnostic — the per-attribute demux
// (e.g. interleaved-vs-per-coordinate columns, quaternion ordering, SH
// layout) is the consumer's responsibility.

// Runtime mirror of the wasm SplatAttribute enum (registered with
// enum_value_type::number, so the values match what arrives in onChunk).
export const SplatAttribute = Object.freeze({
  Positions: 0,
  Alphas:    1,
  Colors:    2,
  Scales:    3,
  Rotations: 4,
  Sh:        5
});

// Copy `bytes` into a freshly allocated WASM-heap region, run `fn(ptr, len)`,
// and free the region in `finally` so the heap doesn't leak on throw. Use this
// when you need to call a heap-ptr binding (e.g. loadSpzFromHeapPtr) directly.
export const withHeapBuffer = (spz, bytes, fn) => {
  const len = bytes.byteLength;
  const ptr = spz._malloc(len);
  if (ptr === 0) {
    throw new Error(`SPZ heap allocation failed (${len} bytes)`);
  }
  spz.HEAPU8.set(bytes, ptr);
  try {
    return fn(ptr, len);
  } finally {
    spz._free(ptr);
  }
};

// Promise wrapper for `loadSpzStreaming`. Copies `bytes` into the WASM heap,
// runs the streaming decode with the caller's `onChunk` (and optional
// `onHeader`), resolves with the header on completion, rejects on any error.
// Frees the heap buffer in all paths.
export const loadSpzStreamingAsync = (spz, bytes, options, onChunk, onHeader) =>
  new Promise((resolve, reject) => {
    let header;
    try {
      withHeapBuffer(spz, bytes, (ptr, len) => {
        spz.loadSpzStreaming(ptr, len, options, {
          onHeader: (info) => {
            header = info;
            if (onHeader) onHeader(info);
          },
          onChunk,
          onDone: () => resolve(header),
          onError: (message) => reject(new Error(message))
        });
      });
    } catch (err) {
      reject(err);
    }
  });
