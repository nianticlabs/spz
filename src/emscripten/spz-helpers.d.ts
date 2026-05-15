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

import type {
  SpzModule,
  SpzStreamCallbacks,
  SpzStreamHeader,
  UnpackOptions
} from './spz';

/**
 * Runtime mirror of the wasm `SplatAttribute` enum. The values match the
 * integer codes delivered to streaming `onChunk` callbacks, so you can do
 * `if (attr === SplatAttribute.Sh) { ... }` against the numeric `attr`.
 */
export const SplatAttribute: {
  readonly Positions: 0;
  readonly Alphas:    1;
  readonly Colors:    2;
  readonly Scales:    3;
  readonly Rotations: 4;
  readonly Sh:        5;
};

/**
 * Copy `bytes` into a wasm-heap region, run `fn(ptr, byteLength)`, and free
 * the region (even if `fn` throws). Use for direct calls into pointer-taking
 * bindings like `loadSpzFromHeapPtr` when you need finer control than
 * `loadSpzStreamingAsync` provides.
 */
export function withHeapBuffer<T>(
  spz: SpzModule,
  bytes: Uint8Array,
  fn: (ptr: number, byteLength: number) => T
): T;

/**
 * Promise wrapper around `loadSpzStreaming`. Copies `bytes` into the wasm
 * heap, drives the streaming decode with the caller's `onChunk` (and optional
 * `onHeader`), resolves with the header on `onDone`, and rejects on
 * `onError` or any thrown exception. The heap buffer is freed in all paths.
 */
export function loadSpzStreamingAsync(
  spz: SpzModule,
  bytes: Uint8Array,
  options: UnpackOptions,
  onChunk: SpzStreamCallbacks['onChunk'],
  onHeader?: SpzStreamCallbacks['onHeader']
): Promise<SpzStreamHeader>;
