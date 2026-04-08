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

#ifndef SPZ_EMSCRIPTEN_UTILS_H_
#define SPZ_EMSCRIPTEN_UTILS_H_

#include <emscripten/val.h>
#include <type_traits>
#include <vector>

// Utility functions for converting between C++ vectors and JavaScript arrays
template <typename T>
inline ::emscripten::val jsArrayFromVector(const std::vector<T>& vec) {
  ::emscripten::val array = ::emscripten::val::array();
  for (const auto& item : vec) {
    array.call<void>("push", item);
  }
  return array;
}

template <typename T>
inline void vectorFromJsArray(const ::emscripten::val& array, std::vector<T>& out) {
  if constexpr (std::is_arithmetic_v<T>) {
    out = ::emscripten::convertJSArrayToNumberVector<T>(array);
  } else {
    const size_t length = array["length"].as<size_t>();
    out.resize(length);
    for (size_t i = 0; i < length; ++i) {
      out[i] = array[i].as<T>();
    }
  }
}

#endif  // SPZ_EMSCRIPTEN_UTILS_H_

