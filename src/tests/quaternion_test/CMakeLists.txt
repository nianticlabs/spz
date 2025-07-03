cmake_minimum_required(VERSION 3.10)
project(quaternion_test)
find_package(spz REQUIRED)

add_executable(quaternion_test quaternion_test.cpp)

set_target_properties(quaternion_test PROPERTIES
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED ON
)

target_link_libraries(quaternion_test PRIVATE spz::spz)
