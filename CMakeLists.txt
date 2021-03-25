cmake_minimum_required(VERSION 3.10)
project(heat)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "-pthread -O3")


add_executable(heat Temperature_map.cpp Temperature_map.h main.cpp)
