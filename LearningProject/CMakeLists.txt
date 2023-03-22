cmake_minimum_required(VERSION 3.24)
project(LearningProject)

set(CMAKE_CXX_STANDARD 17)

add_executable(LearningProject main.cpp Lattice.cpp SolarSystem.cpp point_type.hpp cmake-build-debug/matrices.cpp cmake-build-debug/matrices.h cmake-build-debug/kzm_silicon_2d.cpp cmake-build-debug/rng.cpp cmake-build-debug/rng.h)
