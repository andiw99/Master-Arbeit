cmake_minimum_required(VERSION 3.25)
project(TestCudaProject CUDA)

set(CMAKE_CUDA_STANDARD 23)

add_executable(TestCudaProject main.cu)

set_target_properties(TestCudaProject PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
