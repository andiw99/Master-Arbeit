cmake_minimum_required(VERSION 3.25)

set(CMAKE_CUDA_ARCHITECTURES 75)
set(CMAKE_CUDA_COMPILER /usr/local/cuda-12.1/bin/nvcc)
project(CudaProject CUDA)

set(CMAKE_CUDA_STANDARD 23)


add_executable(CudaProject
        Simulation-cuda.cuh
        SubsystemRelaxation-cuda.cu
        main-cuda.cuh
        systems-cuda.cuh
        observers.cuh
        steppers.cuh
        SubsystemQuench.cu
        AutoQuench-h.cu
        AutoQuench.cu
)

set_target_properties(CudaProject PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON)
