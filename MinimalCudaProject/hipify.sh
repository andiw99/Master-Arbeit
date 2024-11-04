#!/bin/bash

# Get the input and output file names from the command line arguments
input_file=${1:-"main.cu"}
output_file=${2:-"main_script"}
para_set=${3:-"0"}
inter_file=${4:-"main.cpp"}

# transform cu to cpp
hipify-perl systems-cuda.cuh > systems.cuh
hipify-perl main-cuda.cuh > main.cuh
# hipify-perl deprecated-systems-cuda.cuh > deprecated-systems.cuh
echo pearled
hipify-perl $input_file > $inter_file
# ./hipify-clang $input_file --cuda-path=/usr/local/cuda -- -std=c++17

# Compile the c++ file with hipcc
hipcc -I /opt/rocm/include/hiprand -I /opt/rocm/include/hipfft/ -std=c++17  $inter_file -o $output_file -lhipfft

# Check if the compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful"
    # run the simulation
    # save the output in variable output
    output=$(./$output_file ${para_set} | tee /dev/tty)
    # cut only the last line, which is the save path
    rootpath=$(echo "$output" | tail -1)
else
    echo "Compilation failed"
fi
