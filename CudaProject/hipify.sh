#!/bin/bash

# Get the input and output file names from the command line arguments
input_file=${1:-"main.cu"}
output_file=${2:-"main_script"}
inter_file=${3:-"main.cpp"}

# transform cu to cpp
# hipify-perl $input_file > $inter_file
./hipify-clang $input_file --cuda-path=/usr/local/cuda -- -std=c++17

# Compile the c++ file with hipcc
hipcc -I /opt/rocm/include/ -std=c++17  $inter_file -o $output_file

# Check if the compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful"
    # run the simulation
    # save the output in variable output
    output=$(./$output_file | tee /dev/tty)
    # cut only the last line, which is the save path
    rootpath=$(echo "$output" | tail -1)
else
    echo "Compilation failed"
fi
