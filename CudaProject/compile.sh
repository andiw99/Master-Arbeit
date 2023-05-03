#!/bin/bash

# Load modules
module load cuda/11.8
module load gcc/9.2.0
module load gcc/10.2.0
module load nvidia/21.7

# Get the input and output file names from the command line arguments
input_file=${1:-"main.cu"}
output_file=${2:-"main_script"}

# Compile the C++ file with g++
nvcc -I /home/weitze73/Code/boost_1_82_0 -std=c++17  $input_file -o $output_file

# Check if the compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful"
else
    echo "Compilation failed"
fi

