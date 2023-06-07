#!/bin/bash

# Get the input and output file names from the command line arguments
input_file=${1:-"main.cu"}
output_file=${2:-"main_script"}

# Compile the C++ file with g++
g++ -I /home/weitze73/Code/boost_1_82_0 -std=c++17  $input_file -o $output_file

# Check if the compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful, starting..."
    output=$(./$output_file test)
    rootpath=$(echo "$output" | tail -1) 
    echo "$rootpath"
    
    output=$(./$output_file "$rootpath")
    rootpath=$(echo "$output" | tail -1) 
    echo "$rootpath"
else
    echo "Compilation failed"
fi


