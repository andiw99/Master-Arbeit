#!/bin/bash


# Loading modules not necessary on jupyter and somehow not working with ssh
# module load cuda/11.8
# module load gcc/9.2.0
# module load gcc/10.2.0
# module load nvidia/21.7
module load openmpi/4.1.1
module load fftw/3.3.10-ompi411

# Get the input and output file names from the command line arguments
input_file=${1:-"main.cu"}
output_file=${2:-"main_script"}

# Compile the cuda file with nvcc
nvcc -I /home/weitze73/Code/boost_1_82_0 -std=c++17  $input_file -o $output_file
# compile everything else for calcing correlation function and stuff
g++ -I /home/weitze73/Code/boost_1_82_0 -std=c++17 ./../LearningProject/DomainSize/LatticeFFTError.cpp -o fourierlatticeTraf.out -lfftw3 -lm
# Check if the compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful"
    # run the simulation
    # save the output in variable output
    output=$(./$output_file | tee /dev/tty)
    # cut only the last line, which is the save path
    rootpath=$(echo "$output" | tail -1)
    # run calc correlation function with the path as input variable
    # remember to call it with argument since otherwise the hardcoded path is used
    ./fourierlatticeTraf.out "need to call with argument"
    # make the plots
    python3 ./../PythonProject/PlotBathGPU.py "$rootpath"
else
    echo "Compilation failed"
fi

