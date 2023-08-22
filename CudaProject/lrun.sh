module load gcc/10.2.0
module load cuda/12.1
module load openmpi/4.1.1
module load fftw/3.3.10-ompi411
module load boost
# Get the input and output file names from the command line arguments
input_file=${1:-"main.cu"}
output_file=${2:-"main_script"}

# Compile the cuda file with nvcc
nvcc -I /home/weitze73/Code/boost_1_82_0 -std=c++17  $input_file -o $output_file
# Check if the compilation was successful
if [ $? -eq 0 ]; then
    echo "Compilation successful"
    # run the simulation
    # save the output in variable output
    ./${output_file}
    # cut only the last line, which is the save path
else
    echo "Compilation failed"
fi
