#!/bin/bash

# Check if correct number of arguments are provided
if [ "$#" -lt 2 ]; then
    echo "Usage: $0 <start> <end> <n> [template_file]"
    exit 1
fi

dir="Code/Master-Arbeit/CudaProject/parameters/"
# Parse command line arguments
start=$1
end=$2
n=$3
template_file_nr=$4

# Calculate step size
step_size=$(awk "BEGIN {print ($end - $start) / $n}")

# Create output directory if it doesn't exist
mkdir -p "$output_dir"

# Generate and write values to the output file
for ((i=0; i<=$n; i++)); do
    current_value=$(awk "BEGIN {print $start + $i * $step_size}")
    output_file="$output_dir/para_set_nr+${i}.txt"

    # Replace placeholder values in the template file and save to output file
    sed -e "s/min_temp, current_value/min_temp, $current_value/g" -e "s/max_temp, current_value/max_temp, $current_value)/g" "$dir + $template_file" > "$dir + $output_file"

    echo "Values replaced in the template file: $output_file"
done

exit 0
