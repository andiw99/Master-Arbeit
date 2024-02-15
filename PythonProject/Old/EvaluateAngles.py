import math
import numpy as np
import os

def calculate_angle_with_y_axis(point1, point2):
    # Calculate the angle with the y-axis using arctan
    dz = (point2[2] - point1[2])
    dy = (point2[1] - point1[1])
    angle_rad = math.atan(dz / dy)
    angle_deg = math.degrees(angle_rad)
    return angle_deg

def extract_top_n_si_lines(file_path, n):
    si_lines = []

    with open(file_path, 'r') as file:
        for line in file:
            if line.startswith("Si"):
                parts = line.split()
                z_value = float(parts[-1])
                si_lines.append((line, z_value))

    # Sort si_lines based on Z values (last column)
    si_lines.sort(key=lambda x: x[1], reverse=True)

    # Extract the top n lines with the largest Z values
    top_n_si_lines = [line for (line, z_value) in si_lines[:n]]

    return top_n_si_lines
def average_angles_over_lines(file_path, x_tolerance=0.2, y_tolerance=3, nr_dimers=12):
    # Read coordinates from file
    coordinates = []
    with open(file_path, 'r') as file:
        largest_z_lines = extract_top_n_si_lines(file_path, 2 * nr_dimers)
        for line in largest_z_lines:
            if line.startswith("Si"):
                parts = line.split()
                coordinates.append((float(parts[1]), float(parts[2]), float(parts[3])))

    # Group coordinates based on x-values
    grouped_coordinates = []
    lonely_coordinates = []
    for x, y, z in coordinates:
        newly_found = None
        for x_lonely, y_lonely, z_lonely in lonely_coordinates:
            if abs(x - x_lonely) < x_tolerance and abs(y - y_lonely) < y_tolerance:
                # if this is true, the coorinates belong together
                grouped_coordinates.append(((x, y, z), (x_lonely, y_lonely, z_lonely)))
                newly_found = (x_lonely, y_lonely, z_lonely)
                break
        if newly_found:
            lonely_coordinates.remove(newly_found)
        else:
            # if nothing has been found, the current coordinates are lonely
            lonely_coordinates.append((x, y, z))


        # now we somehow need to remove
    # Calculate angles with y-axis for each group of coordinates
    angles = []
    angle_map = {}
    for si1, si2 in grouped_coordinates:
            # Calculate angle for each pair of consecutive points in the group
        angle = calculate_angle_with_y_axis(si1, si2)
        angle_map[(si1, si2)] = angle
        angles.append(angle)
    # sort the angle map for x values
    angle_map = {k: v for k, v in sorted(angle_map.items(), key=lambda item:item[0][0][0])}


    # Average angles over all lines
    if angles:
        print(angles)
        print(len(angle_map))
        for k, v in angle_map.items():
            print(k, v, "\t", 90 - np.abs(v), "\t", (90 - np.abs(v)) / 180 * np.pi)
        average_angle = sum(np.abs(angles)) / len(angles)
        return average_angle
    else:
        return None

def evaluate_file(file_path):
    average_angle = average_angles_over_lines(file_path)
    if average_angle is not None:
        print("Average angle with y-axis:", average_angle)
        theta = 90 - average_angle
        print("Average theta:", theta)
        print("theta_rad:" , theta / 180 * np.pi)
    else:
        print("No data points found to calculate the angle.")

# Example usage:
def main():
    path = "../../../Generated content/DFT/p1-si-relaxed-geom-pbesol-pbe-pz/pbe/"  # Replace with the path to your file
    if os.path.isdir(path):
        # if this is the case we want to go through all files and calculate the angels
        files = os.listdir(path)

        for file in files:
            print(file, ":")
            filepath = os.path.join(path, file)
            evaluate_file(filepath)
            print(" ")
            print("  ")
    else:
        evaluate_file(path)


if __name__ == "__main__":
    main()
