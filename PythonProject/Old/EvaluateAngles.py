import math
import numpy as np

def calculate_angle_with_y_axis(point1, point2):
    # Calculate the angle with the y-axis using arctan
    dz = np.abs(point2[2] - point1[2])
    dy = np.abs(point2[1] - point1[1])
    angle_rad = math.asin(dz / dy)
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
def average_angles_over_lines(file_path, x_tolerance=0.03, y_tolerance=3, nr_dimers=12):
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
    for si1, si2 in grouped_coordinates:
            # Calculate angle for each pair of consecutive points in the group
        angle = calculate_angle_with_y_axis(si1, si2)
        angles.append(angle)

    # Average angles over all lines
    if angles:
        average_angle = sum(angles) / len(angles)
        return average_angle
    else:
        return None

# Example usage:
file_path = "../../../Generated content/p1-si-relaxed-geom-pbesol-pbe-pz/pbesol/sisurc4x2-12dim-pbesol.xyz"  # Replace with the path to your file
average_angle = average_angles_over_lines(file_path)
if average_angle is not None:
    print("Average angle with y-axis:", average_angle)
else:
    print("No data points found to calculate the angle.")
