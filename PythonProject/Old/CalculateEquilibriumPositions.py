import csv
import numpy as np

def average_values_in_last_line(csv_file_path):
    with open(csv_file_path, 'r') as file:
        # Read the CSV file
        csv_reader = csv.reader(file)

        # Get the last line
        last_line = None
        for row in csv_reader:
            last_line = np.array(row)[2:]

        # Filter numeric values from the last line
        numeric_values = [float(value) for value in last_line if value != ""]

        # Separate positive and negative values
        positive_values = [value for value in numeric_values if value > 0]
        negative_values = [value for value in numeric_values if value < 0]

        # Calculate averages
        avg_positive = sum(positive_values) / len(positive_values) if positive_values else None
        avg_negative = sum(negative_values) / len(negative_values) if negative_values else None

        return avg_positive, avg_negative


# Example usage:
csv_file_path = '../../../Generated content/Silicon/Quench/EquilibriumPosition/p=2.15/400/4096.000000/0-19:0-local-andi-B550.csv'  # Replace with the path to your CSV file
result = average_values_in_last_line(csv_file_path)

print("Average of positive values:", result[0])
print("Average of negative values:", result[1])
