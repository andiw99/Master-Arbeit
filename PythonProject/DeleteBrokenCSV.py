import os
import csv
from FunctionsAndClasses import *

def has_rows_of_different_lengths(file_path):
    with open(file_path, 'r') as csvfile:
        reader = csv.reader(csvfile)
        lengths = set([len(row) for row in reader])
        print(lengths)
        return len(lengths) > 1

def delete_files_with_rows_of_different_lengths(directory):
    nr_files = 0
    deleted = 0
    for root, dirs, files in os.walk(directory):
        for filename in files:
            if filename.endswith(".csv"):
                file_path = os.path.join(root, filename)
                nr_files += 1
                if has_rows_of_different_lengths(file_path):
                    os.remove(file_path)
                    deleted += 1
                    print(f"Deleted: {file_path}")
    print(f"Checked {nr_files} files, deleted {deleted}")

# Example: Provide the path to the directory
directory_path = "../../Generated content/Subsystems/Silicon AA Even 2"
delete_files_with_rows_of_different_lengths(directory_path)
