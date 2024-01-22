import os
import time
from datetime import datetime

def delete_old_files(folder_path, days_threshold, hours_threshold):
    current_time = time.time()
    threshold_time = current_time - (days_threshold * 24 * 60 * 60 + hours_threshold * 60 * 60)  # Convert days to seconds
    thres_time = datetime.fromtimestamp(threshold_time).strftime(
        '%Y-%m-%d %H:%M:%S')
    print("threshold_time:", thres_time)
    for root, dirs, files in os.walk(folder_path):
        for file_name in files:
            #print(file_name)
            file_path = os.path.join(root, file_name)
            file_creation_time = os.path.getctime(file_path)
            creation_time = datetime.fromtimestamp(file_creation_time).strftime(
                '%Y-%m-%d %H:%M:%S')
            print(file_name, "created", creation_time)
            if file_creation_time < threshold_time:
                try:
                    os.remove(file_path)
                    print(f"Deleted file: {file_path}")
                except Exception as e:
                    print(f"Error deleting file {file_path}: {e}")

# Example usage
def main():
    folder_to_clean = "../../../Generated content/Silicon/Subsystems/h/Largest h"
    days_threshold = 4  # Files older than 7 days will be deleted
    hours_threshold = 10
    delete_old_files(folder_to_clean, days_threshold, hours_threshold)

if __name__ == "__main__":
    main()