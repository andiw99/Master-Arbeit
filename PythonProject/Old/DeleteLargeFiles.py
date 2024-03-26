import os

def delete_csv_files(directory_path, threshold_size_mb):
    for root, dirs, files in os.walk(directory_path):
        # Check if there are multiple CSV files in the current directory
        csv_files = [file for file in files if file.lower().endswith('.csv')]
        if len(csv_files) > 1:
            # Sort CSV files by size (largest first)
            csv_files.sort(key=lambda x: os.path.getsize(os.path.join(root, x)), reverse=True)

            # Delete CSV files exceeding the threshold size
            total_size_mb = 0
            for csv_file in csv_files[1:]:  # Skip the first (largest) file
                file_path = os.path.join(root, csv_file)
                file_size_mb = os.path.getsize(file_path) / (1024 * 1024)  # Convert to MB

                if file_size_mb > threshold_size_mb:
                    print(f"Deleting {file_path} (Size: {file_size_mb:.2f} MB)")
                    os.remove(file_path)
                else:
                    total_size_mb += file_size_mb

            print(f"Keeping {csv_files[0]} (Total Size of Kept Files: {total_size_mb:.2f} MB)")

def delete_mag_files(directory_path, threshold_size_mb):
    for root, dirs, files in os.walk(directory_path):
        # Check if there are multiple CSV files in the current directory
        csv_files = [file for file in files if file.lower().endswith('.mag')]
        if len(csv_files) > 1:
            # Sort CSV files by size (largest first)
            csv_files.sort(key=lambda x: os.path.getsize(os.path.join(root, x)), reverse=True)

            # Delete CSV files exceeding the threshold size
            total_size_mb = 0
            for csv_file in csv_files[1:]:  # Skip the first (largest) file
                file_path = os.path.join(root, csv_file)
                file_size_mb = os.path.getsize(file_path) / (1024 * 1024)  # Convert to MB

                if file_size_mb > threshold_size_mb:
                    print(f"Deleting {file_path} (Size: {file_size_mb:.2f} MB)")
                    os.remove(file_path)
                else:
                    total_size_mb += file_size_mb

            print(f"Keeping {csv_files[0]} (Total Size of Kept Files: {total_size_mb:.2f} MB)")

def main():
    directory_path = "../../../Generated content/Final/z-measurement-small/h=2"
    threshold_size_mb = 100

    if not os.path.isdir(directory_path):
        print(f"Error: '{directory_path}' is not a valid directory.")
        return

    delete_csv_files(directory_path, threshold_size_mb)
    delete_mag_files(directory_path, threshold_size_mb)

if __name__ == "__main__":
    main()
