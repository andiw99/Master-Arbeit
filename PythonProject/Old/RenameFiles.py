import os

def rename_files(directory):
    for foldername, subfolders, filenames in os.walk(directory):
        for filename in filenames:
            print(filename)
            if len(filename) > 35:
                old_file_path = os.path.join(foldername, filename)

                # Extract file extension
                root, extension = os.path.splitext(filename)

                # Truncate file name to the first 7 and last 21 characters
                truncated_name = root[:7] + root[-22:]

                # Construct the new file name and path
                new_filename = truncated_name + extension
                new_file_path = os.path.join(foldername, new_filename)

                # Rename the file
                os.rename(old_file_path, new_file_path)
                print(f'Renamed: {old_file_path} -> {new_file_path}')

# Replace 'your_directory_path' with the path to the directory you want to process
def main():
    root_directory = '../../../Generated content/Silicon/Subsystems/OBC3'
    rename_files(root_directory)


if __name__ == '__main__':
    main()