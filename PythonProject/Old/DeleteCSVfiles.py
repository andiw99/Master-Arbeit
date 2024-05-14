import os
import glob

def delete_csv_files(directory):
    for root, dirs, files in os.walk(directory):
        for dir in dirs:
            dir_path = os.path.join(root, dir)
            csv_files = glob.glob(os.path.join(dir_path, '*.csv'))
            for csv_file in csv_files[1:]:
                os.remove(csv_file)

# Call the function with the directory path
delete_csv_files('../../../Generated content/Final/Quenches-old')
