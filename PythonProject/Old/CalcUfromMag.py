from FunctionsAndClasses import *

path = "../../../Generated content/Silicon/Subsystems/Suite/Exp/h=3300/Jx_Jy=31-2/3300/Tc/48/30843.936465"

threshold=0.1
file_ending = "mag"
value="m"


print(process_temp_folder(path, threshold, file_ending, value, process_file_func=process_new_mag_file_to_U_L))