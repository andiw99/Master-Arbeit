from FunctionsAndClasses import *

path = "../../../Generated content/Silicon/Subsystems/Suite/Exp/h=2800/Jx_Jy=60/2800/Tc/80/20021.358570"

threshold = 0.1
file_ending = "mag"
value="m"


print(process_temp_folder(path, threshold, file_ending, value, process_file_func=recalculate_mag_file_to_U_L))