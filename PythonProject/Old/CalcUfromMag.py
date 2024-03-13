from FunctionsAndClasses import *

path = "../../../Generated content/Silicon/Subsystems/Suite/h/Large Jx/Jx=10-Lx_Ly=32/0.4161791450287818/Tc/64/1.750317"

threshold = 0.5
file_ending = "mag"
value="m"


print(process_temp_folder(path, threshold, file_ending, value, process_file_func=recalculate_mag_file_to_U_L))