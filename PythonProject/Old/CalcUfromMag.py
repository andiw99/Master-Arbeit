from FunctionsAndClasses import *

#path = "../../../Generated content/Final/CriticalTemperature/Tc/128/43607.941684"
path = "../../../Generated content/Silicon/Subsystems/Suite/Exp/h=40000/Jx_Jy=60-2/40000/Tc/128/20021.358570"
path = "../../../Generated content/Silicon/Subsystems/Suite/Exp/h=930/Jx_Jy=31/930/Tc/128/18506.361879"

threshold = 0.5
file_ending = "mag"
value="m"


print(process_temp_folder(path, threshold, file_ending, value, process_file_func=recalculate_mag_file_to_U_L))