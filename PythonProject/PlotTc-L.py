from Suite import crit_temp_measurement
import matplotlib.pyplot as plt
from FunctionsAndClasses import *
def main():
    equil_cutoff = 0.1
    process_file_func = recalculate_mag_file_to_U_L
    simulation_path = "../../Generated content/Final/Nu-L-old-selected"

    selected_sizes = [128, 16, 48, 72]
    J_para = 3.11
    selected_sizes = None
    results = crit_temp_measurement.construct_results(simulation_path, equil_cutoff, selected_temps=None,
                          selected_sizes=selected_sizes, value_name="U_L", file_ending="mag",
                          process_file_func=process_file_func)

    config = {"labelhorizontalalignment": "right",
              "increasefontsize": 0.3}

    fig, ax = crit_temp_measurement.plot_Tc_L_dependence(results, J_para=J_para)
    plt.show()

if __name__ == "__main__":
    main()