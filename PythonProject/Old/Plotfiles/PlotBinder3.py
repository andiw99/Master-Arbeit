from Suite import crit_temp_measurement
import matplotlib.pyplot as plt
from FunctionsAndClasses import *
def main():
    equil_cutoff = 0.5
    process_file_func = recalculate_mag_file_to_U_L
    simulation_path = "../../Generated content/Final/CriticalTemperature/Tc-copy/Tc"
    simulation_path = "../../Generated content/Final/CriticalTemperature/Nu-L-old/"
    selected_sizes = None
    selected_sizes = [128, 16, 48, 72]
    results = crit_temp_measurement.construct_results(simulation_path, equil_cutoff, selected_temps=None,
                          selected_sizes=selected_sizes, value_name="U_L", file_ending="mag",
                          process_file_func=process_file_func)
    print(results)
    try:
        T_cs, val_intersections = get_first_intersections(results, "U_L")
        T_c = T_cs[0]
        val_intersection = val_intersections[0]
    except:
        pass
    crit_point= None

    config = {"labelhorizontalalignment": "right",
              "increasefontsize": 0.6}

    fig, ax = crit_temp_measurement.plot_value_curve(simulation_path, results, crit_point=crit_point, value_name="U_L", title="Binder Cumulant on T",
                                           plotname="cum_time_avg", equil_error=None, config=config)
    plt.show()
    #T_range = (0.836012, 0.870134)
    # T_range=None
    # fig, ax = crit_temp_measurement.fit_and_plot_nu(simulation_path, results, crit_point=None, T_range=T_range, value_name="U_L")
    # create_directory_if_not_exists(f"{simulation_path}/plots/")
    # fig.savefig(simulation_path + f"/plots/nu.png", format="png",
    #             dpi=300, transparent=False)
    # plt.show()


    selected_temps = np.linspace(26533.333333, 28111.111111, 10)
    T_min = 29071.961123
    T_max = 31494.624550

    selected_temps = np.array(next(iter(results.values()))["T"])
    selected_temps = selected_temps[(selected_temps >= T_min)  & (selected_temps <= T_max)]

    print(selected_temps)
    results = crit_temp_measurement.construct_results(simulation_path, equil_cutoff, selected_temps=selected_temps,
                          selected_sizes=None, value_name="U_L", file_ending="mag",
                          process_file_func=process_file_func)

    T_cs, val_intersections = get_first_intersections(results, "U_L")
    T_c = T_cs[0]
    val_intersection = val_intersections[0]
    crit_point= (T_c, val_intersection)

    config = {"labelhorizontalalignment": "right",
              "increasefontsize": 0.6}

    crit_temp_measurement.plot_value_curve(simulation_path, results, crit_point=crit_point, value_name="U_L", title="Binder Cumulant on T",
                                           plotname="cum_time_avg", equil_error="zoom", config=config)

    plt.show()

def binder_plot_static_sec():
    equil_cutoff = 0.5
    process_file_func = recalculate_mag_file_to_U_L
    simulation_path = "../../Generated content/Final/CriticalTemperature/Tc-copy/Tc"

    selected_sizes = [128, 16, 48, 72]
    selected_sizes = None
    results = crit_temp_measurement.construct_results(simulation_path, equil_cutoff, selected_temps=None,
                          selected_sizes=selected_sizes, value_name="U_L", file_ending="mag",
                          process_file_func=process_file_func)
    print(results)
    try:
        T_cs, val_intersections = get_first_intersections(results, "U_L")
        T_c = T_cs[0]
        val_intersection = val_intersections[0]
    except:
        pass
    crit_point= None

    config = {"labelhorizontalalignment": "right",
              "increasefontsize": 0.3}

    fig, ax = crit_temp_measurement.plot_value_curve(simulation_path, results, crit_point=crit_point, value_name="U_L", title="Binder Cumulant on T",
                                           plotname="cum_time_avg", equil_error=None, config=config)
    plt.show()


    selected_temps = np.linspace(26533.333333, 28111.111111, 10)
    T_min = 29071.961123
    T_max = 31494.624550

    selected_temps = np.array(next(iter(results.values()))["T"])
    selected_temps = selected_temps[(selected_temps >= T_min)  & (selected_temps <= T_max)]

    print(selected_temps)
    results = crit_temp_measurement.construct_results(simulation_path, equil_cutoff, selected_temps=selected_temps,
                          selected_sizes=None, value_name="U_L", file_ending="mag",
                          process_file_func=process_file_func)

    T_cs, val_intersections = get_first_intersections(results, "U_L")
    T_c = T_cs[0]
    val_intersection = val_intersections[0]
    crit_point= (T_c, val_intersection)

    config = {"labelhorizontalalignment": "right",
              "increasefontsize": 0.3}

    crit_temp_measurement.plot_value_curve(simulation_path, results, crit_point=crit_point, value_name="U_L", title="Binder Cumulant on T",
                                           plotname="cum_time_avg", equil_error="zoom", config=config)

    plt.show()

if __name__ == "__main__":
    main()