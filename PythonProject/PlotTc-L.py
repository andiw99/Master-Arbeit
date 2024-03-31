from Suite import crit_temp_measurement
import matplotlib.pyplot as plt
from FunctionsAndClasses import *
def main():
    equil_cutoff = 0.1
    equil_cutoffs = np.linspace(0.02, 0.3, 1)
    for equil_cutoff in equil_cutoffs:
        print("equil_cutoff = ", equil_cutoff)
        process_file_func = recalculate_mag_file_to_U_L
        simulation_path = "../../Generated content/Final/Nu-L-old (copy)"

        J_para = 3.11
        selected_sizes = [128, 16, 48, 72]
        selected_sizes = None
        results = crit_temp_measurement.construct_results(simulation_path, equil_cutoff, selected_temps=None,
                              selected_sizes=selected_sizes, value_name="U_L", file_ending="mag",
                              process_file_func=process_file_func)

        config = {"labelhorizontalalignment": "right",
                  "increasefontsize": 0.6}
        fig, ax = crit_temp_measurement.plot_value_curve(simulation_path, results, crit_point=None, value_name="U_L", title="Binder Cumulant on T",
                                               plotname="cum_time_avg", equil_error=None, config=config)
        plt.show()

        fig, ax = crit_temp_measurement.fit_and_plot_nu(simulation_path, results, crit_point=None, value_name="U_L")
        create_directory_if_not_exists(f"{simulation_path}/plots/")
        fig.savefig(simulation_path + f"/plots/nu-{equil_cutoff}.png", format="png",
                    dpi=300, transparent=False)
        fig.savefig(simulation_path + f"/plots/nu-{equil_cutoff}-svg.svg", format="svg")
        plt.show()

    fig, ax = crit_temp_measurement.plot_Tc_L_dependence(results, J_para=J_para)
    create_directory_if_not_exists(simulation_path + "/plots/")
    fig.savefig(simulation_path + "/plots/Tc_L.png", dpi=400)
    fig.savefig(simulation_path + "/plots/Tc_L-200-dpi.png", dpi=200)
    fig.savefig(simulation_path + "/plots/Tc_L.svg", format="svg")
    plt.show()


if __name__ == "__main__":
    main()