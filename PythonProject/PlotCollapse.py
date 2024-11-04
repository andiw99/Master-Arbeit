from Suite import crit_temp_measurement
import matplotlib.pyplot as plt
from FunctionsAndClasses import *
def main():
    equil_cutoff = 0.5
    process_file_func = process_new_mag_file_to_U_L
    simulation_path = "../../Generated content/Paper content/Binder intersection/m scalar/detailed/Tc/"

    selected_sizes = [64, 96, 128]
    selected_sizes = None
    exclude = None
    exclude = [32]
    selected_temps = None
    #selected_temps = np.linspace(1.91, 2.0, 10)
    min_T_plot = None
    max_T_plot = None
    min_T_plot = 1.92
    max_T_plot = 2.02
    highlight = False

    T_c = 1.969     # calculated from plotBinder skript

    results = crit_temp_measurement.construct_results(simulation_path, equil_cutoff, selected_temps=selected_temps,
                          selected_sizes=selected_sizes, value_name="U_L", file_ending="mag",
                          process_file_func=process_file_func)
    print(results)

    config = {"labelhorizontalalignment": "center",
              "labelverticalalignment": "bottom",
              "increasefontsize": 0.6,
              "labelrotation": 90,}

    fig, ax = plt.subplots(1, 1)

    for i,size in enumerate(sorted(results.keys())):
        x = (results[size]['T'] - T_c) / T_c * size # ** (1/nu) but we assume that nu is one?
        if highlight:
            if i == 0:
                color = "C1"
            else:
                color = colors[2*(i-1)]
        else:
            color = colors[5 * i]

        ax.plot(x, results[size]['U_L'], marker="s", **get_point_kwargs_color(color),  label=rf"$L_\parallel = {size}$")

    ax.set_ylabel(r"Binder cumulant $U_L$")
    ax.set_xlabel(r"rescaled reduced temperature $\varepsilon \, L^{1/\nu}$")

    configure_ax(fig, ax, config)
    create_directory_if_not_exists(f"{simulation_path}/plots/")
    fig.savefig(simulation_path + f"/plots/collapse.png", format="png",
                dpi=300, transparent=False)
    fig.savefig(simulation_path + f"/plots/collapse-svg.svg", format="svg")
    plt.show()

if __name__ == "__main__":
    main()