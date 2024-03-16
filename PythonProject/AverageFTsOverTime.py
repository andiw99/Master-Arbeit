import numpy as np

from FunctionsAndClasses import *
from decimal import Decimal


def plot_process(size_dic, t_eq, t_eq_end=0, quench=True, quench_zoom=1, max_nr_curves=np.infty, y_scale="log", direction="parallel", sim_path=None):
    if quench:
        setting_var = r"$\tau$"
    # first we need to find out which taus we want to plot
    # we use just the messreihe with the most taus
    most_tau_size = 0
    max_nr_taus = 0
    for size in size_dic:
        nr_taus = len(size_dic[size])
        if nr_taus > max_nr_taus:
            most_tau_size = size
            max_nr_taus = nr_taus

    t_xi_largest = size_dic[most_tau_size]
    max_nr_curves = min(max_nr_curves, len(t_xi_largest))
    valid_inds = np.linspace(0, len(t_xi_largest) - 1, max_nr_curves, dtype=int)
    taus = np.array(list(t_xi_largest.keys()))[np.argsort([float(key) for key in t_xi_largest.keys()])]
    taus_plot = taus[valid_inds]
    # now we can basically do what we did below?
    fig, axes = plt.subplots(1, 2, figsize=(10, 7), gridspec_kw={'width_ratios': [1, 3]})
    ax_equil = axes[0]
    ax_quench = axes[1]
    color_ind = 0
    color_dic = {}
    print(taus_plot)
    for marker_ind, size in enumerate(size_dic):
        t_xi = size_dic[size]
        for i, tau in enumerate(taus_plot):
            try:
                t = np.array(list(t_xi[tau].keys()))
                if quench:
                    if t_eq == 0:
                        folder_path = f"{sim_path}/{size}/{tau}"
                        csv_path = find_first_csv_file(folder_path)
                        t_eq = find_time_of_temperature_change(csv_path)
                    # instead of using this complicated function in this case we
                    # only need to separate the equilibration and the quenching and then divide by tau
                    # t, t_q_s = rescale_t(t, float(tau), t_eq, quench_zoom)
                    t_equil = t[t <= t_eq] - t_eq
                    end_quench_t = np.max(t) - t_eq_end
                    t_quench = (t[(t > t_eq) & (t < end_quench_t)] - t_eq) / float(tau)
                    # okay and we set the quench beginning to t = 0
                xi = np.array(list(t_xi[tau].values()))
                xi_equil = xi[t <= t_eq]
                xi_quench = xi[(t > t_eq) & (t < end_quench_t)]
                if i == 0:
                    pre_equil_xi = np.mean(xi_equil)
                if marker_ind == 0:
                    color = colors[color_ind]
                    color_ind += 1
                    color_dic[tau] = color  # we save the color for the next size
                    # for the legend
                    ax_quench.plot([], [], linestyle="", marker="s", color=color, label=f"{setting_var} = {float(tau):.2f}")
                else:
                    try:
                        color = color_dic[tau]
                    except KeyError as e:
                        # if the setting is non existent, we use the next color in colors
                        color = colors[color_ind]
                        color_ind += 1
                ax_equil.plot(t_equil, xi_equil, linestyle="", markersize=5, marker=markers[marker_ind],
                               color=color, markerfacecolor="none")  # label=rf"$\xi_x$  {setting_var} = {float(setting):.2f}",)
                ax_quench.plot(t_quench, xi_quench, linestyle="", markersize=5, marker=markers[marker_ind],
                               color=color, markerfacecolor="none")  # label=rf"$\xi_x$  {setting_var} = {float(setting):.2f}",)
            except KeyError:
                print(f"tau = {tau} not available in size = {size}")
        ax_quench.plot([], [], linestyle="", marker=markers[marker_ind], label=f"L = {size}", color="black",  markerfacecolor="none")
    ax_quench.set_yscale(y_scale)
    ax_equil.set_yscale(y_scale)
    if y_scale == "log":
        ax_quench.set_ylim(0.8 * pre_equil_xi, ax_quench.get_ylim()[1])
        ax_equil.set_ylim(ax_quench.get_ylim())     # same limits for the two axes
    ax_equil.set_xlim(ax_equil.get_xlim()[0], 0)
    ax_quench.set_xlim(0, ax_quench.get_xlim()[1])
    quench_config = {
        "titlesize": 0,
        "ytickfontsize": 0,
        "ylabelsize": 0,
        "y_tickwidth": 0,
        "y_ticklength": 0
    }
    equil_config = {
        "nr_x_major_ticks": 2
    }
    configure_ax(fig, ax_quench, quench_config)
    configure_ax(fig, ax_equil, equil_config)
    fig.subplots_adjust(wspace=0.01)
    return fig, axes

def plot_process2(size_dic, t_eq, quench=True, quench_zoom=1, max_nr_curves=np.infty, y_scale="log"):
    if quench:
        setting_var = r"$\tau$"
    # first we need to find out which taus we want to plot
    # we use just the messreihe with the most taus
    most_tau_size = 0
    max_nr_taus = 0
    for size in size_dic:
        nr_taus = len(size_dic[size])
        if nr_taus > max_nr_taus:
            most_tau_size = size
            max_nr_taus = nr_taus

    t_xi_largest = size_dic[most_tau_size]
    max_nr_curves = min(max_nr_curves, len(t_xi_largest))
    valid_inds = np.linspace(0, len(t_xi_largest) - 1, max_nr_curves, dtype=int)
    taus = np.array(list(t_xi_largest.keys()))[np.argsort([float(key) for key in t_xi_largest.keys()])]
    taus_plot = taus[valid_inds]
    # now we can basically do what we did below?
    # Okay so for this the plan is to have only the equilibration shown and then
    # the quench process, but in two separat plots
    fig, ax = plt.subplots(1, 1)
    color_ind = 0
    color_dic = {}
    print(taus_plot)
    for marker_ind, size in enumerate(size_dic):
        t_xi = size_dic[size]
        for i, tau in enumerate(taus_plot):
            try:
                t = list(t_xi[tau].keys())
                if quench:
                    t, t_q_s = rescale_t(t, float(tau), t_eq, quench_zoom)
                xi = np.array(list(t_xi[tau].values()))
                if i == 0:
                    pre_equil_xi = np.mean(xi[t < 0])
                if marker_ind == 0:
                    color = colors[color_ind]
                    color_ind += 1
                    color_dic[tau] = color  # we save the color for the next size
                    # for the legend
                    ax.plot([], [], linestyle="", marker="s", color=color, label=f"{setting_var} = {float(tau):.2f}")
                else:
                    try:
                        color = color_dic[tau]
                    except KeyError as e:
                        # if the setting is non existent, we use the next color in colors
                        color = colors[color_ind]
                        color_ind += 1
                ax.plot(t, xi, linestyle="", markersize=4, marker=markers[marker_ind], color=color)  # label=rf"$\xi_x$  {setting_var} = {float(setting):.2f}",)
            except KeyError:
                print(f"tau = {tau} not available in size = {size}")
        ax.plot([], [], linestyle="", marker=markers[marker_ind], label=f"L = {size}", color="black")
    ax.set_yscale(y_scale)
    if y_scale == "log":
        ax.set_ylim(0.8 * pre_equil_xi, ax.get_ylim()[1])
    ax.vlines((0, np.max(t) - t_eq), ax.get_ylim()[0], ax.get_ylim()[1], linestyles="dashed", color="grey",
              alpha=0.5)
    configure_ax(fig, ax)
    return fig, ax
def main():
    #simulation_path = "../../Generated content/Silicon/Quench/Meshs/New/StructFactTest/PBC XY Switch/"
    simulation_path = "../../Generated content/Silicon/Subsystems/Suite/h/1.7320508075688776/Quench/256"
    cut_zero_impuls = True
    quench = True
    scale_time = True
    max_nr_curves = 5
    quench_zoom = 50
    y_scale = "log"
    y_lower_lim = 0.1
    set_fts_to_zero = False
    min_tau_fit = 100
    max_tau_fit = np.infty
    plot_struct = True
    cut_around_peak = True
    peak_cut_threshold = 0.1
    min_points_fraction = 0.8
    fitfunc = lorentz_offset
    struct_func_time = 00000
    t_xix = {}
    t_xiy = {}
    size_x_dic = {}
    size_y_dic = {}

    for size in os.listdir(simulation_path):
        t_xix = {}
        t_xiy = {}
        if (size != "plots") & (size[0] != "."):
            sizepath = os.path.join(simulation_path, size)
            if os.path.isdir(sizepath):
                for setting in os.listdir(sizepath):
                    if (setting != "plots") & (setting[0] != "."):
                        settingpath = os.path.join(sizepath, setting)
                        print(settingpath)
                        parapath = find_first_txt_file(simulation_path)
                        parameters = read_parameters_txt(parapath)

                        Lx = parameters["subsystem_Lx"]
                        Ly = parameters["subsystem_Ly"]
                        if os.path.isdir(settingpath):
                            t_xix[setting] = {}
                            t_xiy[setting] = {}
                            ft_k, ft_l = average_ft_unequal_times(settingpath)
                            for t in ft_k:
                                ft_k_fit = ft_k[t]
                                ft_k_fit, p_k = prepare_fit_data(cut_around_peak, cut_zero_impuls, ft_k_fit,
                                                                 peak_cut_threshold, set_fts_to_zero, min_points_fraction)

                                popt_x, perr_x = fit_lorentz(p_k, ft_k_fit,
                                                             fitfunc=fitfunc)
                                # print("offset = ", popt_x[1])
                                xix = np.minimum(np.abs(popt_x[0]), Lx)
                                t_xix[setting][t] = xix
                            for t in ft_l:
                                ft_l_fit = ft_l[t]
                                ft_l_fit, p_l = prepare_fit_data(cut_around_peak, cut_zero_impuls, ft_l_fit,
                                                                 peak_cut_threshold, set_fts_to_zero, min_points_fraction)

                                popt_y, perr_y = fit_lorentz(p_l, ft_l_fit,
                                                                fitfunc=fitfunc)
                                xiy = np.minimum(np.abs(popt_y[0]), Ly)
                                t_xiy[setting][t] = xiy

                            if plot_struct:
                                if struct_func_time:
                                    t = find_closest_key(ft_k, struct_func_time)
                                else:
                                    t = list(ft_k.keys())[len(ft_k.keys())//2]
                                ft_k_fit = ft_k[t]
                                ft_l_fit = ft_l[t]
                                p_k_plot, ft_k_plot = cut_zero_imp(get_frequencies_fftw_order(len(ft_k_fit)), ft_k_fit)
                                p_l_plot, ft_l_plot = cut_zero_imp(get_frequencies_fftw_order(len(ft_l_fit)), ft_l_fit)
                                ft_k_fit, p_k = prepare_fit_data(cut_around_peak, cut_zero_impuls, ft_k_fit,
                                                                 peak_cut_threshold, set_fts_to_zero, min_points_fraction)
                                ft_k_min = np.min(ft_k_fit)
                                ft_l_fit, p_l = prepare_fit_data(cut_around_peak, cut_zero_impuls, ft_l_fit,
                                                                 peak_cut_threshold, set_fts_to_zero, min_points_fraction)
                                ft_l_min = np.min(ft_l_fit)
                                fig, axes = plot_struct_func(p_k_plot, p_l_plot,
                                                             ft_k_plot,
                                                             ft_l_plot)
                                popt_x, perr_x = fit_lorentz(p_k, ft_k_fit,
                                                             fitfunc=fitfunc)
                                popt_y, perr_y = fit_lorentz(p_l, ft_l_fit,
                                                             fitfunc=fitfunc)
                                p = np.linspace(-np.pi, np.pi, 200)
                                # Since for the fit we use the cut_around_peak function which
                                lorentz_x = fitfunc(p, *popt_x)
                                lorentz_y = fitfunc(p, *popt_y)
                                axes[0].plot(p, lorentz_x,
                                             label="Lorentzian fit")
                                axes[1].plot(p, lorentz_y,
                                             label="Lorentzian fit")
                                axes[1].hlines(ft_l_min, -np.pi, np.pi, alpha=0.3, linestyles="--", label=f"fit limit = {ft_l_min:.2E} \n "
                                                                                                          f"min points = {min_points_fraction}")
                                axes[0].hlines(ft_k_min, -np.pi, np.pi, alpha=0.3, linestyles="--", label=f"fit limit = {ft_k_min:.2E}\n"
                                                                                                          f"min points = {min_points_fraction}")
                                axes[0].set_title(
                                    rf"$\xi_x = {xix:.2f}  \quad t = {t} \quad  \tau = {setting}$")
                                axes[1].set_title(
                                    rf"$\xi_y = {xiy:.2f}  \quad t = {t} \quad  \tau = {setting}$")
                                config = {"legendlocation" : "upper right"}
                                configure_ax(fig, axes[0], config)
                                configure_ax(fig, axes[1], config)
                                #axes[0].legend(loc="upper right")
                                #axes[1].legend(loc="upper right")
                                #plt.tight_layout()
                                create_directory_if_not_exists(simulation_path + "/plots")
                                plt.savefig(simulation_path + f"/plots/{setting}-{t}-{len(ft_k[t])}", format="png")
                                plt.show()
                    print("size: ", int(size))
                    print(t_xix.keys())
                    size_x_dic[int(size)] = t_xix.copy()
                    size_y_dic[int(size)] = t_xiy.copy()
    # print(t_xiy)
    # print("\n\n")
    # print(t_xiy.keys())


    t_eq = 0
    t_q_s = 0
    if scale_time:
        parapath = find_first_txt_file(simulation_path)
        parameters = read_parameters_txt(parapath)
        t_eq = parameters["equil_time"]

    setting_var = "T"
    if quench:
        setting_var = r"$\tau$"


    fig, ax = plot_process(size_x_dic, t_eq, t_eq_end=t_eq, quench=quench, quench_zoom=quench_zoom,
                           max_nr_curves=max_nr_curves, y_scale="log", sim_path=simulation_path)
    ax_equil = ax[0]
    ax_quench = ax[1]
    ax_equil.set_title("Equilibration")
    ax_quench.set_title("Quench")
    fig.suptitle("Time during quench is scaled")
    ax_equil.set_ylabel(r"$\xi_\parallel / a_\parallel$")
    ax_equil.set_xlabel("t/ns")
    ax_quench.set_xlabel(
        r"$t / \tau_Q$")
    fig.savefig(simulation_path + f"/xix-process.png", format="png", dpi=300)
    plt.show()
    exit()

    fig, ax = plot_process(size_y_dic, t_eq, quench=quench, quench_zoom=quench_zoom, max_nr_curves=max_nr_curves,
                           y_scale="log", sim_path=simulation_path)
    ax.set_title("Time during quench is scaled")
    ax.set_ylabel(r"$\xi_y$")
    ax.set_xlabel(
        r"$t \qquad \qquad \qquad \qquad \qquad \qquad t/\tau \qquad \qquad  \qquad \qquad \qquad \qquad t $")
    fig.savefig(simulation_path + f"/xiy-process.png", format="png", dpi=300)
    plt.show()

    if quench:
        figx, axx = plt.subplots(1, 1)
        figy, axy = plt.subplots(1, 1)
        axx.set_yscale("log")
        axy.set_yscale("log")
        axx.set_xscale("log")
        axy.set_xscale("log")
        xix_after_quench_dic = {}
        xiy_after_quench_dic = {}
        for marker_ind, size in enumerate(size_y_dic):
            tau_arr = []
            # okay i want to show the scaling so I extract the xi after t_quench
            # do i still know t_quench?
            t_xix = size_x_dic[size]
            t_xiy = size_y_dic[size]

            xix_after_quench_dic[size] = []
            xiy_after_quench_dic[size] = []

            keys = np.array(list(t_xix.keys()))[np.argsort([float(key) for key in t_xix.keys()])]
            for i, tau in enumerate(keys):
                t = np.array(list(t_xix[tau].keys()))
                t_after_quench = np.max(t) - t_eq
                xix_arr = np.array(list(t_xix[tau].values()))      # this is all the xi after time
                xiy_arr = np.array(list(t_xiy[tau].values()))
                xix = xix_arr[t <= t_after_quench][-1]
                xiy = xiy_arr[t <= t_after_quench][-1]
                xix_after_quench_dic[size].append(xix)
                xiy_after_quench_dic[size].append(xiy)
                tau_arr.append(float(tau))
            tau_arr = np.array(tau_arr)
            xix_fit = np.array(xix_after_quench_dic[size])[(tau_arr > min_tau_fit) & (tau_arr < max_tau_fit)]
            xiy_fit = np.array(xiy_after_quench_dic[size])[(tau_arr > min_tau_fit) & (tau_arr < max_tau_fit)]
            tau_fit = tau_arr[(tau_arr > min_tau_fit) & (tau_arr < max_tau_fit)]
            # do the fitting
            print("log tau")
            print(np.log(tau_fit))
            popt_x, _ = curve_fit(linear_fit, np.log(tau_fit), np.log(xix_fit))
            print("xi:")
            print(xiy_fit)
            popt_y, _ = curve_fit(linear_fit, np.log(tau_fit), np.log(xiy_fit))

            axx.plot(tau_arr, xix_after_quench_dic[size], linestyle="",     color="C" + str(marker_ind), marker=markers[marker_ind], label=f"$L_x = ${size}")
            axy.plot(tau_arr, xiy_after_quench_dic[size], linestyle="",     color="C" + str(marker_ind), marker=markers[marker_ind], label=f"$L_x = ${size}")
            axx.plot(tau_fit, poly(tau_fit, popt_x[0], np.exp(popt_x[1])),  color="C" + str(marker_ind), label=r"$\frac{\nu}{1 + \nu z} =$" + f"{popt_x[0]:.2f}" )
            axy.plot(tau_fit, poly(tau_fit, popt_y[0], np.exp(popt_y[1])),  color="C" + str(marker_ind), label=r"$\frac{\nu}{1 + \nu z} =$" + f"{popt_y[0]:.2f}" )

    axx.set_xlabel(r"$\tau$")
    axx.set_ylabel(r"$\xi_x$")
    axy.set_xlabel(r"$\tau$")
    axy.set_ylabel(r"$\xi_y$")

    configure_ax(figx, axx)
    figx.savefig(simulation_path + "/xix-quench-scaling.png", format="png", dpi=300)

    configure_ax(figy, axy)
    figy.savefig(simulation_path + "/xiy-quench-scaling.png", format="png", dpi=300)
    plt.show()




if __name__ == "__main__":
    main()