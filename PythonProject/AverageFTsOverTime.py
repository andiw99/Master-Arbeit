import numpy as np

from FunctionsAndClasses import *


def plot_process(size_dic, t_eq, quench=True, quench_zoom=1, max_nr_curves=np.infty, y_scale="log"):
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
    simulation_path = "../../Generated content/Silicon/Quench/h/High h"
    cut_zero_impuls = True
    quench = True
    scale_time = True
    max_nr_curves = 5
    quench_zoom = 20
    y_scale = "log"
    y_lower_lim = 0.05
    min_tau_fit = 2
    max_tau_fit = np.infty

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
                            ft_k, ft_l = average_ft(settingpath)
                            for t in ft_k:
                                p_k = get_frequencies_fftw_order(len(ft_k[t]))
                                if cut_zero_impuls:
                                    p_k, ft_k[t] = cut_zero_imp(p_k, ft_k[t])
                                popt_x, perr_x = fit_lorentz(p_k, ft_k[t])
                                xix = np.minimum(np.abs(popt_x[0]), Lx)
                                t_xix[setting][t] = xix
                            for t in ft_l:
                                p_l = get_frequencies_fftw_order(len(ft_l[t]))
                                if cut_zero_impuls:
                                    p_l, ft_l[t] = cut_zero_imp(p_l, ft_l[t])
                                popt_y, perr_y = fit_lorentz(p_l, ft_l[t])
                                xiy = np.minimum(np.abs(popt_y[0]), Ly)
                                t_xiy[setting][t] = xiy
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


    fig, ax = plot_process(size_x_dic, t_eq, quench=quench, quench_zoom=quench_zoom, max_nr_curves=max_nr_curves, y_scale="log")
    ax.set_title("Time during quench is scaled")
    ax.set_ylabel(r"$\xi_x$")
    ax.set_xlabel(
        r"$t \qquad \qquad \qquad \qquad \qquad \qquad t/\tau \qquad \qquad  \qquad \qquad \qquad \qquad t $")
    fig.savefig(simulation_path + f"/xix-process.png", format="png", dpi=300)
    plt.show()

    fig, ax = plot_process(size_y_dic, t_eq, quench=quench, quench_zoom=quench_zoom, max_nr_curves=max_nr_curves, y_scale="log")
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
            popt_x, _ = curve_fit(linear_fit, np.log(tau_fit), np.log(xix_fit))
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