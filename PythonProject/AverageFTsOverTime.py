import numpy as np

from FunctionsAndClasses import *

def main():
    simulation_path = "../../Generated content/Silicon/Quench/OBC2"
    cut_zero_impuls = True
    quench = True
    scale_time = True
    t_xix = {}
    t_xiy = {}
    max_nr_curves = 5
    quench_zoom = 20
    y_scale = "log"
    y_lower_lim = 0.05

    size_x_dic = {}
    size_y_dic = {}

    for size in os.listdir(simulation_path):
        if (size != "plots") & (size[0] != "."):
            sizepath = os.path.join(simulation_path, size)
            if os.path.isdir(sizepath):
                for setting in os.listdir(sizepath):
                    if (setting != "plots") & (setting[0] != "."):
                        settingpath = os.path.join(sizepath, setting)

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

    fig, ax = plt.subplots(1, 1)
    ax.set_yscale(y_scale)
    for marker_ind, size in enumerate(size_x_dic):
        t_xi = size_x_dic[size]
        valid_inds = np.linspace(0, len(t_xi) - 1, len(t_xi), dtype=int)
        if max_nr_curves:
            valid_inds = np.linspace(0, len(t_xi) - 1, max_nr_curves, dtype=int)
        keys = np.array(list(t_xi.keys()))[np.argsort([float(key) for key in t_xi.keys()])]
        color_ind = 0
        for i, setting in enumerate(keys):
            if i in valid_inds:
                t = list(t_xi[setting].keys())
                if quench & scale_time:
                    tau = float(setting)
                    t, t_q_s = rescale_t(t, tau, t_eq, quench_zoom)
                xix = np.array(list(t_xi[setting].values()))
                if i == 0:
                    pre_equil_xi = np.mean(xix[t < 0])
                ax.plot(t, xix, linestyle="", markersize=4, marker=markers[marker_ind],  color=colors[color_ind]) # label=rf"$\xi_x$  {setting_var} = {float(setting):.2f}",)
                if marker_ind == 0:
                    ax.plot([], [], linestyle="", marker="s", color=colors[color_ind], label = f"{setting_var} = {float(setting):.2f}")
                color_ind += 1
        ax.plot([], [], linestyle="", marker=markers[marker_ind], label=f"L = {size}", color="black")
    configure_ax(fig, ax)
    if y_scale == "log":
        ax.set_ylim(0.8 * pre_equil_xi, ax.get_ylim()[1])
    ax.vlines((0, np.max(t) - t_eq), ax.get_ylim()[0], ax.get_ylim()[1], linestyles="dashed", color="grey", alpha=0.5)
    ax.set_title("Time during quench is scaled")
    ax.set_ylabel(r"$\xi_x$")
    ax.set_xlabel(r"$t \qquad \qquad \qquad \qquad \qquad \qquad t/\tau \qquad \qquad  \qquad \qquad \qquad \qquad t $")
    fig.savefig(simulation_path + "/xix-process.png", format="png", dpi=300)
    plt.show()

    fig, ax = plt.subplots(1, 1)
    ax.set_yscale(y_scale)
    for marker_ind, size in enumerate(size_y_dic):
        t_xi = size_y_dic[size]
        valid_inds = np.linspace(0, len(t_xi) - 1, len(t_xi), dtype=int)
        if max_nr_curves:
            valid_inds = np.linspace(0, len(t_xi) - 1, max_nr_curves, dtype=int)
        keys = np.array(list(t_xi.keys()))[np.argsort([float(key) for key in t_xi.keys()])]
        color_ind = 0
        for i, setting in enumerate(keys):
            if i in valid_inds:
                t = list(t_xi[setting].keys())
                if quench & scale_time:
                    tau = float(setting)
                    t, t_q_s = rescale_t(t, tau, t_eq, quench_zoom)
                xi = np.array(list(t_xi[setting].values()))
                if i == 0:
                    pre_equil_xi = np.mean(xi[t < 0])
                ax.plot(t, xi, linestyle="", markersize=4, marker=markers[marker_ind],  color=colors[color_ind]) # label=rf"$\xi_x$  {setting_var} = {float(setting):.2f}",)
                if marker_ind == 0:
                    ax.plot([], [], linestyle="", marker="s", color=colors[color_ind], label = f"{setting_var} = {float(setting):.2f}")
                color_ind += 1
        ax.plot([], [], linestyle="", marker=markers[marker_ind], label=f"L = {size}", color="black")
    configure_ax(fig, ax)
    if y_scale == "log":
        ax.set_ylim(pre_equil_xi * 0.8, ax.get_ylim()[1])
    ax.set_title("Time during quench is scaled")
    ax.vlines((0, t_q_s), ax.get_ylim()[0], ax.get_ylim()[1], linestyles="dashed", color="grey", alpha=0.5)
    ax.set_ylabel(r"$\xi_y$")
    ax.set_xlabel(r"$t \qquad \qquad \qquad \qquad \qquad \qquad t/\tau \qquad \qquad  \qquad \qquad \qquad \qquad t $")
    fig.savefig(simulation_path + "/xiy-process.png", format="png", dpi=300)
    plt.show()

    if quench:
        figx, axx = plt.subplots(1, 1)
        figy, axy = plt.subplots(1, 1)
        axx.set_yscale("log")
        axy.set_yscale("log")
        axx.set_xscale("log")
        axy.set_xscale("log")
        tau_arr = []
        xix_after_quench_dic = {}
        xiy_after_quench_dic = {}
        for marker_ind, size in enumerate(size_y_dic):
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
                xix = xix_arr[t > t_after_quench][0]
                xiy = xiy_arr[t > t_after_quench][0]
                xix_after_quench_dic[size].append(xix)
                xiy_after_quench_dic[size].append(xiy)
                if marker_ind == 0:
                    tau_arr.append(float(tau))
            tau_arr = np.array(tau_arr)
            # do the fitting
            popt_x, _ = curve_fit(linear_fit, np.log(tau_arr), np.log(np.array(xix_after_quench_dic[size])))
            popt_y, _ = curve_fit(linear_fit, np.log(tau_arr), np.log(np.array(xiy_after_quench_dic[size])))

            axx.plot(tau_arr, xix_after_quench_dic[size], linestyle="",     color="C" + str(marker_ind), marker=markers[marker_ind], label=f"{size}")
            axy.plot(tau_arr, xiy_after_quench_dic[size], linestyle="",     color="C" + str(marker_ind), marker=markers[marker_ind], label=f"{size}")
            axx.plot(tau_arr, poly(tau_arr, popt_x[0], np.exp(popt_x[1])),  color="C" + str(marker_ind), label=r"$\frac{\nu}{1 + \nu z} =$" + f"{popt_x[0]:.2f}" )
            axy.plot(tau_arr, poly(tau_arr, popt_y[0], np.exp(popt_y[1])),  color="C" + str(marker_ind), label=r"$\frac{\nu}{1 + \nu z} =$" + f"{popt_y[0]:.2f}" )

    axx.set_xlabel(r"$\tau$")
    axx.set_ylabel(r"$\xi_x$")
    axy.set_xlabel(r"$\tau$")
    axy.set_ylabel(r"$\xi_y$")

    configure_ax(figx, axx)
    configure_ax(figy, axy)

    figx.savefig(simulation_path + "/xix-quench-scaling.png", format="png", dpi=300)
    figy.savefig(simulation_path + "/xiy-quench-scaling.png", format="png", dpi=300)


    plt.show()
if __name__ == "__main__":
    main()