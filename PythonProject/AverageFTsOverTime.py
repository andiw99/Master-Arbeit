import numpy as np

from FunctionsAndClasses import *

def main():
    simulation_path = "../../Generated content/Silicon/Quench/Left-Right"
    cut_zero_impuls = True
    quench = True
    scale_time = True
    t_xix = {}
    t_xiy = {}
    max_nr_curves = 5
    quench_zoom = 10
    y_scale = "log"
    y_lower_lim = 0.05

    for size in os.listdir(simulation_path):
        if size != "plots":
            sizepath = os.path.join(simulation_path, size)
            if os.path.isdir(sizepath):
                for setting in os.listdir(sizepath):
                    settingpath = os.path.join(sizepath, setting)
                    if os.path.isdir(settingpath):
                        t_xix[setting] = {}
                        t_xiy[setting] = {}
                        ft_k, ft_l = average_ft(settingpath)
                        if setting == "128.000000":
                            print(ft_l)
                        for t in ft_k:
                            p_k = get_frequencies_fftw_order(len(ft_k[t]))
                            if cut_zero_impuls:
                                p_k, ft_k[t] = cut_zero_imp(p_k, ft_k[t])
                            popt_x, perr_x = fit_lorentz(p_k, ft_k[t])
                            xix = np.abs(popt_x[0])
                            t_xix[setting][t] = xix
                        for t in ft_l:
                            p_l = get_frequencies_fftw_order(len(ft_l[t]))
                            if cut_zero_impuls:
                                p_l, ft_l[t] = cut_zero_imp(p_l, ft_l[t])
                            popt_y, perr_y = fit_lorentz(p_l, ft_l[t])
                            xiy = np.abs(popt_y[0])
                            t_xiy[setting][t] = xiy

    print(t_xiy)
    print("\n\n")
    print(t_xiy.keys())


    t_eq = 0
    if scale_time:
        parapath = find_first_txt_file(simulation_path)
        parameters = read_parameters_txt(parapath)
        t_eq = parameters["equil_time"]

    setting_var = "T"
    if quench:
        setting_var = r"$\tau$"

    fig, ax = plt.subplots(1, 1)
    ax.set_yscale(y_scale)
    valid_inds = np.linspace(0, len(t_xix) - 1, len(t_xix), dtype=int)
    if max_nr_curves:
        valid_inds = np.linspace(0, len(t_xix) - 1, max_nr_curves, dtype=int)
    keys = np.array(list(t_xix.keys()))[np.argsort([float(key) for key in t_xix.keys()])]
    for i, setting in enumerate(keys):
        print(setting, " with i = ", i , " in valid inds:", i in valid_inds)
        if i in valid_inds:
            t = list(t_xix[setting].keys())
            if quench & scale_time:
                tau = float(setting)
                t = rescale_t(t, tau, t_eq, quench_zoom)
            print(setting)
            xix = np.array(list(t_xix[setting].values()))
            if i == 0:
                pre_equil_xi = np.mean(xix[t < 0])
            ax.plot(t, xix, linestyle="", marker=".", ms=2, label=rf"$\xi_x$  {setting_var} = {float(setting):.2f}")
    configure_ax(fig, ax)
    if y_scale == "log":
        ax.set_ylim(0.8 * pre_equil_xi, ax.get_ylim()[1])
    plt.show()

    fig, ax = plt.subplots(1, 1)
    ax.set_yscale(y_scale)
    valid_inds = np.linspace(0, len(t_xiy) - 1, len(t_xiy), dtype=int)
    if max_nr_curves:
        valid_inds = np.linspace(0, len(t_xiy) - 1, max_nr_curves, dtype=int)
    keys = np.array(list(t_xiy.keys()))[np.argsort([float(key) for key in t_xiy.keys()])]
    for i, setting in enumerate(keys):
        t = list(t_xix[setting].keys())
        if i in valid_inds:
            if quench & scale_time:
                tau = float(setting)
                t = rescale_t(t, tau, t_eq, quench_zoom)
            xiy = np.array(list(t_xiy[setting].values()))
            if i == 0:
                pre_equil_xi = np.mean(xiy[t < 0])
            ax.plot(t, xiy, linestyle="", marker=".", ms=3, label=rf"$\xi_y$  {setting_var} = {float(setting):.2f}")
    configure_ax(fig, ax)
    if y_scale == "log":
        ax.set_ylim(pre_equil_xi * 0.8, ax.get_ylim()[1])
    plt.show()



if __name__ == "__main__":
    main()