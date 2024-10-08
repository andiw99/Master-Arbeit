from matplotlib import pyplot as plt
from FunctionsAndClasses import *
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.ticker as ticker


def analyze(df, parameters=None, savepath="./structfact.png", cutoff=np.pi/2, fitfunc=lorentzian, errors_for_fit=True,
            plot_struct = False, cut_zero_impuls = False, cut_around_peak=False):

    ft_avg_x, ft_avg_y, px, py, x_error, y_error = prepare_data(cut_zero_impuls, cutoff, df, cut_around_peak)
    # sorting
    #ft_avg_x = ft_avg_x[np.argsort(px)]
    #ft_avg_x = ft_avg_x[np.argsort(px)]
    #px = np.sort(px)
    #ft_avg_y = ft_avg_y[np.argsort(py)]
    #py = np.sort(py)
    #px = ((px + 2 * max(px)) % (2 * max(px))) - max(px)
    #py = ((py + 2 * max(px)) % (2 * max(py))) - max(py)

    print(px)
    print(ft_avg_y)

    if errors_for_fit:
        popt_x, perr_x = fit_lorentz(px, ft_avg_y, fitfunc, y_error)
        popt_y, perr_y = fit_lorentz(py, ft_avg_x, fitfunc, x_error)
    else:
        popt_x, perr_x = fit_lorentz(px, ft_avg_y, fitfunc, None)
        popt_y, perr_y = fit_lorentz(py, ft_avg_x, fitfunc, None)
    #print("a = %g" % popt_x[0])
    #print("x0 = %g" % popt_x[1])
    #print("gamma = %g" % popt_x[2])


    # xix = 1 / (np.abs(popt_x[2]) * 2)
    # xiy = 1/ (np.abs(popt_y[2]) * 2)
    xix = np.abs(popt_x[0])
    xiy = np.abs(popt_y[0])

    # plotting
    if not parameters:
        T = 0
    else:
        T = parameters["T"]
    if plot_struct:
        print("struct.fact file ")
        print(ft_avg_x)
        fig, axes = plot_struct_func(px, py,ft_avg_y, ft_avg_x, y_error, x_error)
        axk, axl = axes
        axk_ylims = axk.get_ylim()
        axl_ylims = axl.get_ylim()
        p = np.linspace(min(px), max(px), np.array(px).size)
        try:
            lorentz_x = fitfunc(p, *popt_x)
            lorentz_y = fitfunc(p, *popt_y)
        except:
            popt_x = (0, 0, 0)
            popt_y = (0, 0, 0)
            lorentz_x = fitfunc(p, *(popt_x))
            lorentz_y = fitfunc(p, *(popt_y))
        axes[0].plot(p, lorentz_x, label="Lorentzian fit")
        axes[1].plot(p, lorentz_y, label="Lorentzian fit")
        axk.set_ylim(axk_ylims)
        axl.set_ylim(axl_ylims)
        axes[0].set_title(rf"$\xi_x = {xix:.2f} \quad T = {T:2f}$")
        axes[1].set_title(rf"$\xi_y = {xiy:.2f}\quad T = {T:2f}$")
        configure_ax(fig, axes[0])
        configure_ax(fig, axes[1])
        plt.savefig(savepath, format="png")
        plt.show()
    #print("FWHM x:", np.abs(popt_x[2]) * 2)
    #print("FWHM y:", np.abs(popt_y[2]) * 2)
    #print("Corr Length x:", xix)
    #print("Corr Length y:", xiy)
    try:
        if perr_x.any():
            xix_err = perr_x[0]
            xiy_err = perr_y[0]
        else:
            xix_err = None
            xiy_err = None
    except:
        xix_err = None
        xiy_err = None

    return xix, xiy, xix_err, xiy_err


def prepare_data(cut_zero_impuls, cutoff, df, cut_around_peak=False, peak_cut_threshold=0.2, min_points_fraction=0.2, set_Fts_to_zero=False):
    px = df["px"]
    px = np.array(px)
    indices = [(-cutoff < x) & (x < cutoff) for x in px]
    #cutoff
    px = px[indices]
    px = px[~np.isnan(px)]
    ft_avg_y = np.array(df["ft_avg_y"])[indices]
    ft_avg_y = ft_avg_y[~np.isnan(ft_avg_y)]
    py = df["py"]
    py = np.array(py)[indices]
    py = py[~np.isnan(py)]
    ft_avg_x = np.array(df["ft_avg_x"])[indices]
    ft_avg_x = ft_avg_x[~np.isnan(ft_avg_x)]

    # reorder...
    ft_avg_x = np.concatenate((ft_avg_x[len(ft_avg_x)//2:], ft_avg_x[:len(ft_avg_x)//2]))
    ft_avg_y = np.concatenate((ft_avg_y[len(ft_avg_y)//2:], ft_avg_y[:len(ft_avg_y)//2]))

    ft_avg_x, py = prepare_fit_data(cut_around_peak, cut_zero_impuls, ft_avg_x,
                                     peak_cut_threshold, set_Fts_to_zero,
                                 min_points_fraction)
    ft_avg_y, px = prepare_fit_data(cut_around_peak, cut_zero_impuls, ft_avg_y,
                                     peak_cut_threshold, set_Fts_to_zero,
                                 min_points_fraction)


    y_error = None
    x_error = None
    #if cut_zero_impuls:
        #ft_avg_y = ft_avg_y[px != 0]
        #y_error = y_error[px != 0]
        #px = px[px != 0]
        #ft_avg_x = ft_avg_x[py != 0]
        #x_error = x_error[py != 0]
        #py = py[py != 0]
    return ft_avg_x, ft_avg_y, px, py, x_error, y_error


def main():
    # parameters
    root = "../../Generated content/Silicon/Subsystems/Suite/L_xi/Check-OBC/0.4161791450287818/Tc/OBC/160"
    root = "../../Generated content/Silicon/Subsystems/Suite/L_xi/scan-more-flips-more-vals/0.4161791450287818/Tc/160"
    name = "struct.fact"
    png_name = "struct.fact-fit2"
    root_dirs = os.listdir(root)
    cutoff =  np.pi
    fitfunc = lorentz_offset
    errors_for_fit = False
    plot_struct = True
    cut_zero_impuls = True
    cut_around_peak = False
    nu_est = 0.8
    T_c_est = 0.7
    print(root_dirs)
    # arrays to save the xi corrsponding to T
    T_arr = []
    xix_arr = []
    xix_err_arr = []
    xiy_arr = []
    xiy_err_arr = []
    xi_arr = []
    xi_err_arr = []
    fwhm = False
    # Loop through the directory contents and print the directories
    for item in root_dirs:

        if (item != "plots"):
            # Create the full path to the item
            dir_path = os.path.join(root, item)

            # Check if the item is a directory
            if os.path.isdir(dir_path) & (dir_path != root + "plots"):
                filename = dir_path + "/" + name
                print("reading: ", filename)
                files = os.listdir(dir_path)
                parameters = {}
                for f in files:
                    # we take the first file to be the parameters
                    if(os.path.splitext(f)[1] == ".txt"):
                        parameters = read_parameters_txt(os.path.join(dir_path, f))
                        T = parameters["T"]
                df = read_struct_func(filename)
                csv_file = find_first_csv_file(dir_path)
                struct_fact_k, struct_fact_l = calc_structure_factor(csv_file)
                st_fact_k = struct_fact_k[np.max(list(struct_fact_k.keys()))]
                st_fact_l = struct_fact_l[np.max(list(struct_fact_l.keys()))]
                print(st_fact_l)
                print(st_fact_k)

                xix_2nd = get_2nd_moment_corr_length(st_fact_l)
                xiy_2nd = get_2nd_moment_corr_length(st_fact_k)

                ft_k, k = prepare_fit_data(cut_around_peak, cut_zero_impuls, st_fact_k,
                                     0.2, False,
                                           0.2, nr_additional_cuts=0)
                ft_k = ft_k[np.argsort(k)]
                k = np.sort(k)
                #k = 4 * k
                ft_l, l = prepare_fit_data(cut_around_peak, cut_zero_impuls, st_fact_l,
                                     0.2, False,0.2, nr_additional_cuts=0)
                ft_l = ft_l[np.argsort(l)]
                l = np.sort(l)
                #l = 4 * l
                popt_x, perr_x = fit_lorentz(k, ft_k, fitfunc)
                popt_y, perr_y = fit_lorentz(l, ft_l, fitfunc)
                xix = np.abs(popt_x[0])
                xiy = np.abs(popt_y[0])
                print("\n", st_fact_k)
                print(k)
                fig, ax = plot_struct_func(l, k, ft_l, ft_k)
                #ax[0].set_ylim(0, 1.5)
                try:
                    lorentz_x = fitfunc(k, *popt_x)
                    lorentz_y = fitfunc(l, *popt_y)
                    ax[1].plot(k, lorentz_x, label="Lorentzian fit")
                    ax[0].plot(l, lorentz_y, label="Lorentzian fit")
                except:
                    pass
                ax[0].set_title(rf"$\xi_x = {xiy:.2f} \quad T = {T:2f} \quad second \xi_x = {xix_2nd}$")
                ax[1].set_title(rf"$\xi_y = {xix:.2f}\quad T = {T:2f} \quad second \xi_y = {xiy_2nd}$")
                configure_ax(fig, ax[0])
                configure_ax(fig, ax[1])
                plt.show()
                try:
                    t_ft_l, t_ft_k = average_ft(dir_path)
                    st_fact_k = t_ft_k[np.max(list(t_ft_k.keys()))]
                    st_fact_l = t_ft_l[np.max(list(t_ft_l.keys()))]

                    ft_k, k = prepare_fit_data(cut_around_peak, cut_zero_impuls, st_fact_k,
                                         0.2, False,0.2)
                    ft_l, l = prepare_fit_data(cut_around_peak, cut_zero_impuls, st_fact_l,
                                         0.2, False,0.2)
                    print(".ft file;:")
                    print("\n", st_fact_k)
                    print(k)
                    fig, ax = plot_struct_func(l, k, ft_l, ft_k)
                    #ax[0].set_ylim(0, 100)
                    configure_ax(fig, ax[0])
                    configure_ax(fig, ax[1])
                    plt.show()
                except:
                    pass


                if not parameters:
                    T = 0
                else:
                    T = parameters["T"]

                if fwhm:
                    ft_avg_x, ft_avg_y, px, py, x_error, y_error = prepare_data(cut_zero_impuls, cutoff, df)

                    xix = 1 / find_fwhm(px, ft_avg_y)
                    xiy = 1 / find_fwhm(py, ft_avg_x)
                    xix_err = 0
                    xiy_err = 0
                else:
                    xix, xiy, xix_err,\
                        xiy_err = analyze(df, parameters,
                                                      savepath=dir_path + png_name,
                                                      cutoff=cutoff, fitfunc=fitfunc,
                                                      errors_for_fit=errors_for_fit,
                                                      plot_struct=plot_struct,
                                                      cut_zero_impuls=cut_zero_impuls,
                                                        cut_around_peak=cut_around_peak)

                xi = np.abs(1 / 2 * (xix + xiy))
                if xix_err == None:
                    xi_err = None
                else:
                    xi_err = 1 / 2 * (xix_err + xiy_err)

                T_arr.append(T)
                xix_arr.append(xix)
                xiy_arr.append(xiy)
                xix_err_arr.append(xix_err)
                xiy_err_arr.append(xiy_err)
                xi_arr.append(xi)
                xi_err_arr.append(xi_err)


    xix_sorted = np.array(xix_arr)[np.argsort(T_arr)]
    xiy_sorted = np.array(xiy_arr)[np.argsort(T_arr)]
    xi_sorted = np.array(xi_arr)[np.argsort(T_arr)]
    xix_err_sorted = np.array(xix_err_arr)[np.argsort(T_arr)]
    xiy_err_sorted = np.array(xiy_err_arr)[np.argsort(T_arr)]
    xi_err_sorted = np.array(xi_err_arr)[np.argsort(T_arr)]
    T_arr = np.sort(T_arr)

    # okay fitting the amplitudes
    # lets say criticical Temperature is just the maximum of the correlation length

    T_c = T_arr[np.argmax(xi_sorted)]
    T_c = 0.95

    eps_array = (T_arr - T_c) / T_c
    print(eps_array)
    xix_fit = xix_sorted[eps_array > 0]
    xiy_fit = xiy_sorted[eps_array > 0]
    eps_fit = eps_array[eps_array > 0]
    popt_x, _ = curve_fit(critical_amplitude, eps_fit, xix_fit)
    popt_y, _ = curve_fit(critical_amplitude, eps_fit, xiy_fit)

    xi0_x = popt_x[0]
    xi0_y = popt_y[0]

    print(popt_x)
    print(popt_y)

    print("xi_x / xi_y = ", xi0_x / xi0_y)
    # plotting
    fig, ax = plt.subplots(1, 1)
    # Setze Tickmarken und Labels
    ax.tick_params(direction='in', which='both', length=6, width=2, labelsize=9)
    ax.tick_params(direction='in', which='minor', length=3, width=1, labelsize=9)
    span = np.maximum(np.max(T_arr) - np.min(T_arr), 0.05)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=span / 4))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(base=span / 4 / 5))
    # TODO minor locator muss
    #ax.yaxis.set_minor_locator((plt.MultipleLocator(0.2)))
    # Füge Gitterlinien hinzu
    ax.grid(which='major', linestyle='--', alpha=0.5)
    print(T_arr)
    eps = 10 ** (-4)
    ax.errorbar(T_arr, xix_sorted, yerr=xix_err_sorted, ls="", marker="x", color="C0", ecolor="black", capsize=3)
    ax.errorbar(T_arr, xiy_sorted, yerr=xiy_err_sorted, ls="", marker="x", color="C1", ecolor="black", capsize=3)
    print(critical_amplitude(eps_array, xi0_x))
    pre_ylim_min = ax.get_ylim()[0]
    T_plot = np.linspace(T_c + eps, np.max(T_arr), 200)
    eps_plot = (T_plot - T_c) / T_c
    ax.plot(T_plot, critical_amplitude(eps_plot, xi0_x), color="C0", label=rf"$\xi^+_x = {xi0_x:.2f}$")
    ax.plot(T_plot, critical_amplitude(eps_plot, xi0_y), color="C1", label=rf"$\xi^+_y = {xi0_y:.2f}$")
    ax.vlines(T_c, pre_ylim_min, np.max(xix_sorted) * 1.05, linestyles="dashed", alpha=0.2, colors="black", label=rf"$T_c  = ${T_c}")
    ax.plot([], [], label = rf"$\xi^+_x / \xi^+_y = {(xi0_x / xi0_y):.2f}$", linestyle="")
    ax.set_xlabel("T")
    ax.set_ylabel(r"$\xi(T)$")
    ax.set_ylim((pre_ylim_min, np.max(xix_sorted) * 1.05))
    ax.set_title("Corr Length depending on T")
    ax.legend()
    save_plot(root, "/xix-xiy.png")
    # plotting xi
    fig, ax = plt.subplots(1, 1)
    # Setze Tickmarken und Labels
    ax.tick_params(direction='in', which='both', length=6, width=2, labelsize=9)
    ax.tick_params(direction='in', which='minor', length=3, width=1, labelsize=9)

    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=span / 4))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(base=span / 4 / 5))
    # TODO minor locator muss
    # ax.yaxis.set_minor_locator((plt.MultipleLocator(0.2)))
    # Füge Gitterlinien hinzu
    ax.grid(which='major', linestyle='--', alpha=0.5)
    ax.errorbar(T_arr, xi_sorted, yerr=xi_err_sorted, ls="", marker="x", color="C0", ecolor="black", capsize=3)
    ax.set_xlabel("T")
    ax.set_ylabel(r"$\xi(T)$")
    ax.set_title("Corr Length depending on T")
    save_plot(root, "/xi.png")

    eps = (T_arr - T_c_est) / T_c_est
    Tg = np.linspace(T_c_est, np.max(T_arr), 100)[1:]
    Tl = np.linspace(T_c_est, np.min(T_arr), 100)[1:]
    scaling_right = corr_scaling_right(Tg, T_c_est, nu_est, 0.2)
    scaling_left = corr_scaling_left(Tl, T_c_est, nu_est, 2)
    ax.set_ylim(0, np.max(xi_arr) + 0.1 * np.max(xi_arr))
    # ax.plot(Tg, scaling_right)
    # ax.plot(Tl, scaling_left)

    fig, ax = plt.subplots(1, 1)
    xix_inv = 1 / xix_sorted
    xiy_inv = 1 / xiy_sorted
    ax.plot(T_arr, xix_inv, ls="", marker="o", ms=4, fillstyle="none", color="C0")
    ax.plot(T_arr, xiy_inv, ls="", marker="o", ms=4, fillstyle="none", color="C1")
    ax.set_xlabel("T")
    ax.set_ylabel(r"$\frac{1}{\xi(T)}$")
    configure_ax(fig, ax)
    plt.savefig(root + "/1_xi.png")

    fig, ax = plt.subplots(1, 1)
    xi_inv = 1 / xi_sorted
    ax.plot(np.log(T_arr), np.log(xix_inv), ls="", marker="o", ms=4, fillstyle="none", color="C0")
    ax.set_xlabel("T")
    ax.set_ylabel(r"$\frac{1}{\xi(T)}$")
    configure_ax(fig, ax)
    plt.savefig(root + "/1_xi_avg_log.png")

    plt.show()



if __name__ == "__main__":
    main()