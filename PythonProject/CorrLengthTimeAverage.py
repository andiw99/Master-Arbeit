from matplotlib import pyplot as plt
from FunctionsAndClasses import *
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.ticker as ticker



def inv_xi_fit(T, xi, Tc):
    eps = (T - Tc) / Tc
    return eps / xi


def main():
    # parameters
    simulation_folder = ("../../Generated content/Silicon/Amplitude/Small/")
    name = "struct.fact"
    png_name = "struct.fact-fit2"
    root_dirs = os.listdir(simulation_folder)
    nu_est = 0.8
    T_c_est_min = 0.95
    T_c_est_max = 0.956
    T_min_fit_min = 0.96
    T_min_fit_max = 0.9625
    res = 100
    cut_zero_impuls = True
    print(root_dirs)

    threshold = 50000  # Example threshold value, adjust as needed
    max_L_fit = 100
    transparent_plots = False
    linewidth = 1

    results_x = {}
    results_y = {}

    results_ft_x = {}
    results_ft_y = {}

    for size_folder in os.listdir(simulation_folder):
        if (size_folder != "plots") & (size_folder[0] != "."):
            size_folder_path = os.path.join(simulation_folder, size_folder)
            if os.path.isdir(size_folder_path):
                #size_result_x = process_size_folder(size_folder_path, threshold, key='T', value="xix", file_ending=".corr")
                #size_result_y = process_size_folder(size_folder_path, threshold, key='T', value="xiy", file_ending=".corr")

                T_arr, xix_ft_arr, xiy_ft_arr = process_size_folder_ft(cut_zero_impuls,
                                                                       size_folder_path, threshold)

                try:
                    #results_x[int(size_folder)] = size_result_x
                    #results_y[int(size_folder)] = size_result_y
                    results_ft_x[int(size_folder)] = (T_arr, xix_ft_arr)
                    results_ft_y[int(size_folder)] = (T_arr, xiy_ft_arr)

                except ValueError:
                    pass
    print(results_x)
    print(results_y)

    # for ((keyx, valuex), (keyy, valuey)) in zip(results_x.items(), results_y.items()):
    #     xix_sorted = np.array(valuex["xix"])[np.argsort(valuex['T'])]
    #     xiy_sorted = np.array(valuey["xiy"])[np.argsort(valuey['T'])]
    #     xix_err_sorted = np.zeros_like(xix_sorted)
    #     xiy_err_sorted = np.zeros_like(xiy_sorted)
    #     T_arr = np.sort(valuex['T'])
    #
    #
    #     # okay fitting the amplitudes
    #     # lets say criticical Temperature is just the maximum of the correlation length
    #
    #     T_c = 1/2 * (T_c_est_min + T_c_est_max)
    #
    #     eps_array = (T_arr - T_c) / T_c
    #     print(eps_array)
    #     xix_fit = xix_sorted[eps_array > 0]
    #     print(xix_fit)
    #     xiy_fit = xiy_sorted[eps_array > 0]
    #     eps_fit = eps_array[eps_array > 0]
    #     popt_x, _ = curve_fit(critical_amplitude, eps_fit, xix_fit)
    #     popt_y, _ = curve_fit(critical_amplitude, eps_fit, xiy_fit)
    #
    #     xi0_x = popt_x[0]
    #     xi0_y = popt_y[0]
    #
    #     print(popt_x)
    #     print(popt_y)
    #
    #     print("xi_x / xi_y = ", xi0_x / xi0_y)
    #     # plotting
    #     fig, ax = plt.subplots(1, 1)
    #     # Setze Tickmarken und Labels
    #     ax.tick_params(direction='in', which='both', length=6, width=2, labelsize=9)
    #     ax.tick_params(direction='in', which='minor', length=3, width=1, labelsize=9)
    #     span = np.maximum(np.max(T_arr) - np.min(T_arr), 0.05)
    #     ax.xaxis.set_major_locator(ticker.MultipleLocator(base=span / 4))
    #     ax.xaxis.set_minor_locator(ticker.MultipleLocator(base=span / 4 / 5))
    #     # TODO minor locator muss
    #     #ax.yaxis.set_minor_locator((plt.MultipleLocator(0.2)))
    #     # Füge Gitterlinien hinzu
    #     ax.grid(which='major', linestyle='--', alpha=0.5)
    #     print(T_arr)
    #     eps = 10 ** (-4)
    #     ax.errorbar(T_arr, xix_sorted, yerr=xix_err_sorted, ls="", marker="x", color="C0", ecolor="black", capsize=3)
    #     ax.errorbar(T_arr, xiy_sorted, yerr=xiy_err_sorted, ls="", marker="x", color="C1", ecolor="black", capsize=3)
    #     print(critical_amplitude(eps_array, xi0_x))
    #     pre_ylim_min = ax.get_ylim()[0]
    #     T_plot = np.linspace(T_c + eps, np.max(T_arr), 200)
    #     eps_plot = (T_plot - T_c) / T_c
    #     ax.plot(T_plot, critical_amplitude(eps_plot, xi0_x), color="C0", label=rf"$\xi^+_x = {xi0_x:.2f}$")
    #     ax.plot(T_plot, critical_amplitude(eps_plot, xi0_y), color="C1", label=rf"$\xi^+_y = {xi0_y:.2f}$")
    #     ax.vlines(T_c, pre_ylim_min, np.max(xix_sorted) * 1.05, linestyles="dashed", alpha=0.2, colors="black", label=rf"$T_c  = ${T_c}")
    #     ax.plot([], [], label = rf"$\xi^+_x / \xi^+_y = {(xi0_x / xi0_y):.2f}$", linestyle="")
    #     ax.set_xlabel("T")
    #     ax.set_ylabel(r"$\xi(T)$")
    #     ax.set_ylim((pre_ylim_min, np.max(xix_sorted) * 1.05))
    #     ax.set_title(r"Corr Length depending on T, $\xi$ on run")
    #     ax.plot([], [], linestyle="", marker="", label=rf"$L_x = {keyx}$")
    #     ax.legend()
    #     save_plot(simulation_folder, f"/xix-xiy-time-average-{keyx}.png")
    #
    #     plt.show()

    for ((sizex, valuex), (sizey, valuey)) in zip(results_ft_x.items(),
                                                  results_ft_y.items()):
        # I want to plot the 1 / correlation lengths
        # sizex and sizey are what the names suggest
        # value_x should be tuples of (T_arr, xix_array)
        T = valuex[0]
        xix_inv_sorted = 1 / np.array(valuex[1])[np.argsort(T)]
        xiy_inv_sorted = 1 / np.array(valuey[1])[np.argsort(T)]
        T_arr = np.sort(T)

        # custom interval to get a feeling of r² value:
        # T_min_fit = 0.95
        # T_max_fit = 1.05
        # xix_inv_fit = xix_inv_sorted[(T_arr > T_min_fit)  & (T_arr < T_max_fit)]
        # T_fit = T_arr[(T_arr > T_min_fit)  & (T_arr < T_max_fit)]

        # reg = linregress(T_fit, xix_inv_fit)
        # print(f"custom r² {reg.rvalue ** 2}")


        # We look for the most linear region lol
        Tc_guess = 0.95     # this will be mandatory in the future
        reg_x, T_include_start_x, T_include_end_x = best_fit_inv(T_arr, xix_inv_sorted, Tc_guess, tolerance=0.05)
        print(f"x: Included Temperatures between {T_arr[T_include_start_x]} and {T_arr[T_include_end_x]}.")
        print(f"The r²-value is r² = {reg_x.rvalue ** 2}")
        reg_y, T_include_start_y, T_include_end_y = best_fit_inv(T_arr,
                                                                 xiy_inv_sorted, Tc_guess, tolerance=0.05)
        print(
            f"y: Included Temperatures between {T_arr[T_include_start_y]} and {T_arr[T_include_end_y]}.")
        print(f"The r²-value is r² = {reg_y.rvalue ** 2}")
        # from the regression parameter we can get the estimated Tc
        xix_ampl = 1 / reg_x.slope
        xiy_ampl = 1 / reg_y.slope
        Tc_x = - reg_x.intercept * xix_ampl
        Tc_y = - reg_y.intercept * xiy_ampl

        #T_min = 0.96        # for now we have a manual T_min for the fit
        # cut the plotting data
        # xix_inv_fit = xix_inv_sorted[T_arr > T_min]
        # xiy_inv_fit = xiy_inv_sorted[T_arr > T_min]
        # T_fit = T_arr[T_arr > T_min]
        #
        # # Fit the data
        # popt_x, _ = curve_fit(inv_xi_fit, T_fit, xix_inv_fit, p0=(1, 0.95), maxfev=10000)
        # xix_ampl = popt_x[0]
        # Tc_x = popt_x[1]
        # popt_y, _ = curve_fit(inv_xi_fit, T_fit, xiy_inv_fit, p0=(1, 0.95), maxfev=10000)
        # xiy_ampl = popt_y[0]
        # Tc_y = popt_y[1]

        fig, axx = plt.subplots(1, 1)
        # Plot the data
        axx.plot(T_arr, xix_inv_sorted, label=rf"$1 / \xi_x$", linestyle="", marker="x")
        # Plot the fit
        axx.plot(T_arr, reg_x.intercept + reg_x.slope * T_arr,
                 label=rf"$\xi_x^+ = {xix_ampl:.2f}, T_c = {Tc_x:.3f}$")

        #axx.plot(T_arr, inv_xi_fit(T_arr, xix_ampl, Tc_x), label=rf"$\xi_x^+ = {xix_ampl:.2f}, T_c = {Tc_x:.3f}$")
        configure_ax(fig, axx)
        create_directory_if_not_exists(simulation_folder + "/plots")
        plt.savefig(simulation_folder + f"/plots/xix-inv-{sizex}.png", format="png")
        plt.show()

        fig, axy = plt.subplots(1, 1)
        axy.plot(T_arr, xiy_inv_sorted, label=rf"$1 / \xi_y$", linestyle="", marker="x")
        axy.plot(T_arr, reg_y.intercept + reg_y.slope * T_arr,
                 label=rf"$\xi_y^+ = {xiy_ampl:.2f}, T_c = {Tc_y:.3f}$")
        #axy.plot(T_arr, inv_xi_fit(T_arr, xiy_ampl, Tc_y),
        #         label=rf"$\xi_y^+ = {xiy_ampl:.2f}, T_c = {Tc_y:.3f}$")
        print(f"xi_x / xi_y = {xix_ampl / xiy_ampl}")
        configure_ax(fig, axy)
        plt.savefig(simulation_folder + f"/plots/xiy-inv-{sizey}.png", format="png")
        plt.show()

    exit()
    for ((sizex, valuex), (sizey, valuey)) in zip(results_ft_x.items(), results_ft_y.items()):
        # the values here are the T_arr, xix_arr pairs
        T = valuex[0]
        xix_sorted = np.array(valuex[1])[np.argsort(T)]
        xiy_sorted = np.array(valuey[1])[np.argsort(T)]
        xix_err_sorted = np.zeros_like(xix_sorted)
        xiy_err_sorted = np.zeros_like(xiy_sorted)
        T_arr = np.sort(T)

        # okay fitting the amplitudes
        # lets say criticical Temperature is just the maximum of the correlation length

        T_c_x, T_min_fit_x, eps_array_x, eps_fit_x, popt_x, xix_fit = best_fit(T_arr, T_c_est_max, T_c_est_min, T_min_fit_max,
                                                                       T_min_fit_min, res, xix_sorted)
        print("T_c_x = ", T_c_x)
        print("T_min_fit_x = ", T_min_fit_x)
        T_c_y, T_min_fit_y, eps_array_y, eps_fit_y, popt_y, xiy_fit = best_fit(T_arr, T_c_est_max, T_c_est_min, T_min_fit_max,
                                                                       T_min_fit_min, res, xiy_sorted)
        print("T_c_y = ", T_c_y)
        print("T_min_fit_y = ", T_min_fit_y)
        xi0_x = popt_x[0]
        xi0_y = popt_y[0]

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
        # ax.yaxis.set_minor_locator((plt.MultipleLocator(0.2)))
        # Füge Gitterlinien hinzu
        ax.grid(which='major', linestyle='--', alpha=0.5)
        print(T_arr)
        eps = 10 ** (-4)
        ax.errorbar(T_arr, xix_sorted, yerr=xix_err_sorted, ls="", marker="x", color="C0", ecolor="black",
                    capsize=3)
        ax.errorbar(T_arr, xiy_sorted, yerr=xiy_err_sorted, ls="", marker="x", color="C1", ecolor="black",
                    capsize=3)
        print(critical_amplitude(eps_array_x, xi0_x))
        pre_ylim_min = ax.get_ylim()[0]
        T_plot_x = np.linspace(T_c_x + eps, np.max(T_arr), 200)
        eps_plot_x = (T_plot_x - T_c_x) / T_c_x
        ax.plot(T_plot_x, critical_amplitude(eps_plot_x, xi0_x),
                color="C0", label=rf"$\xi^+_x = {xi0_x:.2f}$")
        ax.plot(T_plot_x, critical_amplitude(eps_plot_x, xi0_x + 0.5 * xi0_x), color="C0", label=rf"$\xi^+_x = {xi0_x:.2f}$")
        ax.plot(T_plot_x, critical_amplitude(eps_plot_x, xi0_x + 0.25 * xi0_x),
                color="C0", label=rf"$\xi^+_x = {xi0_x:.2f}$")
        T_plot_y = np.linspace(T_c_y + eps, np.max(T_arr), 200)
        eps_plot_y = (T_plot_y - T_c_y) / T_c_y
        ax.plot(T_plot_y, critical_amplitude(eps_plot_y, xi0_y), color="C1", label=rf"$\xi^+_y = {xi0_y:.2f}$")
        ax.vlines(T_c_x, pre_ylim_min, np.max(xix_sorted) * 1.05, linestyles="dashed", alpha=0.2, colors="C0",
                  label=rf"$T_c  = ${T_c_x}")
        ax.vlines(T_c_y, pre_ylim_min, np.max(xix_sorted) * 1.05, linestyles="dashed", alpha=0.2, colors="C1",
                  label=rf"$T_c  = ${T_c_y:.3f}")
        ax.plot([], [], label=rf"$\xi^+_x / \xi^+_y = {(xi0_x / xi0_y):.2f}$", linestyle="")
        ax.set_xlabel("T")
        ax.set_ylabel(r"$\xi(T)$")
        ax.set_ylim((pre_ylim_min, np.max(xix_sorted) * 1.05))
        ax.set_title("Corr Length depending on T, (ft) on run")
        ax.plot([], [], linestyle="", marker="", label=rf"$L_x = {sizex}$")
        ax.legend()
        save_plot(simulation_folder, f"/xix-xiy-ft-time-average-{sizex}.png")

        plt.show()




def best_fit(T_arr, T_c_est_max, T_c_est_min, T_min_fit_max, T_min_fit_min, res, xi_sorted):
    T_c_ests = np.linspace(T_c_est_min, T_c_est_max, num=res, endpoint=True)
    T_mins_fit = np.linspace(T_min_fit_min, T_min_fit_max, num=res, endpoint=True)
    min_mse = np.infty
    T_c_best = 0
    T_min_best = 0
    popt = ()
    for T_c in T_c_ests:
        for T_min_fit in T_mins_fit:
            if T_min_fit > T_c:
                eps_array = (T_arr - T_c) / T_c
                xi_fit = xi_sorted[(eps_array > 0) & (T_arr > T_min_fit)]
                eps_fit = eps_array[(eps_array > 0) & (T_arr > T_min_fit)]
                cur_popt, cur_pcov, infodict, _, _ = curve_fit(critical_amplitude, eps_fit, xi_fit, full_output=True)
                mse = np.mean(infodict["fvec"] ** 2)
                if mse < min_mse:
                    T_c_best = T_c
                    T_min_best = T_min_fit
                    popt = cur_popt
                    min_mse = mse
    return T_c_best, T_min_best, eps_array, eps_fit, popt, xi_fit


def process_size_folder_ft(cut_zero_impuls, size_folder_path, threshold):
    xix_ft_arr = []
    xiy_ft_arr = []
    T_arr = []
    for temp in os.listdir(size_folder_path):
        if (temp != "plots") & (temp[0] != "."):
            settingpath = os.path.join(size_folder_path, temp)
            if os.path.isdir(settingpath):
                parapath = find_first_txt_file(settingpath)
                parameters = read_parameters_txt(parapath)

                Lx = parameters["subsystem_Lx"]
                Ly = parameters["subsystem_Ly"]

                xix_arr = []
                xiy_arr = []

                ft_k, ft_l = average_ft(settingpath)
                for t in ft_k:
                    if t > threshold:
                        p_k = get_frequencies_fftw_order(len(ft_k[t]))
                        if cut_zero_impuls:
                            p_k, ft_k[t] = cut_zero_imp(p_k, ft_k[t])
                        popt_x, perr_x = fit_lorentz(p_k, ft_k[t])
                        xix = np.minimum(np.abs(popt_x[0]), Lx)
                        xix_arr.append(xix)
                for t in ft_l:
                    if t > threshold:
                        p_l = get_frequencies_fftw_order(len(ft_l[t]))
                        if cut_zero_impuls:
                            p_l, ft_l[t] = cut_zero_imp(p_l, ft_l[t])
                        popt_y, perr_y = fit_lorentz(p_l, ft_l[t])
                        xiy = np.minimum(np.abs(popt_y[0]), Ly)
                        xiy_arr.append(xiy)
                # time average of the correlation length for this specific temperature
                xix_ft = np.mean(xix_arr)
                xiy_ft = np.mean(xiy_arr)
                # add the current temp to the array
                T_arr.append(float(temp))
                # add the correlation lenghts to the arrays
                xix_ft_arr.append(xix_ft)
                xiy_ft_arr.append(xiy_ft)
    return T_arr, xix_ft_arr, xiy_ft_arr


if __name__ == "__main__":
    main()