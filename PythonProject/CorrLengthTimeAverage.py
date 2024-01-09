from matplotlib import pyplot as plt
from FunctionsAndClasses import *
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.ticker as ticker





def main():
    # parameters
    simulation_folder = ("../../Generated content/Silicon/Amplitude/Small/")
    name = "struct.fact"
    png_name = "struct.fact-fit2"
    root_dirs = os.listdir(simulation_folder)
    nu_est = 0.8
    T_c_est = 0.9
    cut_zero_impuls = False
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
                size_result_x = process_size_folder(size_folder_path, threshold, key='T', value="xix", file_ending=".corr")
                size_result_y = process_size_folder(size_folder_path, threshold, key='T', value="xiy", file_ending=".corr")

                T_arr, xix_ft_arr, xiy_ft_arr = process_size_folder_ft(cut_zero_impuls, simulation_folder,
                                                                       size_folder_path, threshold)

                try:
                    results_x[int(size_folder)] = size_result_x
                    results_y[int(size_folder)] = size_result_y
                    results_ft_x[int(size_folder)] = (T_arr, xix_ft_arr)
                    results_ft_y[int(size_folder)] = (T_arr, xiy_ft_arr)

                except ValueError:
                    pass
    print(results_x)
    print(results_y)

    for ((keyx, valuex), (keyy, valuey)) in zip(results_x.items(), results_y.items()):
        xix_sorted = np.array(valuex["xix"])[np.argsort(valuex['T'])]
        xiy_sorted = np.array(valuey["xiy"])[np.argsort(valuey['T'])]
        xix_err_sorted = np.zeros_like(xix_sorted)
        xiy_err_sorted = np.zeros_like(xiy_sorted)
        T_arr = np.sort(valuex['T'])


        # okay fitting the amplitudes
        # lets say criticical Temperature is just the maximum of the correlation length

        T_c = T_c_est

        eps_array = (T_arr - T_c) / T_c
        print(eps_array)
        xix_fit = xix_sorted[eps_array > 0]
        print(xix_fit)
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
        save_plot(simulation_folder, f"/xix-xiy-time-average-{keyx}.png")

        plt.show()

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

        T_c = T_c_est

        eps_array = (T_arr - T_c) / T_c
        print(eps_array)
        xix_fit = xix_sorted[eps_array > 0]
        print(xix_fit)
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
        # ax.yaxis.set_minor_locator((plt.MultipleLocator(0.2)))
        # Füge Gitterlinien hinzu
        ax.grid(which='major', linestyle='--', alpha=0.5)
        print(T_arr)
        eps = 10 ** (-4)
        ax.errorbar(T_arr, xix_sorted, yerr=xix_err_sorted, ls="", marker="x", color="C0", ecolor="black",
                    capsize=3)
        ax.errorbar(T_arr, xiy_sorted, yerr=xiy_err_sorted, ls="", marker="x", color="C1", ecolor="black",
                    capsize=3)
        print(critical_amplitude(eps_array, xi0_x))
        pre_ylim_min = ax.get_ylim()[0]
        T_plot = np.linspace(T_c + eps, np.max(T_arr), 200)
        eps_plot = (T_plot - T_c) / T_c
        ax.plot(T_plot, critical_amplitude(eps_plot, xi0_x), color="C0", label=rf"$\xi^+_x = {xi0_x:.2f}$")
        ax.plot(T_plot, critical_amplitude(eps_plot, xi0_y), color="C1", label=rf"$\xi^+_y = {xi0_y:.2f}$")
        ax.vlines(T_c, pre_ylim_min, np.max(xix_sorted) * 1.05, linestyles="dashed", alpha=0.2, colors="black",
                  label=rf"$T_c  = ${T_c}")
        ax.plot([], [], label=rf"$\xi^+_x / \xi^+_y = {(xi0_x / xi0_y):.2f}$", linestyle="")
        ax.set_xlabel("T")
        ax.set_ylabel(r"$\xi(T)$")
        ax.set_ylim((pre_ylim_min, np.max(xix_sorted) * 1.05))
        ax.set_title("Corr Length depending on T")
        ax.legend()
        save_plot(simulation_folder, f"/xix-xiy-ft-time-average-{keyx}.png")

        plt.show()

def process_size_folder_ft(cut_zero_impuls, simulation_folder, size_folder_path, threshold):
    xix_ft_arr = []
    xiy_ft_arr = []
    T_arr = []
    for temp in os.listdir(size_folder_path):
        if (temp != "plots") & (temp[0] != "."):
            settingpath = os.path.join(size_folder_path, temp)
            if os.path.isdir(settingpath):
                parapath = find_first_txt_file(simulation_folder)
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