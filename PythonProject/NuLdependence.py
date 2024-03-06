import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from FunctionsAndClasses import *
from scipy.optimize import curve_fit


def main():
    simulation_folder = '../../Generated content/Silicon/Subsystems/Suite/L/'
    threshold = 0.1
    transparent_plots = False
    linewidth = 1

    results = {}

    for size_folder in os.listdir(simulation_folder):
        if (size_folder[0] != ".") & (size_folder != "plots"):
            size_folder_path = os.path.join(simulation_folder, size_folder)
            if os.path.isdir(size_folder_path):
                size_result = process_size_folder(size_folder_path, threshold)
                results[int(size_folder)] = size_result

    sizes = sorted(list(results.keys()))
    # How do we know which sizes at the beginning were pairs? Do we just multiply it by two?
    results_pairs = {}   # this will be a dictionary that has a pair of sizes as keys
    for size in sizes:
        pair_size = 2 * size        # this is arbitrary
        if pair_size in sizes:
            pair = (size, pair_size)
            result_pair = {}        # this should be a dict like the big result but only for the two sizes so that we can use our usual functions
            for s in pair:
                result_pair[s] = results[s]
            results_pairs[pair] = result_pair

    # Now we should have a fancy dictionary
    T_cs = {}        # again dictionaries
    U_Ls = {}
    dU_dT = {}

    for pair in results_pairs:
        size1, size2 = pair
        results_pair = results_pairs[pair]
        T_c, U_L = get_first_intersections(results_pair)
        # I think even if it only returns one value it is still a list
        T_c = T_c[0]
        U_L = U_L[0]
        T_cs[pair] = T_c
        U_Ls[pair] = U_L
        T_1 = results_pair[size1]["T"]
        T_2 = results_pair[size2]["T"]
        U_L_1 = results_pair[size1]["U_L"]
        U_L_2 = results_pair[size2]["U_L"]
        T, U_L_1, U_L_2 = find_common_range(T_1, T_2, U_L_1, U_L_2)
        cum_dic = {}
        cum_dic[size1] = U_L_1
        cum_dic[size2] = U_L_2
        print(pair)
        diff_1 = simple_diff(T_c, T_1, U_L_1)
        diff_2 = simple_diff(T_c, T_2, U_L_2)
        #diff_arr, _ = calc_diff_at(T_c, T, cum_dic)
        # diff_arr should now only have two values I think
        # We get most sizes twice is the thing
        if size1 in dU_dT.keys():
            dU_dT[size1].append(diff_1)
        else:
            dU_dT[size1] = [diff_1]
        if size2 in dU_dT.keys():
            dU_dT[size2].append(diff_2)
        else:
            dU_dT[size2] = [diff_2]

    for size in dU_dT:
        dU_dT[size] = np.mean(dU_dT[size])

    diff_fit_arr = np.array(list(dU_dT.values()))
    size_arr_fit = np.array(list(dU_dT.keys()))

    popt, _ = curve_fit(linear_fit, np.log(size_arr_fit), np.log(diff_fit_arr))
    nu = 1 / popt[0]

    print("FITTING RESULTS:")
    print("nu = ", nu)

    fig, ax = plt.subplots(1, 1)
    ax.plot(size_arr_fit, poly(size_arr_fit, 1 / nu, np.exp(popt[1])), label=rf"$\nu = {nu:.2f}$", color=colors[0])
    ax.plot(size_arr_fit, diff_fit_arr, linestyle="", marker="x", color=colors[0])
    ax.set_xlabel("L")
    ax.set_ylabel(r"$\frac{d U_L}{d \varepsilon}$")
    ax.legend()
    ax.set_title(r"$\frac{d U_L}{d \varepsilon}$ for different System sizes $L$")
    configure_ax(fig, ax)
    # save_plot(root, "/critical_exponent.pdf", format="pdf")
    fig.savefig(simulation_folder + "/critical_exponent_time_avg.png", format="png", dpi=250, transparent=transparent_plots)
    plt.show()

    # okay we will use the results_dic real quick to calculate the intersection in dependence of L
    # We want to use at least three (for now two) sizes to calculate the critical temperature.
    if len(sizes) > 5:
        min_sizes = []
        nus = []
        for size in size_arr_fit:
            min_size = size

            remaining_diffs = diff_fit_arr[size_arr_fit > min_size]
            remaining_sizes = size_arr_fit[size_arr_fit > min_size]
            print(min_size)
            print(remaining_sizes)
            if len(remaining_sizes) > 2:
                min_sizes.append(min_size)

                # them i want to fit
                popt, _ = curve_fit(linear_fit, np.log(remaining_sizes), np.log(remaining_diffs))
                nu = 1 / popt[0]
                nus.append(nu)

        # that seems to be it already
        fig, ax = plt.subplots(1, 1)
        print(min_sizes)
        print(nus)
        ax.errorbar(min_sizes, nus, np.zeros_like(nus), capsize=2, label="", linestyle=None, marker="s")
        ax.set_title(r"Dependence of $\nu$ on the minimum system size")
        ax.set_ylabel(r"$\nu$")
        ax.set_xlabel(r"$L_{min}$")
        configure_ax(fig, ax)
        plt.savefig(simulation_folder + "/nu-L_dependence.png", format="png", dpi=250)
        plt.show()




if __name__ == "__main__":
    main()
