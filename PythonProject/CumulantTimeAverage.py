import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from FunctionsAndClasses import *
from scipy.optimize import curve_fit


def main():
    simulation_folder = '../../Generated content/Silicon/Subsystems/Dampening/Large eta'
    threshold = 25000  # Example threshold value, adjust as needed
    max_L_fit = 100
    transparent_plots = False
    linewidth = 1
    min_T = 0
    max_T = 0
    nr_curves = 6

    results = {}

    for size_folder in os.listdir(simulation_folder):
        if (size_folder[0] != ".") & (size_folder != "plots"):
            size_folder_path = os.path.join(simulation_folder, size_folder)
            if os.path.isdir(size_folder_path):
                size_result = process_size_folder(size_folder_path, threshold)
                results[int(size_folder)] = size_result


    x_range, U_L_intersection, T_intersection, U_L_interpolated = interpolate_and_minimize(results)
    print("Critical Temperature T_c = ", T_intersection)

    fig, ax = plt.subplots(1, 1)

    y_upper_lim = 0
    y_lower_lim = np.infty
    shown_inds = np.linspace(0, nr_curves, nr_curves+1, endpoint=True, dtype=np.int64)
    ind = 0
    for i, size in enumerate(sorted(results.keys())):
        if i in shown_inds:
            T = np.array(results[size]["T"])
            U_L = np.array(results[size]["U_L"])
            print(T)
            print(U_L)
            ax.plot(T, U_L, linestyle="", marker="x", color=colors[ind])
            ax.plot(x_range, U_L_interpolated[i], color=colors[ind], label=rf"L = {size}", linewidth=linewidth)
            ind += 1
            if max_T:
                y_upper_lim = np.maximum(np.max(U_L[(min_T < T) & (T < max_T)]), y_upper_lim)
                y_lower_lim = np.minimum(np.min(U_L[(min_T < T) & (T < max_T)]), y_lower_lim)

    y_span = y_upper_lim - y_lower_lim
    print(y_upper_lim, ", ", y_lower_lim)

    ax.set_xlabel("T")
    ax.set_ylabel(r"$U_L$")
    ax.set_title("Binder Cumulant on T")
    if min_T:
        ax.set_xlim(min_T, ax.get_xlim()[1])
        ax.set_ylim(y_lower_lim - 0.2 * y_span, y_upper_lim + 0.2 * y_span)
    if max_T:
        ax.set_xlim(ax.get_xlim()[0], max_T)
        ax.set_ylim(y_lower_lim - 0.2 * y_span, y_upper_lim + 0.2 * y_span)
    mark_point(ax, T_intersection, U_L_intersection, label = rf"$T_c = {T_intersection:.4f}$")
    configure_ax(fig, ax)
    fig.savefig(simulation_folder + "/cum_time_avg.png", format="png", dpi=300, transparent=transparent_plots)
    plt.show()

    # constructing cum dic
    cum_dic = {}
    for size in results:
        cum_dic[size] = results[size]["U_L"]


    diff_arr, size_arr = calc_diff_at(T_intersection, list(results.values())[0]["T"], cum_dic)
    diff_fit_arr = diff_arr[size_arr < max_L_fit]
    size_arr_fit = size_arr[size_arr < max_L_fit]

    popt, _ = curve_fit(linear_fit, np.log(size_arr_fit), np.log(diff_fit_arr))
    nu = 1 / popt[0]

    print("FITTING RESULTS:")
    print("nu = ", nu)

    fig, ax = plt.subplots(1, 1)
    L_fit = np.linspace(0, np.max(size_arr) + 0.2 * np.max(size_arr), 101)
    ax.plot(L_fit, poly(L_fit, 1 / nu, np.exp(popt[1])), label=rf"$\nu = {nu:.2f}$", color=colors[0])
    ax.plot(size_arr , diff_arr, linestyle="", marker="x", color=colors[0])
    ax.set_xlabel("L")
    ax.set_ylabel(r"$\frac{d U_L}{d \varepsilon}$")
    ax.legend()
    ax.set_title(r"$\frac{d U_L}{d \varepsilon}$ for different System sizes $L$")
    configure_ax(fig, ax)
    # save_plot(root, "/critical_exponent.pdf", format="pdf")
    fig.savefig(simulation_folder + "/critical_exponent_time_avg.png", format="png", dpi=250, transparent=transparent_plots)
    plt.show()

if __name__ == "__main__":
    main()
