import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from FunctionsAndClasses import *
from scipy.optimize import curve_fit

def process_file(file_path, threshold):
    """
    Process a single file and calculate the average after the given threshold.
    """
    df = pd.read_csv(file_path)
    df = df[df['t'] >= threshold]
    nr_values = df.size // 2
    average_value = df['U_L'].mean()
    return average_value, nr_values


def process_size_folder(size_folder, threshold):
    """
    Process all files in a size folder and return a dictionary with temperature and average values.
    """
    result = {'T': [], 'U_L': []}

    for temp_folder in os.listdir(size_folder):
        temp_folder_path = os.path.join(size_folder, temp_folder)
        if os.path.isdir(temp_folder_path):
            temp_average = []
            nr_avg_values = []
            for file_name in os.listdir(temp_folder_path):
                if file_name.endswith('.cum'):
                    file_path = os.path.join(temp_folder_path, file_name)
                    average_value, nr_values = process_file(file_path, threshold)
                    print(file_path)
                    print(average_value)
                    temp_average.append(average_value)
                    nr_avg_values.append(nr_values)
            if temp_average:
                print("mean:", np.mean(temp_average))
                result['T'].append(float(temp_folder))
                U_L = np.sum(np.array(temp_average) * np.array(nr_avg_values)) / np.sum(nr_avg_values)
                result['U_L'].append(U_L)

    result['U_L'] = np.array(result['U_L'])[np.argsort(result['T'])]
    result['T'] = np.sort(result['T'])

    return result




def main():
    simulation_folder = '../../Generated content/Silicon/Subsystems/Time Integral2/'
    threshold = 10000  # Example threshold value, adjust as needed
    max_L_fit = 100
    transparent_plots = False

    results = {}

    for size_folder in os.listdir(simulation_folder):
        size_folder_path = os.path.join(simulation_folder, size_folder)
        if os.path.isdir(size_folder_path):
            size_result = process_size_folder(size_folder_path, threshold)
            results[int(size_folder)] = size_result

    print(results)

    x_range, U_L_intersection, T_intersection, U_L_interpolated = interpolate_and_minimize(results)
    print("Critical Temperature T_c = ", T_intersection)

    fig, ax = plt.subplots(1, 1)

    for i,size in enumerate(sorted(results.keys())):
        ax.plot(results[size]["T"], results[size]["U_L"], linestyle="", marker="x", color=colors[i])
        ax.plot(x_range, U_L_interpolated[i], color=colors[i], label=rf"L = {size}")


    mark_point(ax, T_intersection, U_L_intersection)
    ax.set_xlabel("T")
    ax.set_ylabel(r"$U_L$")
    ax.set_title("Binder Cumulant on T")
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
