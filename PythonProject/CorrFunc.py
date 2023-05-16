from FunctionsAndClasses import *
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import os

def lorentz_wiki(x, gamma, a, b, x0):
    return  a + 1/ np.pi * b * (gamma / ((x - x0) ** 2 + gamma ** 2))


def exp(x, xi):
    return np.exp(-x/xi)

def fit_exp(r, C, func=exp):
    popt, pcov = curve_fit(exp, r, C)
    return popt


def main():
    root = "../../Generated content/New Scan/"
    name = "corr.func"
    # get directories of detaileder
    root_dirs = os.listdir(root)

    # define a cutoff since the exponential decay is only correct for large distances
    cutoff = 10

    # arrays to save the xi corrsponding to T

    T_arr = []
    xi_arr = []

    # Loop through the directory contents and print the directories
    for item in root_dirs:
        # Create the full path to the item
        dir_path = os.path.join(root, item)

        # Check if the item is a directory
        if os.path.isdir(dir_path) & (dir_path != root + "plots"):
            # we need to read in the parmeters
            files = os.listdir(dir_path)
            parameters = {"T": "not avail"}
            for file in files:
                if os.path.splitext(file)[1] == ".txt":
                    # we assume that every file has the same parameters but they are different realizations
                    parameters = read_parameters_txt(os.path.join(dir_path, file))
                    break
            T = parameters["T"]

            filepath = os.path.join(dir_path, name)
            df = pd.read_csv(filepath, index_col=None)

            dists = np.arange(0, df.shape[0])
            C_x = df.loc[:, "C_x"] / np.max(df.loc[:, "C_x"])
            C_y = df.loc[:, "C_y"] / np.max(df.loc[:, "C_y"])
            C = 1/2 * (C_x + C_y)
            S = np.fft.fft(C)
            k = np.fft.fftfreq(len(dists), d=dists[1] - dists[0])

            popt=fit_exp(dists[cutoff:], C[cutoff:])
            xi = popt[0]

            # save
            T_arr.append(T)
            xi_arr.append(xi)
            # plotting

            #fig, axes = plt.subplots(2, 1)
            #exp_fit = exp(dists, popt[0])
            #ms=1.5
            #axes[0].plot(dists, C_x, label="x direction", ls="",    marker=".", ms=ms)
            #axes[0].plot(dists, C_y, label="y direction", ls="",    marker=".", ms=ms)
            #axes[0].plot(dists, C, label="averaged", ls="",         marker=".", ms=ms)
            #axes[0].set_xlabel("r")
            #axes[0].set_ylabel("C(r)")
            #axes[0].set_title(rf"   $T = {T:.2f} \qquad \xi = {popt[0]:.2f}$")
            #axes[0].plot(dists, exp_fit, label="fit")
            #axes[1].plot(k, S, label="Structure factor", ls="", marker=".")
            #axes[1].set_xlabel("k")
            #axes[1].set_ylabel("S(k)")
            #axes[0].legend()
            #axes[1].legend()
            #try:
            #    plt.savefig(root + "plots/corr/" + str(T), format="png")
            #except FileNotFoundError:
            #    os.makedirs(root + "plots/corr/")
            #    plt.savefig(root + "plots/corr/" + str(T), format="png")
            #plt.show()

# sorting

    xi_sorted = np.array(xi_arr)[np.argsort(T_arr)]
    T_arr = np.sort(T_arr)

    fig, ax = plt.subplots(1, 1)
    ax.plot(T_arr, xi_sorted, ls="", marker="o")
    ax.set_xlabel("T")
    ax.set_ylabel(r"$\xi(T)$")
    ax.set_title("Corr Length depending on T")
    plt.show()




if __name__ == "__main__":
    main()