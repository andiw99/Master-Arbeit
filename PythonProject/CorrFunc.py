from FunctionsAndClasses import *
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
import os

#import matplotlib; matplotlib.use("TkAgg")

def lorentz_wiki(x, gamma, a, b, x0):
    return  a + 1/ np.pi * b * (gamma / ((x - x0) ** 2 + gamma ** 2))


def exp(x, xi):
    return np.exp(-x/xi)

def fit_exp(r, C, func=exp):
    popt, pcov = curve_fit(exp, r, C)
    return popt


def corr_scaling_right(T, Tc, nu, xi0):
    eps = (Tc - T) / Tc         # so negative temps are above Tc
    return xi0 / (-eps ** nu)

def corr_scaling_left(T, Tc, nu, xi0):
    eps = (Tc - T) / Tc         # so negative temps are above Tc
    return xi0 / (eps ** nu)

def fit_scaling(T, xi, func=corr_scaling_right):
    popt, pcov = curve_fit(func, T, xi)
    return popt

def corr_eps_scaling(eps, nu, xi0):
    return xi0 / (eps ** nu)


def fit_eps_scaling(eps, xi):
    popt, pcov = curve_fit(corr_eps_scaling, eps,  xi)
    return popt

def log_scaling(T, Tc, nu, xi0):
    return np.log(xi0) + nu * (np.log(Tc) - np.log(Tc-T))

def fit_log_scaling(T, lnxi):
    popt, pcov = curve_fit(log_scaling, T,  lnxi)
    return popt

def main():
    root = "../../Generated content/Coulomb/Detailed/"
    name = "corr.func"
    # get directories of detaileder
    root_dirs = os.listdir(root)
    plot_fits = True
    side_func = corr_scaling_left
    ms = 5

    # define a cutoff since the exponential decay is only correct for large distances
    cutoff_range = [0]
    cutoff_end = -1
    for cutoff in cutoff_range:


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
                C_x = np.array(df.loc[:, "C_x"]) / np.max(df.loc[:, "C_x"])
                C_y = np.array(df.loc[:, "C_y"]) / np.max(df.loc[:, "C_y"])
                C = 1/2 * (C_x + C_y)
                S = np.fft.fft(C)
                k = np.fft.fftfreq(len(dists), d=dists[1] - dists[0])
                print(filepath)
                popt=fit_exp(dists[cutoff:cutoff_end], C[cutoff:cutoff_end])
                xi = popt[0]

                # save
                T_arr.append(T)
                xi_arr.append(xi)
                # plotting
                if plot_fits:
                    fig, axes = plt.subplots(2, 1)
                    exp_fit = exp(dists, popt[0])

                    print(dists)
                    print(C_x)
                    axes[0].set_yscale('log')
                    lower_y_lim = np.minimum(np.min(C_x), np.min(C_y))
                    print("lower y lim: ", lower_y_lim)
                    upper_y_lim = np.maximum(np.max(C_x), np.max(C_y))
                    print("upper y lim: ", upper_y_lim)
                    axes[0].set_ylim(0.00001, 3 )
                    axes[0].plot(dists, C_x, label="x direction", ls="",    marker=".", linewidth=2, ms=ms)
                    axes[0].plot(dists, C_y, label="y direction", ls="",    marker=".", linewidth=2, ms=ms)
                    axes[0].plot(dists, C, label="averaged", ls="",         marker="+", linewidth=2, ms=2*ms)
                    axes[0].set_xlabel("r")
                    axes[0].set_ylabel("C(r)")
                    axes[0].set_title(rf"   $T = {T:.2f} \qquad \xi = {popt[0]:.2f}$")
                    axes[0].plot(dists, exp_fit, label="fit", alpha=0.2)
                    axes[1].plot(k, S, label="Structure factor", ls="", marker=".")
                    axes[1].set_xlabel("k")
                    axes[1].set_ylabel("S(k)")
                    axes[0].legend()
                    axes[1].legend()
                    try:
                        plt.savefig(root + "plots/corr/" + str(T), format="png")
                    except FileNotFoundError:
                        os.makedirs(root + "plots/corr/")
                        plt.savefig(root + "plots/corr/" + str(T), format="png")
                    plt.show()

    # sorting

        detailed_T = np.linspace(0, 200, 201)

        xi_sorted = np.append(np.array(xi_arr)[np.argsort(T_arr)], [])
        T_arr = np.append(np.sort(T_arr), [])


        fig, ax = plt.subplots(1, 1)
        ax.plot(T_arr, xi_sorted, ls="", marker="o")
        ax.set_xlabel("T")
        ax.set_ylabel(r"$\xi(T)$")
        ax.set_title("Corr Length depending on T")
        plt.savefig(root + "plots/corr/xi-corr-fit", format="png")

        # popt = fit_scaling(T_arr, xi_sorted, side_func)
        # fit = side_func(T_arr, *popt)
        # nu = popt[1]
        # Tc = popt[0]
#
        # print("critical exponent nu = ", nu)
        # print("critical temperature Tc = ", Tc)
#
        # T_c_guess = 85
        # eps_arr = (T_c_guess - T_arr) / T_c_guess
        # print(eps_arr)
        # print(xi_sorted)
#
        # popt2 = fit_eps_scaling(eps_arr, xi_sorted)
        # nu2 = popt[1]
        # Tc2 = popt[0]
#
        # eps_fit = corr_eps_scaling(eps_arr, nu2, Tc2)
        # print("critical exponent nu = ", nu2)
        # print("critical temperature Tc = ", Tc2)
#
#
        # popt2 = fit_log_scaling(T_arr, np.log(xi_sorted))
        # nu2 = popt[1]
        # Tc2 = popt[0]
#
        # log_fit = log_scaling(eps_arr, *popt2)
        # print("critical exponent nu = ", nu2)
        # print("critical temperature Tc = ", Tc2)


        plt.show()




if __name__ == "__main__":
    main()