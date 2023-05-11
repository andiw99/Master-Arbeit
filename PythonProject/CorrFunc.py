from FunctionsAndClasses import *
import pandas as pd
from matplotlib import pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

def lorentz_wiki(x, gamma, a, b, x0):
    return  a + 1/ np.pi * b * (gamma / ((x - x0) ** 2 + gamma ** 2))


def exp(x, xi):
    return np.exp(-x/xi)

def fit_exp(r, C, func=exp):
    popt, pcov = curve_fit(exp, r, C)
    return popt


def main():
    filepath = "../../Generated content/Quadratic/corr.func.csv"
    df = pd.read_csv(filepath, index_col=None)

    dists = np.arange(0, df.shape[0])
    print(df)
    C_x = df.loc[:, "C_x"] / np.max(df.loc[:, "C_x"])
    C_y = df.loc[:, "C_y"] / np.max(df.loc[:, "C_y"])
    C = 1/2 * (C_x + C_y)
    S = np.fft.fft(C)
    k = np.fft.fftfreq(len(dists), d=dists[1] - dists[0])
    fig, axes = plt.subplots(2, 1)

    popt=fit_exp(dists, C)
    print(popt)
    exp_fit = exp(dists, popt[0])
    ms=1.5
    axes[0].plot(dists, C_x, label="x direction", ls="",    marker=".", ms=ms)
    axes[0].plot(dists, C_y, label="y direction", ls="",    marker=".", ms=ms)
    axes[0].plot(dists, C, label="averaged", ls="",         marker=".", ms=ms)
    axes[0].set_xlabel("r")
    axes[0].set_ylabel("C(r)")
    axes[0].set_title(rf"Direct fit of exponential function   $\xi = {popt[0]:.2f}$")
    axes[0].plot(dists, exp_fit, label="fit")
    axes[1].plot(k, S, label="Structure factor", ls="", marker=".")
    axes[1].set_xlabel("k")
    axes[1].set_ylabel("S(k)")
    axes[0].legend()
    axes[1].legend()
    plt.show()

if __name__ == "__main__":
    main()