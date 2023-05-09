import numpy as np
from matplotlib import pyplot as plt
from scipy.optimize import curve_fit
from FunctionsAndClasses import *
import pandas as pd

def lorentzian(x, a, b, x0, gamma):
    return a + b * (gamma**2 / ((x - x0)**2 + gamma**2))


def lorentz_wiki(x, gamma, a, b, x0):
    return  a + 1/ np.pi * b * (gamma / ((x - x0) ** 2 + gamma ** 2))


def lorentz_ft(x, xi, a, b):
    return b + a * 2 * xi / (1 + 4 * np.pi ** 2 *  (x) ** 2 * xi ** 2)

def lorentz_fta(x, a):
    return 2 * a / (a ** 2 + 4 * np.pi ** 2 * x ** 2)


def fit_lorentz(p, ft, lorentz):
    popt, pcov = curve_fit(lorentz, p, ft)
    return popt

def exp(x, xi):
    return np.exp(-x/xi)

def gaussian(x, a, x0, sigma):
    # Define a Gaussian function
    return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

def fwhm(x, y):
    # Fit a Gaussian distribution to the data and extract the FWHM
    p0 = [np.max(y), x[np.argmax(y)], np.std(y)]
    popt, _ = curve_fit(gaussian, x, y, p0=p0)
    return 2 * np.sqrt(2 * np.log(2)) * np.abs(popt[2])


def plot_decay(df, fig, ax):

    x = df.loc[:, "x"]
    f = df.loc[:, "f"]

    ax.plot(x, f, label="exponential decay")
    ax.set_title(label=r"$\xi = 1$")
    ax.set_ylabel(r"$e^{ - r / \xi}$")
    ax.set_xlabel("r")


def plot_lorentz(df, fig, ax):
    p = np.array(df.loc[:, "p"])
    ft = np.array(df.loc[:, "ft"])

    ft = ft[np.argsort(p)]
    p = sorted(p)


    ax.plot(p, ft, label="fourier trafo")
    ax.set_title(label="")
    ax.set_ylabel(r"$S(p)$")
    ax.set_xlabel("p")

def make_fit(df, fig, axes, lorentz=lorentz_wiki):
    x = np.array(df.loc[:, "x"])
    f = np.array(df.loc[:, "f"])
    p = np.array(df.loc[:, "p"])
    ft = np.array(df.loc[:, "ft"])

    ft = ft[np.argsort(p)]
    p = np.sort(p)

    popt = fit_lorentz(p, ft, lorentz)
    print(popt)
    ft_fit = lorentz(p, *popt)

    axes[1].plot(p, ft_fit, label="fit", ls="dashed")

    N = x.size
    xmax = np.max(x)

    xi = 1 / (2* popt[0])
    xi = xi / (N/(2 * xmax))
    #xi = 1/popt[0]
    actual_xi = 1

    f_fit = exp(x , xi)
    f_act = exp(x, actual_xi)
    axes[0].plot(x, f_fit, label=rf"fit with $\xi = {xi:.2f}$", ls = "dashed")


def make_plot(df):

    fig, axes = plt.subplots(2, 1)


    df.loc[:, "ft"] /= df.shape[0]
    plot_decay(df, fig, axes[0])
    plot_lorentz(df, fig, axes[1])


    make_fit(df, fig, axes)

    axes[0].legend()
    axes[1].legend()

    plt.tight_layout()
    plt.show()


def main():

    root = "../../Generated content/LatticeTrafo/1D/"
    name = "1000+20+1.2"
    filename = root + name
    df = pd.read_csv(filename, index_col = None, header=0)
    make_plot(df)




    exit()
    ns = [11, 100, 1000]
    end_x = 20
    for n in ns:
        x = np.linspace(0, end_x, n)
        xi = 1

        y = exp(x, xi)

        y_fft = np.fft.fft(y)
        k = np.fft.fftfreq(len(x), d=x[1] - x[0])

        y_fft = np.abs(y_fft[np.argsort(k)])
        k = sorted(k)

        fig, axes = plt.subplots(2, 1, figsize=(6, 7))
        axes[0].set_ylabel("C(r)")
        axes[0].set_xlabel("r")
        axes[0].plot(x, y, label="Function")
        axes[1].plot(k, y_fft, label="FFT Function")
        axes[1].set_ylabel("S(q)")
        axes[1].set_xlabel("q")
        axes[0].set_title("Correlation Function with Manhattan Metrik")

        width = fwhm(k, y_fft)
        xi_estimate = 1/width
        axes[0].plot(x, exp(x, xi_estimate), label="fit")
        axes[0].legend()

        plt.show()

        print("FWHM: ", width)

if __name__ == "__main__":
    main()