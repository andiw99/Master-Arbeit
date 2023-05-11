from matplotlib import pyplot as plt
from FunctionsAndClasses import *
import numpy as np
from scipy.optimize import curve_fit

def lorentzian(x, a, x0, gamma):
    return a * (gamma**2 / ((x - x0)**2 + gamma**2))


def fit_lorentz(p, ft):
    popt, pcov = curve_fit(lorentzian, p, ft)
    return popt

def plot_struct_func(px, py, fx, fy):
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    axx = axes[0]
    axy = axes[1]


    axx.plot(px, fx, ls=" ", marker=".", label="Structure Func")
    axx.set_xlabel(r"$p_x$")
    axx.set_ylabel(r"$S(p_x)$")
    axy.plot(py, fy, ls=" ", marker=".", label="Structure Func")
    axy.set_xlabel(r"$p_y$")
    axy.set_ylabel(r"$S(p_y)$")
    plt.tight_layout()
    axx.legend()
    axy.legend()
    return fig, axes

def analyze(df):

    px = df["px"]
    ft_avg_y = df.iloc[:, 1]
    py = df.iloc[:, 2]
    ft_avg_x = df.iloc[:, 3]
    print(py)
    print(ft_avg_x)


    px = np.array(px)
    px = ((px + 2 * max(px)) % (2 * max(px))) - max(px)
    py = np.array(py)
    py = ((py + 2 * max(px)) % (2 * max(py))) - max(py)

    fig, axes = plot_struct_func(px, py,ft_avg_y, ft_avg_x)


    popt_x = fit_lorentz(px, ft_avg_y)
    popt_y = fit_lorentz(py, ft_avg_x)
    print("a = %g" % popt_x[0])
    print("x0 = %g" % popt_x[1])
    print("gamma = %g" % popt_x[2])

    p = np.linspace(min(px), max(px), px.size)
    lorentz_x = lorentzian(p, popt_x[0], popt_x[1], popt_x[2])
    lorentz_y = lorentzian(p, popt_y[0], popt_y[1], popt_y[2])

    axes[0].plot(p, lorentz_x, label="Lorentzian fit")
    axes[1].plot(p, lorentz_y, label="Lorentzian fit")

    xix = 1 / (np.abs(popt_x[2]) * 2)
    xiy = 1/ (np.abs(popt_y[2]) * 2)
    print("FWHM x:", np.abs(popt_x[2]) * 2)
    print("FWHM y:", np.abs(popt_y[2]) * 2)
    axes[0].set_title(rf"$\xi = {xix:.2f}$")
    axes[1].set_title(rf"$\xi = {xiy:.2f}$")
    print("Corr Length x:", xix)
    print("Corr Length y:", xiy)


def main():
    filepath = "../../Generated content/Repeat Cooling/structfunc"

    df = read_struct_func(filepath)
    analyze(df)
    plt.show()


if __name__ == "__main__":
    main()