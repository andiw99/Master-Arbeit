from matplotlib import pyplot as plt
from FunctionsAndClasses import *
import numpy as np
from scipy.optimize import curve_fit

def lorentzian(x, a, x0, gamma):
    return a * (gamma**2 / ((x - x0)**2 + gamma**2))


def fit_lorentz(p, ft):
    try:
        popt, pcov = curve_fit(lorentzian, p, ft)
    except RuntimeError:
        # function has form of a delta peak
        # we delete the largest value
        ind = np.argmax(ft)
        ft = np.delete(ft, ind)
        p = np.delete(p, ind)
        return fit_lorentz(p, ft)
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

def analyze(df, parameters=None):

    if not parameters:
        T = 0
    else:
        T = parameters["T"]

    px = df["px"]
    ft_avg_y = df.iloc[:, 1]
    py = df.iloc[:, 2]
    ft_avg_x = df.iloc[:, 3]


    px = np.array(px)
    px = ((px + 2 * max(px)) % (2 * max(px))) - max(px)
    py = np.array(py)
    py = ((py + 2 * max(px)) % (2 * max(py))) - max(py)

    popt_x = fit_lorentz(px, ft_avg_y)
    popt_y = fit_lorentz(py, ft_avg_x)
    print("a = %g" % popt_x[0])
    print("x0 = %g" % popt_x[1])
    print("gamma = %g" % popt_x[2])


    xix = 1 / (np.abs(popt_x[2]) * 2)
    xiy = 1/ (np.abs(popt_y[2]) * 2)
    xi = 1 / 2 * (xix + xiy)

    # plotting
    #fig, axes = plot_struct_func(px, py,ft_avg_y, ft_avg_x)
    #p = np.linspace(min(px), max(px), px.size)
    #lorentz_x = lorentzian(p, popt_x[0], popt_x[1], popt_x[2])
    #lorentz_y = lorentzian(p, popt_y[0], popt_y[1], popt_y[2])
#
    #axes[0].plot(p, lorentz_x, label="Lorentzian fit")
    #axes[1].plot(p, lorentz_y, label="Lorentzian fit")
#
    #print("FWHM x:", np.abs(popt_x[2]) * 2)
    #print("FWHM y:", np.abs(popt_y[2]) * 2)
    #axes[0].set_title(rf"$\xi_x = {xix:.2f} \quad T = {T:2f}$")
    #axes[1].set_title(rf"$\xi_y = {xiy:.2f}\quad T = {T:2f}$")
    #print("Corr Length x:", xix)
    #print("Corr Length y:", xiy)
    return xi, T


def main():
    root = "../../Generated content/New Scan/"
    name = "struct.fact"
    root_dirs = os.listdir(root)

    # arrays to save the xi corrsponding to T

    T_arr = []
    xi_arr = []

    # Loop through the directory contents and print the directories
    for item in root_dirs:
        # Create the full path to the item
        dir_path = os.path.join(root, item)

        # Check if the item is a directory
        if os.path.isdir(dir_path) & (dir_path != root + "plots"):
            filename = dir_path + "/" + name
            files = os.listdir(dir_path)
            parameters = {}
            for f in files:
                # we take the first file to be the parameters
                if(os.path.splitext(f)[1] == ".txt"):
                    parameters = read_parameters_txt(os.path.join(dir_path, f))
            print("reading: ", filename)
            df = read_struct_func(filename)
            xi, T = analyze(df, parameters)
            T_arr.append(T)
            xi_arr.append(xi)


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