from matplotlib import pyplot as plt
from FunctionsAndClasses import *
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.ticker as ticker

def lorentzian(x, a, x0, gamma):
    return a * (gamma**2 / ((x - x0)**2 + gamma**2))

def lorentz_ft(x, xi, a, b):
    return b + a * xi ** 2 / (1 + (x) ** 2 * xi ** 2)

def MF_lorentz(x, xi, a):
    return a * xi / (1 + (x) ** 2 * xi ** 2)

def MF_lorentz_shift(x, xi, a, b):
    return b + a * xi / (1 + (x) ** 2 * xi ** 2)

def linear_fit(x, m, a):
    return a + m*x


def fit_lorentz(p, ft, fitfunc=lorentzian, errors=None):
    try:
        popt, pcov = curve_fit(fitfunc, p, ft, sigma=errors)
        perr = np.sqrt(np.diag(pcov))

    except RuntimeError:
        exit()
        # function has form of a delta peak
        # we delete the largest value
        ind = np.argmax(ft)

        ft = np.insert(ft, ind + 1, 1/2 * ft[ind])
        ft = np.insert(ft, ind, 1 / 2 * ft[ind])
        p = np.insert(p, ind + 1, p[ind] + 1/2 * p[ind + 1])
        p = np.insert(p, ind, (p[ind] + 1 / 2 * p[ind-1]))

        #ft = np.delete(ft, ind)
        #p = np.delete(p, ind)
        print("had to insert values")
        print(p)

        # maybe not cut off the maximum but insert values?

        return fit_lorentz(p, ft)
    return popt, perr

def plot_struct_func(px, py, fx, fy, error_x=np.array([]), error_y=np.array([])):
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    axx = axes[0]
    axy = axes[1]

    axx.errorbar(px, fx, yerr=error_x, ls=" ", marker=".", label="Structure Func", color="C1", ecolor="black", capsize=3)
    axx.set_xlabel(r"$p_x$")
    axx.set_ylabel(r"$S(p_x)$")
    axy.errorbar(py, fy, yerr=error_y, ls=" ", marker=".", label="Structure Func", color="C1", ecolor="black", capsize=3)
    axy.set_xlabel(r"$p_y$")
    axy.set_ylabel(r"$S(p_y)$")
    axx.legend()
    axy.legend()
    return fig, axes

def analyze(df, parameters=None, savepath="./structfact.png", cutoff=np.pi/2, fitfunc=lorentzian, errors_for_fit=True, plot=False):


    px = df["px"]
    px = np.array(px)
    indices = [(-cutoff < x) & (x < cutoff) for x in px]
    # cutoff
    px = px[indices]
    ft_avg_y = np.array(df["ft_avg_y"])[indices]


    py = df["py"]
    py = np.array(py)[indices]
    ft_avg_x = np.array(df["ft_avg_x"])[indices]
    try:
        y_error = np.array(df["stddev_y"])
        x_error = np.array(df["stddev_x"])
    except:
        y_error = None
        x_error = None

    # sorting
    #ft_avg_x = ft_avg_x[np.argsort(px)]
    #ft_avg_x = ft_avg_x[np.argsort(px)]
    #px = np.sort(px)
    #ft_avg_y = ft_avg_y[np.argsort(py)]
    #py = np.sort(py)
    #px = ((px + 2 * max(px)) % (2 * max(px))) - max(px)
    #py = ((py + 2 * max(px)) % (2 * max(py))) - max(py)
    if errors_for_fit:
        popt_x, perr_x = fit_lorentz(px, ft_avg_y, fitfunc, y_error)
        popt_y, perr_y = fit_lorentz(py, ft_avg_x, fitfunc, x_error)
    else:
        popt_x, perr_x = fit_lorentz(px, ft_avg_y, fitfunc, None)
        popt_y, perr_y = fit_lorentz(py, ft_avg_x, fitfunc, None)

    #print("a = %g" % popt_x[0])
    #print("x0 = %g" % popt_x[1])
    #print("gamma = %g" % popt_x[2])


    # xix = 1 / (np.abs(popt_x[2]) * 2)
    # xiy = 1/ (np.abs(popt_y[2]) * 2)
    xix = np.abs(popt_x[0])
    xix_err = perr_x[0]
    xiy = np.abs(popt_y[0])
    xiy_err = perr_y[0]
    xi = np.abs(1 / 2 * (xix + xiy))
    # TODO I think gaussian error propagation just yields for the
    # error of xi:
    xi_err = 1/2 * (xix_err + xiy_err)
    # plotting
    p = np.linspace(min(px), max(px), px.size)
    lorentz_x = fitfunc(p, *popt_x)
    lorentz_y = fitfunc(p, *popt_y)
    if plot:
        fig, axes = plot_struct_func(px, py, ft_avg_y, ft_avg_x, y_error, x_error)
        axes[0].set_title(rf"$\xi_x = {xix:.2f} \quad T = {parameters['tau']:2f}$")
        axes[1].set_title(rf"$\xi_y = {xiy:.2f}\quad T = {parameters['tau']:2f}$")
        axes[0].plot(p, lorentz_x, label="Lorentzian fit")
        axes[1].plot(p, lorentz_y, label="Lorentzian fit")
        plt.tight_layout()
        plt.savefig(savepath, format="png")
        plt.show()
    #print("FWHM x:", np.abs(popt_x[2]) * 2)
    #print("FWHM y:", np.abs(popt_y[2]) * 2)
    #print("Corr Length x:", xix)
    #print("Corr Length y:", xiy)
    return xix, xiy, xix_err, xiy_err, xi, xi_err


def main():
    # parameters
    #root = "../../Generated content/Defense/Quench Small"
    root = "../../Generated content/Defense2/Large Quench/"
    #root = "../../Generated content/Trash/New/Overdamped Quenching 2"

    name = "struct.fact"
    png_name = "struct.fact-fit2"
    root_dirs = os.listdir(root)
    cutoff = np.pi
    fitfunc = MF_lorentz
    errors_for_fit = False
    min_tau = 100
    plot_struct = False
    print(root_dirs)
    config = {
        "ylabelsize": 14,
        "xlabelsize": 14,
        "legendfontsize": 14,

    }
    # arrays to save the xi corrsponding to T
    T_arr = []
    xix_arr = []
    xix_err_arr = []
    xiy_arr = []
    xiy_err_arr = []
    xi_arr = []
    xi_err_arr = []
    tau_arr = []
    Jx_Jy = 0
    # Loop through the directory contents and print the directories
    for item in root_dirs:

        if (item != "plots"):
            # Create the full path to the item
            dir_path = os.path.join(root, item)

            # Check if the item is a directory
            if os.path.isdir(dir_path) & (dir_path != root + "plots"):
                filename = dir_path + "/" + name
                print("reading: ", filename)
                files = os.listdir(dir_path)
                parameters = {}
                for f in files:
                    # we take the first file to be the parameters
                    if(os.path.splitext(f)[1] == ".txt"):
                        parameters = read_parameters_txt(os.path.join(dir_path, f))
                df = read_struct_func(filename)

                xix, xiy, xix_err,\
                    xiy_err, xi, xi_err = analyze(df, parameters,
                                                  savepath=dir_path + png_name,
                                                  cutoff=cutoff, fitfunc=fitfunc,
                                                  errors_for_fit=errors_for_fit,
                                                  plot=plot_struct)

                tau_arr.append(parameters["tau"])
                xix_arr.append(xix)
                xiy_arr.append(xiy)
                xix_err_arr.append(xix_err)
                xiy_err_arr.append(xiy_err)
                xi_arr.append(xi)
                xi_err_arr.append(xi_err)
                Jx_Jy = 1 / (np.sqrt(parameters["J"] / parameters["Jy"]))

    print("Where ticks")

    xix_sorted = np.array(xix_arr)[np.argsort(tau_arr)]
    xiy_sorted = np.array(xiy_arr)[np.argsort(tau_arr)]
    xi_sorted = np.array(xi_arr)[np.argsort(tau_arr)]
    xix_err_sorted = np.array(xix_err_arr)[np.argsort(tau_arr)]
    xiy_err_sorted = np.array(xiy_err_arr)[np.argsort(tau_arr)]
    xi_err_sorted = np.array(xi_err_arr)[np.argsort(tau_arr)]
    tau_arr = np.sort(tau_arr)

    # fitting linear fit
    print(xi_err_sorted)
    xi_fit = xi_sorted[tau_arr > min_tau]
    tau_fit_arr = tau_arr[tau_arr > min_tau]
    if not errors_for_fit:
        xi_err_sorted = None
    popt, _ = curve_fit(linear_fit, np.log(tau_fit_arr), np.log(xi_fit), sigma=xi_err_sorted)
    print("FITTING RESULTS")

    print(popt)
    # plotting
    fig, ax = plt.subplots(1, 1, figsize=(1.5 * 6.4, 1.5 * 4.8))
    # Setze Tickmarken und Labels
    ax.tick_params(direction='in', which='both', length=6, width=2, labelsize=9)
    ax.tick_params(direction='in', which='minor', length=3, width=1, labelsize=9)
    span = np.max(tau_arr) - np.min(tau_arr)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=span / 4))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(base=span / 4 / 5))
    ax.set_yscale("log")
    ax.set_xscale("log")
    # TODO minor locator muss
    #ax.yaxis.set_minor_locator((plt.MultipleLocator(0.2)))
    # Füge Gitterlinien hinzu
    ax.grid(which='major', linestyle='--', alpha=0.5)
    ax.errorbar(tau_arr, xix_sorted, yerr=xix_err_sorted, ls="", marker="x", color="C0", ecolor="black", capsize=3)
    ax.errorbar(tau_arr, xiy_sorted, yerr=xiy_err_sorted, ls="", marker="x", color="C1", ecolor="black", capsize=3)
    ax.plot(tau_arr, poly(tau_arr, popt[0], np.exp(popt[1])), color="C3", label="fit")
    ax.set_xlabel(r"$\tau_Q$")
    ax.set_ylabel(r"$\xi(\tau_Q)$")
    ax.set_title("Corr Length depending on Quench time")
    save_plot(root, "/xix-xiy.png")
    # plotting xi
    fig, ax = plt.subplots(1, 1)
    # Setze Tickmarken und Labels
    ax.tick_params(direction='in', which='both', length=6, width=2, labelsize=9)
    ax.tick_params(direction='in', which='minor', length=3, width=1, labelsize=9)


    ax.set_yscale("log")
    ax.set_xscale("log")
    # TODO minor locator muss
    # ax.yaxis.set_minor_locator((plt.MultipleLocator(0.2)))
    # Füge Gitterlinien hinzu
    ax.grid(which='major', linestyle='--', alpha=0.5)
    ax.errorbar(tau_arr, xi_sorted, yerr=xi_err_sorted, ls="", marker="x", color="C0", ecolor="black", capsize=3)
    tau_arr = np.append(tau_arr, tau_arr[-1] * 1.5)
    print(tau_arr)
    ax.plot(tau_arr, poly(tau_arr, popt[0], np.exp(popt[1])), color=colors[2],
            label=r"$\frac{\nu}{1 + \nu z} =$" + f"{popt[0]:.2f}")
    ax.set_xlabel(r"$\tau_Q$")
    ax.set_ylabel(r"$\xi(\tau_Q)$")
    ax.set_title(r"Corr Length depending on $\tau_Q$")
    ax.legend()
    ax.set_ylim((0.8, ax.get_ylim()[1]))
    configure_ax(fig, ax, config)
    save_plot(root, "/xi.png")
    plt.show()

    # Now we also want to plot the ratio of xi_x / xi_y vs tau_Q and J_x / J_y

    xi_ratio = xix_sorted / xiy_sorted
    print(xi_ratio)
    print(tau_arr)
    fig, ax = plt.subplots(1, 1)
    ax.set_xscale("log")
    ax.set_xlabel(r"$\tau_Q$")
    ax.set_ylabel(r"$\frac{\xi_y}{\xi_x}$")
    ax.set_title("Ratio of correlation lengths vs ratio of coupling constants")
    ax.plot(tau_arr[:-1], xi_ratio, ls="", marker="x", label=r"$\frac{\xi_y}{\xi_x}$")
    ax.plot((tau_arr[0], tau_arr[-2]), (Jx_Jy, Jx_Jy), label=r"$\sqrt{\frac{J_y}{J_x}}$")
    configure_ax(fig, ax)
    plt.show()


if __name__ == "__main__":
    main()