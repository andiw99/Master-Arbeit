from matplotlib import pyplot as plt
from FunctionsAndClasses import *
import numpy as np
from scipy.optimize import curve_fit

def lorentzian(x, a, x0, gamma):
    return a * (gamma**2 / ((x - x0)**2 + gamma**2))

def lorentz_ft(x, xi, a, b):
    return 0 + a * 2 * xi ** 2 / (1 + 4 * np.pi ** 2 *  (x) ** 2 * xi ** 2)
def fit_lorentz(p, ft, fitfunc=lorentzian):
    try:
        popt, pcov = curve_fit(fitfunc, p, ft)

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

def analyze(df, parameters=None, savepath="./structfact.png", cutoff=np.pi/2, fitfunc=lorentzian):

    if not parameters:
        T = 0
    else:
        T = parameters["T"]

    px = df["px"]
    px = np.array(px)
    indices = [(-cutoff < x) & (x < cutoff) for x in px]
    # cutoff
    px = px[indices]
    ft_avg_y = np.array(df.iloc[:, 1])[indices]


    py = df.iloc[:, 2]
    py = np.array(py)[indices]
    ft_avg_x = np.array(df.iloc[:, 3])[indices]

    # sorting
    #ft_avg_x = ft_avg_x[np.argsort(px)]
    #ft_avg_x = ft_avg_x[np.argsort(px)]
    #px = np.sort(px)
    #ft_avg_y = ft_avg_y[np.argsort(py)]
    #py = np.sort(py)
    #px = ((px + 2 * max(px)) % (2 * max(px))) - max(px)
    #py = ((py + 2 * max(px)) % (2 * max(py))) - max(py)

    popt_x = fit_lorentz(px, ft_avg_y, fitfunc)
    popt_y = fit_lorentz(py, ft_avg_x, fitfunc)
    #print("a = %g" % popt_x[0])
    #print("x0 = %g" % popt_x[1])
    #print("gamma = %g" % popt_x[2])


    # xix = 1 / (np.abs(popt_x[2]) * 2)
    # xiy = 1/ (np.abs(popt_y[2]) * 2)
    xix = np.abs(popt_x[0])
    xiy = np.abs(popt_y[0])
    xi = np.abs(1 / 2 * (xix + xiy))

    # plotting
    fig, axes = plot_struct_func(px, py,ft_avg_y, ft_avg_x)
    p = np.linspace(min(px), max(px), px.size)
    lorentz_x = fitfunc(p, popt_x[0], popt_x[1], popt_x[2])
    lorentz_y = fitfunc(p, popt_y[0], popt_y[1], popt_y[2])
    axes[0].plot(p, lorentz_x, label="Lorentzian fit")
    axes[1].plot(p, lorentz_y, label="Lorentzian fit")
    axes[0].set_title(rf"$\xi_x = {xix:.2f} \quad T = {T:2f}$")
    axes[1].set_title(rf"$\xi_y = {xiy:.2f}\quad T = {T:2f}$")
    plt.savefig(savepath, format="png")
    #print("FWHM x:", np.abs(popt_x[2]) * 2)
    #print("FWHM y:", np.abs(popt_y[2]) * 2)
    #print("Corr Length x:", xix)
    #print("Corr Length y:", xiy)
    return xi, T


def main():
    root = "../../Generated content/Coulomb/Random Init Test"
    name = "struct.fact"
    png_name = "struct.fact-fit2"
    root_dirs = os.listdir(root)

    # arrays to save the xi corrsponding to T

    T_arr = []
    xi_arr = []
    cutoff = 3/8 * np.pi
    fitfunc = lorentz_ft
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
                xi, T = analyze(df, parameters, savepath=dir_path + png_name, cutoff=cutoff, fitfunc=fitfunc)
                T_arr.append(T)
                xi_arr.append(xi)


    xi_sorted = np.array(xi_arr)[np.argsort(T_arr)]
    T_arr = np.sort(T_arr)
    fig, ax = plt.subplots(1, 1)
    ax.plot(T_arr, xi_sorted, ls="", marker="o")
    ax.set_xlabel("T")
    ax.set_ylabel(r"$\xi(T)$")
    ax.set_title("Corr Length depending on T")
    save_plot(root, "/xi.png")
    plt.show()


if __name__ == "__main__":
    main()