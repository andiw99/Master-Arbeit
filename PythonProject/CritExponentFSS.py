from matplotlib import pyplot as plt
from FunctionsAndClasses import *
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.ticker as ticker

def lorentzian(x, a, x0, gamma):
    return a * (gamma**2 / ((x - x0)**2 + gamma**2))

def lorentz_ft(x, xi, a, b):
    return b + a * xi ** 2 / (1 +  (x) ** 2 * xi ** 2)

def MF_lorentz(x, xi, a):
    return a * xi / (1 + (x) ** 2 * xi ** 2)
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
    axx.legend()
    axy.legend()
    return fig, axes

def analyze(df, parameters=None, savepath="./structfact.png",
            cutoff=np.pi/2, fitfunc=lorentzian, errors_for_fit=False):

    if not parameters:
        T = 0
    else:
        T = parameters["T"]

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
    xiy = np.abs(popt_y[0])
    xi = np.abs(1 / 2 * (xix + xiy))

    return xi, T


def main():
    root = "../../Generated content/New/Coulomb/Critical Exponent"
    name = "struct.fact"
    png_name = "struct.fact-fit2"
    root_dirs = os.listdir(root)
    print(root_dirs)
    errors_for_fit = False
    # arrays to save the xi corrsponding to T

    T_dic = {}
    xi_dic = {}
    L_xi_dic = {}
    cutoff =  np.pi
    fitfunc = MF_lorentz
    # Loop through the directory contents and print the directories
    for size_dir in root_dirs:
        size_path = os.path.join(root, size_dir)
        if(os.path.isdir(size_path)):
            temp_dirs = os.listdir(size_path)
            # need to find out the size
            # os.listdir(tempdirs[0]) gives files like 0.csv
            n = get_size(size_path, temp_dirs)

            # create Arrays for xi and T of this size directory tree
            T_dic[n] = []
            xi_dic[n] = []
            L_xi_dic[n] = []
            for temp_dir in temp_dirs:

                if (temp_dir != "plots"):
                    # Create the full path to the item
                    dir_path = os.path.join(size_path, temp_dir)

                    # Check if the item is a directory
                    if os.path.isdir(dir_path) & (dir_path != root + "plots"):
                        print("Why you not doing anything")
                        filename = dir_path + "/" + name
                        print(filename)
                        files = os.listdir(dir_path)
                        parameters = {}
                        for f in files:
                            # we take the first file to be the parameters
                            if(os.path.splitext(f)[1] == ".txt"):
                                parameters = read_parameters_txt(os.path.join(dir_path, f))
                        df = read_struct_func(filename)

                        xi, T = analyze(df, parameters,
                                        savepath=dir_path + png_name,
                                        cutoff=cutoff, fitfunc=fitfunc,
                                        errors_for_fit=errors_for_fit)
                        # append the values to the arrays for size n
                        T_dic[n].append(T)
                        xi_dic[n].append(xi)

            # we need to do the sorting for size n
            xi_dic[n] = np.array(xi_dic[n])[np.argsort(T_dic[n])]
            T_dic[n] = np.sort(T_dic[n])
            L_xi_dic[n] = n / xi_dic[n]


    fig, ax = plt.subplots(1, 1)

    # Setze Tickmarken und Labels
    ax.tick_params(direction='in', which='both', length=6, width=2, labelsize=9)
    ax.tick_params(direction='in', which='minor', length=3, width=1, labelsize=9)


    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=0.2))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(base=0.04))
    # TODO minor locator muss
    ax.yaxis.set_minor_locator((plt.MultipleLocator(0.2)))
    # Füge Gitterlinien hinzu
    ax.grid(which='major', linestyle='--', alpha=0.5)
    for size in T_dic.keys():
        ax.plot(T_dic[size], xi_dic[size], ls="", marker="+")
    ax.set_xlabel("T")
    ax.set_ylabel(r"$\xi(T)$")
    ax.set_title("Corr Length depending on T")
    save_plot(root, "/xi.png")

    fig, ax = plt.subplots(1, 1)

    # Setze Tickmarken und Labels
    ax.tick_params(direction='in', which='both', length=6, width=2, labelsize=9)
    ax.tick_params(direction='in', which='minor', length=3, width=1, labelsize=9)


    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=0.2))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(base=0.04))
    # TODO minor locator muss
    ax.yaxis.set_major_locator(ticker.MultipleLocator(base=50))
    ax.yaxis.set_minor_locator((plt.MultipleLocator(10)))
    # Füge Gitterlinien hinzu
    ax.grid(which='major', linestyle='--', alpha=0.5)
    for size in T_dic.keys():
        ax.plot(T_dic[size], L_xi_dic[size], ls="", marker="+")
    ax.set_xlabel("T")
    ax.set_ylabel(r"$\frac{L}{\xi(T, L^{-1})}$")
    ax.set_title("Corr Length depending on T")
    save_plot(root, "/xi.png")


    plt.show()


def get_size(size_path, temp_dirs):
    parameters = {}
    for dir in temp_dirs:
        if (os.path.isdir(os.path.join(size_path, dir))):
            dir_path = os.path.join(size_path, dir)
            file_paths = os.listdir(dir_path)
            for f in file_paths:
                if (os.path.splitext(f)[1] == ".txt"):
                    parameters = read_parameters_txt(os.path.join(dir_path, f))
    return np.sqrt(parameters["n"])


if __name__ == "__main__":
    main()