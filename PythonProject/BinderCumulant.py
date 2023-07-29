from matplotlib import pyplot as plt
from FunctionsAndClasses import *
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.ticker as ticker


def det_intersection(x, y_dic):
    min_total_dist = np.max(list(y_dic.values())[0])        # starting value for the min
    x_inter = 0
    y_inter = 0
    i_inter = 0
    for (i, x_val) in enumerate(x):
        total_dist = 0
        for key1 in y_dic.keys():
            for key2 in y_dic.keys():
                total_dist += np.abs(y_dic[key1][i] - y_dic[key2][i])
        total_dist /= len(y_dic.keys()) ** 2
        if total_dist < min_total_dist:
            min_total_dist = total_dist
            x_inter = x_val
            y_inter = list(y_dic.values())[0][i]
            i_inter = i
    return x_inter, y_inter, i_inter

def main():
    root = "../../Generated content/Coulomb/Binder4/"
    name = "binder.cumulants"
    root_dirs = os.listdir(root)
    print(root_dirs)
    errors_for_fit = False
    # arrays to save the xi corrsponding to T

    T_dic = {}
    cum_dic = {}
    cum_err_dic = {}
    m_dic = {}
    interpol_dic = {}
    exclude_dists = 2

    # Loop through the directory contents and print the directories
    for size_dir in root_dirs:
        size_path = os.path.join(root, size_dir)
        if os.path.isdir(size_path):
            cum_path = size_path + "/"  + name
            print(size_path)
            n = get_size(size_path, temp_dirs = os.listdir(size_path))

            df = pd.read_csv(cum_path, delimiter=",", index_col=False)
            print(df)
            labels = df.columns.delete(0)[exclude_dists:]
            print(labels)
            T = df["T"]
            for size in labels[::2]:
                print(size)
                cum_dic[size] = np.array(df[size])[np.argsort(T)]
                cum_err_dic[size] = np.array(df[size + "_err"])[np.argsort(T)]
            T = np.sort(T)
            # m = np.array(m)[np.argsort(T)]
            # m_dic[n] = m


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
    T_inter_arr = np.linspace(np.min(T), np.max(T), 300)
    for i,size in enumerate(cum_dic.keys()):
        print(cum_dic[size])
        # interpolation
        interpol_dic[size] = np.interp(T_inter_arr, T, cum_dic[size])
        ax.errorbar(T, cum_dic[size], yerr = cum_err_dic[size], ls="",
                    marker="x", color="C" + str(i), ecolor="black", capsize=3)
        ax.plot(T_inter_arr, interpol_dic[size], color="C" + str(i),
                label=rf"L = {size}", linewidth=0.5)

    # get intersection
    T_inter, U_inter, i_inter = det_intersection(T_inter_arr, interpol_dic)
    print((T_inter, U_inter))
    ax.scatter(T_inter, U_inter, marker="x", c="black", s=100)
    ax.plot([T_inter, T_inter], [0, U_inter], c="black", ls='--', lw=2, alpha=0.75)
    ax.plot([0, T_inter], [U_inter, U_inter], c="black", ls='--', lw=2, alpha=0.75)
    plt.xticks(list(plt.xticks()[0]) + [T_inter])

    ax.set_xlim(np.min(T) - 0.02 * np.mean(T), np.max(T) + 0.02 * np.mean(T))
    ax.set_ylim(0.95, 3.1)
    ax.set_xlabel("T")
    ax.set_ylabel(r"$U_L$")
    ax.set_title("Binder Cumulant on T")
    ax.legend()
    save_plot(root, "/cum.png")

    # Now we got to make a numerical diff at Tc
    size_arr = []
    diff_arr = []
    for size in cum_dic.keys():
        dU_dT = np.gradient(interpol_dic[size], T_inter_arr)
        diff = dU_dT[i_inter]

        size_arr.append(int(size))
        diff_arr.append(diff)
        print(size, diff)
    diff_arr = np.array(diff_arr)[np.argsort(size_arr)]
    size_arr = np.sort(size_arr)

    # fitting
    popt, _ = curve_fit(linear_fit, np.log(size_arr), np.log(diff_arr))
    nu = 1 / popt[0]
    print("FITTING RESULTS:")
    print("nu = ", nu)
    fig, ax = plt.subplots(1, 1)
    # Setze Tickmarken und Labels
    ax.tick_params(direction='in', which='both', length=6, width=2, labelsize=9)
    ax.tick_params(direction='in', which='minor', length=3, width=1, labelsize=9)


    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=2))
    # TODO minor locator muss
    span = int(np.max(diff_arr) - np.min(diff_arr))
    ax.yaxis.set_major_locator((plt.MultipleLocator(span / 4)))
    ax.yaxis.set_minor_locator((plt.MultipleLocator(span / 4 / 5)))
    # Füge Gitterlinien hinzu
    ax.grid(which='major', linestyle='--', alpha=0.5)
    ax.plot(size_arr, diff_arr, linestyle="", marker="+")
    ax.plot(size_arr, poly(size_arr, 1 / nu, np.exp(popt[1])))
    plt.show()


#
    ## Setze Tickmarken und Labels
    #ax.tick_params(direction='in', which='both', length=6, width=2, labelsize=9)
    #ax.tick_params(direction='in', which='minor', length=3, width=1, labelsize=9)
#
#
    #ax.xaxis.set_major_locator(ticker.MultipleLocator(base=0.2))
    #ax.xaxis.set_minor_locator(ticker.MultipleLocator(base=0.04))
    ## TODO minor locator muss
    #ax.yaxis.set_minor_locator((plt.MultipleLocator(0.2)))
    ## Füge Gitterlinien hinzu
    #ax.grid(which='major', linestyle='--', alpha=0.5)
    #for size in T_dic.keys():
    #    ax.plot(T_dic[size], m_dic[size], ls="", marker="+")
    #ax.set_xlabel("T")
    #ax.set_ylabel(r"$\xi(T)$")
    #ax.set_title("Magnetization on T")
    #save_plot(root, "/mag.png")

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