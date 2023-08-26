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
    root = "../../Generated content/Defense2/Binder Small Detailed"
    name = "binder.cumulants"
    name2 = "corr.lengths"
    root_dirs = os.listdir(root)
    print(root_dirs)
    errors_for_fit = False
    # arrays to save the xi corrsponding to T

    T_dic = {}
    cum_dic = {}
    xix_dic = {}
    xiy_dic = {}
    L_xi_dic = {}
    cum_err_dic = {}
    m_dic = {}
    interpol_dic = {}
    interpol_L_xi_dic = {}
    exclude_large_dists = 2
    exclude_small_dists = 0
    min_temp = 0.62
    max_temp = 0.9
    xi_exclude_large_dists = 2
    xi_exclude_small_dists = 0
    r = 2


    cum_path = root + "/" + name
    xi_path = root + "/" + name2
    n = get_size(root, temp_dirs=os.listdir(root))

    df = pd.read_csv(cum_path, delimiter=",", index_col=False)
    df_xi = pd.read_csv(xi_path, delimiter=",", index_col=False)

    labels = df.columns.delete(0)
    labels = labels[2*exclude_large_dists:len(labels)-2*exclude_small_dists]

    xi_labels = df_xi.columns.delete(0)
    xi_labels = xi_labels[2*xi_exclude_large_dists:len(xi_labels)-2*xi_exclude_small_dists]
    T = df["T"]
    for size in labels[::2]:
        cum_dic[size] = np.array(df[size])[np.argsort(T)]
        cum_err_dic[size] = np.array(df[size + "_err"])[np.argsort(T)]

    for size in xi_labels[::2]:
        xix_dic[size] = np.array(df_xi[size])[np.argsort(T)]
        xiy_dic[size] = np.array(df_xi[size + "_y"])[np.argsort(T)]
    T = np.sort(T)
    # sort out temps
    for size in labels[::2]:
        cum_dic[size] = cum_dic[size][(min_temp < T) & (T < max_temp)]
        cum_err_dic[size] = cum_err_dic[size][(min_temp < T) & (T < max_temp)]
        xix_dic[size] = xix_dic[size][(min_temp < T) & (T < max_temp)]
        xiy_dic[size] = xiy_dic[size][(min_temp < T) & (T < max_temp)]
    T = T[(min_temp < T) & (T < max_temp)]

    T_xi_inter, L_xi_dic = plot_intersection(L_xi_dic, T, interpol_L_xi_dic, r, root, xix_dic, xiy_dic)

    plt.show()
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
    print(cum_dic.keys())
    for i,size in enumerate(cum_dic.keys()):
        # only plot every n-th
        interpol_dic[size] = np.interp(T_inter_arr, T, cum_dic[size])
        if i%r == 0:
            # print(size, cum_dic[size])
            # interpolation
            ax.errorbar(T, cum_dic[size], yerr = cum_err_dic[size], ls="",
                        marker="x", color="C" + str(i), ecolor="black", capsize=3)
            ax.plot(T_inter_arr, interpol_dic[size], color="C" + str(i),
                    label=rf"L = {size}", linewidth=0.5)

    # get intersection
    T_inter, U_inter, i_inter = det_intersection(T_inter_arr, interpol_dic)
    ax.scatter(T_inter, U_inter, marker="x", c="black", s=100)
    ax.plot([T_inter, T_inter], [0, U_inter], c="black", ls='--', lw=2, alpha=0.75)
    ax.plot([0, T_inter], [U_inter, U_inter], c="black", ls='--', lw=2, alpha=0.75)
    plt.xticks(list(plt.xticks()[0]) + [T_inter])

    ax.set_xlim(np.min(T) - 0.02 * np.mean(T), np.max(T) + 0.02 * np.mean(T))
    # ax.set_ylim(0.95, 2.5)
    ax.set_xlabel("T")
    ax.set_ylabel(r"$U_L$")
    ax.set_title("Binder Cumulant on T")
    ax.legend()
    save_plot(root, "/cum.pdf", format="pdf")

    # Now we got to make a numerical diff for the reduced temp at Tc
    eps = (T_inter_arr - T_inter) / T_inter
    diff_arr, size_arr = calc_diff_at(T_inter, T, cum_dic)
    xi_num_diff_arr, xi_size_arr = calc_diff_at(T_xi_inter, T, L_xi_dic)
    # fig, ax = plt.subplots(1, 1)
    # ax.plot(size_arr, diff_beta_arr)
    plt.show()



    # fitting
    popt, _ = curve_fit(linear_fit, np.log(size_arr), np.log(diff_arr))
    popt_ising, _ = curve_fit(ising_corr_poly_fit, size_arr, diff_arr, maxfev=1000000)
    popt_poly_corr, _ = curve_fit(crit_poly_fit_corr, size_arr, diff_arr, p0=(1, popt_ising[0], popt_ising[1],
                                                                              popt_ising[2]), maxfev=10000000)
    nu = 1 / popt[0]
    nu_poly_corr = popt_poly_corr[0]
    print("FITTING RESULTS:")
    print("nu = ", nu)
    print("nu_poly_corr = ", popt_poly_corr)
    fig, ax = plt.subplots(1, 1)
    # Setze Tickmarken und Labels
    ax.tick_params(direction='in', which='both', length=6, width=2, labelsize=9)
    ax.tick_params(direction='in', which='minor', length=3, width=1, labelsize=9)




    # TODO minor locator muss
    span = np.max(diff_arr) - np.min(diff_arr)
    print(diff_arr)
    print(span)
    ax.yaxis.set_major_locator((plt.MultipleLocator(span / 4)))
    ax.yaxis.set_minor_locator((plt.MultipleLocator(span / 4 / 5)))
    x_span = np.max(xi_size_arr) - np.min(xi_size_arr)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=np.maximum(int(x_span/5), 0.5)))
    # Füge Gitterlinien hinzu
    ax.grid(which='major', linestyle='--', alpha=0.5)
    ax.plot(size_arr, diff_arr, linestyle="", marker="+")
    L_fit = np.linspace(0, np.max(size_arr) + 0.2 * np.max(size_arr), 101)
    ax.plot(L_fit, poly(L_fit, 1 / nu, np.exp(popt[1])), label=rf"$\nu = {nu:.2f}$")
    ax.plot(L_fit, crit_poly_fit_corr(L_fit, *popt_poly_corr), label=rf"$\nu = {nu_poly_corr:.2f}$")
    ax.plot(L_fit, ising_corr_poly_fit(L_fit, *popt_ising), label=rf"$\omega = {popt_ising[2]:.2f}$", linestyle="-.")
    ax.set_xlabel("L")
    ax.set_ylabel(r"$\frac{d U_L}{d \varepsilon}$")
    ax.legend()
    save_plot(root, "/critical_exponent.pdf", format="pdf")
    plt.show()


    # fitting L/xi

    popt, _ = curve_fit(linear_fit, np.log(xi_size_arr), np.log(xi_num_diff_arr))
    nu = 1 / popt[0]

    fig, ax = plt.subplots(1, 1)
    # Setze Tickmarken und Labels
    ax.tick_params(direction='in', which='both', length=6, width=2, labelsize=9)
    ax.tick_params(direction='in', which='minor', length=3, width=1, labelsize=9)


    x_span = np.max(xi_size_arr) - np.min(xi_size_arr)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=np.maximum(int(x_span/5), 0.5)))
    # TODO minor locator muss
    span = np.max(xi_num_diff_arr) - np.min(xi_num_diff_arr)
    ax.yaxis.set_major_locator((plt.MultipleLocator(span / 4)))
    ax.yaxis.set_minor_locator((plt.MultipleLocator(span / 4 / 5)))
    # Füge Gitterlinien hinzu
    ax.grid(which='major', linestyle='--', alpha=0.5)
    ax.plot(xi_size_arr, xi_num_diff_arr, linestyle="", marker="+")
    L_fit = np.linspace(0, np.max(xi_size_arr) + 0.2 * np.max(xi_size_arr), 101)
    ax.plot(L_fit, poly(L_fit, 1 / nu, np.exp(popt[1])), label=rf"$\nu = {nu:.2f}$")
    ax.set_xlabel("L")
    ax.set_ylabel(r"$\frac{d (L/\xi)}{d \varepsilon}$")
    ax.legend()
    save_plot(root, "/critical_exponent_xi.pdf", format="pdf")
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


def calc_diff_at(T_inter, T, value_dic):
    size_arr = []

    diff_beta_arr = []
    num_diff_arr = []
    for size in value_dic.keys():
        # okay we now think of a better method for numerical differentiation
        # first we look for the data point that is the closest to the intersection
        nearest_T_for_size, index_nearest_T = \
            find_nearest_value_and_index(T, T_inter)
        # now we calculate a central difference with the nearest value
        # being in the center
        num_diff = (value_dic[size][index_nearest_T + 1] -
                    value_dic[size][index_nearest_T - 1]) \
                   / (2 * (T[index_nearest_T + 1] - T[index_nearest_T]))
        num_diff_arr.append(num_diff)
        size_arr.append(int(size))
    num_diff_arr = np.array(num_diff_arr)[np.argsort(size_arr)]
    size_arr = np.sort(size_arr)
    return num_diff_arr, size_arr


def plot_intersection(L_xi_dic, T, interpol_L_xi_dic, r, root, xix_dic, xiy_dic):
    fig, ax = plt.subplots(1, 1)
    # Setze Tickmarken und Labels
    ax.tick_params(direction='in', which='both', length=6, width=2, labelsize=9)
    ax.tick_params(direction='in', which='minor', length=3, width=1, labelsize=9)
    x_span = np.max(T) - np.min(T)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=np.round(x_span / 2, 2)))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(base=np.round(x_span / 2, 2) / 5))
    # TODO minor locator muss
    # Füge Gitterlinien hinzu
    ax.grid(which='major', linestyle='--', alpha=0.5)
    T_inter_arr = np.linspace(np.min(T), np.max(T), 300)
    for i, size in enumerate(xix_dic.keys()):
        # only plot every n-th
        L_xi_dic[size] = int(size) / (0.5 * (xix_dic[size] + xiy_dic[size]))
        interpol_L_xi_dic[size] = np.interp(T_inter_arr, T, L_xi_dic[size])
        print(i)
        if i % r == 0:
            # print(size, cum_dic[size])
            # interpolation
            ax.errorbar(T, L_xi_dic[size], yerr=None, ls="",
                        marker="x", color="C" + str(i), ecolor="black", capsize=3)
            ax.plot(T_inter_arr, interpol_L_xi_dic[size], color="C" + str(i),
                    label=rf"L = {size}", linewidth=0.5)
    # get intersection
    y_span = 0
    for key in L_xi_dic.keys():
        y_span = np.maximum(np.max(L_xi_dic[key]) - np.min(L_xi_dic[key]), y_span)
    ax.yaxis.set_major_locator((plt.MultipleLocator(y_span / 4)))
    ax.yaxis.set_minor_locator((plt.MultipleLocator(y_span / 4 / 5)))
    T_xi_inter, xi_inter, xi_inter_ind = det_intersection(T_inter_arr, interpol_L_xi_dic)
    ax.scatter(T_xi_inter, xi_inter, marker="x", c="black", s=100)
    ax.plot([T_xi_inter, T_xi_inter], [0, xi_inter], c="black", ls='--', lw=2, alpha=0.75)
    ax.plot([0, T_xi_inter], [xi_inter, xi_inter], c="black", ls='--', lw=2, alpha=0.75)
    plt.xticks(list(plt.xticks()[0]) + [T_xi_inter])
    ax.set_xlim(np.min(T) - 0.02 * np.mean(T), np.max(T) + 0.02 * np.mean(T))
    ax.set_xlabel("T")
    ax.set_ylabel(r"$L/\xi$")
    ax.set_title(r" $L/\xi$ on T")
    ax.legend()
    save_plot(root, "/L_xi.pdf", format="pdf")
    return T_xi_inter, L_xi_dic


def get_size(size_path, temp_dirs):
    parameters = {}
    for dir in temp_dirs:
        if (os.path.isdir(os.path.join(size_path, dir))):
            dir_path = os.path.join(size_path, dir)
            file_paths = os.listdir(dir_path)
            for f in file_paths:
                if (os.path.splitext(f)[1] == ".txt"):
                    print(os.path.join(dir_path, f))
                    parameters = read_parameters_txt(os.path.join(dir_path, f))
                    return np.sqrt(parameters["n"])


if __name__ == "__main__":
    main()