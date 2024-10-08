from matplotlib import pyplot as plt
from FunctionsAndClasses import *
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.ticker as ticker


def main():
    root = "../../Generated content/Silicon/Subsystems/Suite/Test7"
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
    cum_err_dic = {}
    m_dic = {}
    interpol_dic = {}
    exclude_large_dists = 100
    exclude_small_dists = 10
    min_temp = 0
    max_temp = 12
    xi_exclude_large_dists = 0
    xi_exclude_small_dists = 0
    max_L_fit = 1000

    r = 5

    figsize = (1.2 *  6.4, 4.8)
    L_max_lower = 20
    L_max_upper = 50
    transparent_plots = False

    cum_path = root + "/" + name

    df = pd.read_csv(cum_path, delimiter=",", index_col=False)
    print(df)
    labels = df.columns.delete(0)
    labels = labels[::2]
    labels = np.array([int(label) for label in labels])
    labels = labels[(labels > exclude_small_dists) & (labels < exclude_large_dists)]
    labels = [str(label) for label in labels] # lol
    # labels = labels[2*exclude_large_dists:len(labels)-2*exclude_small_dists]
    T = df["T"]
    for size in labels:
        cum_dic[size] = np.array(df[size])[np.argsort(T)]
        cum_err_dic[size] = np.array(df[size + "_err"])[np.argsort(T)]
    T = np.sort(T)

    # sort out temps
    for size in labels:
        cum_dic[size] = cum_dic[size][(min_temp < T) & (T < max_temp)]
        cum_err_dic[size] = cum_err_dic[size][(min_temp < T) & (T < max_temp)]
    T = T[(min_temp < T) & (T < max_temp)]
    try:
        xi_path = root + "/" + name2
        df_xi = pd.read_csv(xi_path, delimiter=",", index_col=False)

        T_xi = df_xi["T"]

        xi_labels = df_xi.columns.delete(0)
        xi_labels = xi_labels[2*xi_exclude_large_dists:len(xi_labels)-2*xi_exclude_small_dists]
        for size in xi_labels[::2]:
            xix_dic[size] = np.array(df_xi[size])[np.argsort(T_xi)]
        for size in xi_labels[1::2]:
            # remove the _y thingy
            xiy_dic[size[:-2]] = np.array(df_xi[size])[np.argsort(T_xi)]
        T_xi = np.sort(T_xi)
        for size in xi_labels[::2]:
            xix_dic[size] = xix_dic[size][(min_temp < T_xi) & (T_xi < max_temp)]
        for size in xi_labels[1::2]:
            xiy_dic[size[:-2]] = xiy_dic[size[:-2]][(min_temp < T_xi) & (T_xi < max_temp)]
        T_xi = T_xi[(min_temp < T_xi) & (T_xi < max_temp)]

        #T_xi_inter, L_xi_dic = plot_intersection(T, r, root, xix_dic, xiy_dic)

        # Okay we want to switch to deal with xi_x and xi_y seperately, so we need to interpolate xix for every size
        # and then use our det_intersection method


        #plt.show()

        # make xix to L_xix
        L_xix_dic = {}
        L_xiy_dic = {}
        for size in xix_dic.keys():
            L_xix_dic[size] = int(size) / xix_dic[size]
        for size in xiy_dic.keys():
            L_xiy_dic[size] = int(size) / xiy_dic[size]

        # plot for x direction
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax.set_title(r"$L/\xi_x$ on $\varepsilon$")
        ax.set_xlabel(r"$\varepsilon$")
        ax.set_ylabel(r"$L/\xi_x$")
        T_x_intersec, L_xix_intersec = plt_inter(ax, fig, T, L_xix_dic, 300, r)
        configure_ax(fig, ax)
        fig.savefig(root + "/L_xix.png", format="png", dpi=300, transparent=transparent_plots)
        plt.show()

        # plot for y direction
        fig, ax = plt.subplots(1, 1, figsize=figsize)
        ax.set_title(r"$L/\xi_y$ on $\varepsilon$")
        ax.set_xlabel(r"$\varepsilon$")
        ax.set_ylabel(r"$L/\xi_y$")
        T_y_intersec, L_xiy_intersec = plt_inter(ax, fig, T, L_xiy_dic, 300, r)
        configure_ax(fig, ax)
        fig.savefig(root + "/L_xiy.png", format="png", dpi=300, transparent=transparent_plots)
        plt.show()

        xix_num_diff_arr, xix_size_arr = calc_diff_at(T_x_intersec, T, L_xix_dic)
        xiy_num_diff_arr, xiy_size_arr = calc_diff_at(T_y_intersec, T, L_xiy_dic)

        # xi_num_diff_arr_fit = xi_num_diff_arr[xi_size_arr < max_L_fit]
        xix_num_diff_arr = xix_num_diff_arr[xix_size_arr < max_L_fit]
        xiy_num_diff_arr = xiy_num_diff_arr[xiy_size_arr < max_L_fit]
        # xi_size_fit_arr = xi_size_arr[xi_size_arr < max_L_fit]
        xix_size_arr = xix_size_arr[xix_size_arr < max_L_fit]
        xiy_size_arr = xiy_size_arr[xiy_size_arr < max_L_fit]
        # popt, _ = curve_fit(linear_fit, np.log(xi_size_fit_arr), np.log(xi_num_diff_arr_fit))

        try:
            popt_x, _ = curve_fit(linear_fit, np.log(xix_size_arr), np.log(xix_num_diff_arr))
            nu_x = 1 / popt_x[0]

            popt_y, _ = curve_fit(linear_fit, np.log(xiy_size_arr), np.log(xiy_num_diff_arr))
            nu_y = 1 / popt_y[0]

            # plotting derivatives for both directions

            fig, ax = plt.subplots(1, 1, figsize=figsize)

            ax.plot(xix_size_arr, xix_num_diff_arr, linestyle="", marker="+", color=colors[0])
            ax.plot(xiy_size_arr, xiy_num_diff_arr, linestyle="", marker="x", color=colors[4])
            L_x_fit = np.linspace(0, np.max(xix_size_arr) + 0.2 * np.max(xix_size_arr), num=101)
            L_y_fit = np.linspace(0, np.max(xiy_size_arr) + 0.2 * np.max(xix_size_arr), num=101)
            ax.plot(L_x_fit, poly(L_x_fit, 1 / nu_x, np.exp(popt_x[1])), label=rf"$\nu_x = {nu_x:.2f}$", color=colors[0])
            ax.plot(L_y_fit, poly(L_y_fit, 1 / nu_y, np.exp(popt_y[1])), label=rf"$\nu_y = {nu_y:.2f}$", color=colors[4])

            ax.set_xlabel("L")
            ax.set_ylabel(r"$\frac{d (L/\xi)}{d \varepsilon}$")
            ax.set_title(r"$\frac{d (L/\xi)}{d \varepsilon}$ for different System sizes $L$ in x- and y-direction")
            configure_ax(fig, ax)
            fig.savefig(root + "/critical_exponent_xi.png", format="png", dpi=250, transparent=transparent_plots)
            plt.show()
        except:
            print("Fitting didnt work")
    except Exception as error:
        print("Corr Length didnt work")
        print(error)


    fig, ax = plt.subplots(1, 1, figsize = figsize)
    T_inter_arr = np.linspace(np.min(T), np.max(T), 300)
    line_nr = 0

    shown_inds = [int(i) for i in np.linspace(0, len(cum_dic.keys()) - 1, r)]
    for i,size in enumerate(cum_dic.keys()):
        # only plot every n-th
        interpol_dic[size] = np.interp(T_inter_arr, T, cum_dic[size])
        if i in shown_inds:
            # print(size, cum_dic[size])
            # interpolation
            print("Cumulant:", size)
            print(T)
            print(cum_dic[size])
            ax.errorbar(T, cum_dic[size], yerr = cum_err_dic[size], ls="",
                        marker="x", color=colors[line_nr % len(colors)], ecolor="black", elinewidth=0, capsize=0)
            ax.plot(T_inter_arr, interpol_dic[size], color=colors[(line_nr) % len(colors)],
                    label=rf"L = {size}", linewidth=0.5)
            line_nr += 1

    # get intersection
    T_inter, U_inter, i_inter = det_intersection(T_inter_arr, interpol_dic)
    print("Critical Temperature T_c = ", T_inter)
    mark_point(ax, T_inter, U_inter, label=rf"$T_c = $ {T_inter:.3f}")
    ax.set_xlabel("T")
    ax.set_ylabel(r"$U_L$")
    ax.set_title("Binder Cumulant on T")
    configure_ax(fig, ax)
    fig.savefig(root + "/cum.png", format="png", dpi=300, transparent=transparent_plots)
    #save_plot(root, "/cum.pdf", format="pdf")

    # Now we got to make a numerical diff for the reduced temp at Tc
    eps = (T_inter_arr - T_inter) / T_inter
    diff_arr, size_arr = calc_diff_at(T_inter, T, cum_dic)

    # xi_num_diff_arr, xi_size_arr = calc_diff_at(T_xi_inter, T, L_xi_dic)
    # numercial diffs for x and y direction aswell as the size arrays


    # fig, ax = plt.subplots(1, 1)
    # ax.plot(size_arr, diff_beta_arr)
    plt.show()


    diff_fit_arr = diff_arr[size_arr < max_L_fit]
    size_arr_fit = size_arr[size_arr < max_L_fit]
    # fitting
    popt, _ = curve_fit(linear_fit, np.log(size_arr_fit), np.log(diff_fit_arr))
    popt_ising, _ = curve_fit(ising_corr_poly_fit, size_arr, diff_arr, maxfev=1000000)
    popt_poly_corr, _ = curve_fit(crit_poly_fit_corr, size_arr, diff_arr, p0=(1, popt_ising[0], popt_ising[1],
                                                                              popt_ising[2]), maxfev=10000000)
    nu = 1 / popt[0]
    nu_poly_corr = popt_poly_corr[0]
    print("FITTING RESULTS:")
    print("nu = ", nu)
    print("nu_poly_corr = ", popt_poly_corr)

    config = {
        "ylabelsize": 14,
        "legendfontsize": 14
    }

    fig, ax = plt.subplots(1, 1, figsize = figsize)
    # Setze Tickmarken und Labels
    ax.tick_params(direction='in', which='both', length=6, width=2, labelsize=9)
    ax.tick_params(direction='in', which='minor', length=3, width=1, labelsize=9)




    # TODO minor locator muss
    span = np.max(diff_arr) - np.min(diff_arr)

    ax.yaxis.set_major_locator((plt.MultipleLocator(span / 4)))
    ax.yaxis.set_minor_locator((plt.MultipleLocator(span / 4 / 5)))
    # x_span = np.max(xi_size_arr) - np.min(xi_size_arr)
    #ax.xaxis.set_major_locator(ticker.MultipleLocator(base=np.maximum(int(x_span/5), 0.5)))
    # Füge Gitterlinien hinzu
    ax.grid(which='major', linestyle='--', alpha=0.5)
    ax.plot(size_arr, diff_arr, linestyle="", marker="+")
    L_fit = np.linspace(0, np.max(size_arr) + 0.2 * np.max(size_arr), 101)
    ax.plot(L_fit, poly(L_fit, 1 / nu, np.exp(popt[1])), label=rf"$\nu = {nu:.2f}$", color=colors[0])
    # ax.plot(L_fit, crit_poly_fit_corr(L_fit, *popt_poly_corr), label=rf"$\nu = {nu_poly_corr:.2f}$")
    # ax.plot(L_fit, ising_corr_poly_fit(L_fit, *popt_ising), label=rf"$\omega = {popt_ising[2]:.2f}$", linestyle="-.")
    ax.set_xlabel("L")
    ax.set_ylabel(r"$\frac{d U_L}{d \varepsilon}$")
    ax.legend()
    ax.set_title(r"$\frac{d U_L}{d \varepsilon}$ for different System sizes $L$")
    configure_ax(fig, ax, config)
    # save_plot(root, "/critical_exponent.pdf", format="pdf")
    fig.savefig(root + "/critical_exponent.png", format="png", dpi=250, transparent=transparent_plots)
    plt.show()


    # fitting L/xi



    # show nu vs L_max

    L_max_upper = np.minimum(L_max_upper, np.max(size_arr))
    nu_arr = []
    L_max_arr = np.arange(L_max_lower + 2, L_max_upper+1)
    for L_max in L_max_arr:
        diff_fit_arr = diff_arr[(size_arr <= L_max) & (size_arr >= L_max_lower)]
        size_fit_arr = size_arr[(size_arr <= L_max) & (size_arr >= L_max_lower)]
        popt, _ = curve_fit(linear_fit, np.log(size_fit_arr), np.log(diff_fit_arr))
        nu = 1 / popt[0]
        nu_arr.append(nu)


    fig, ax = plt.subplots(1, 1)
    ax.plot(L_max_arr, nu_arr, label=r"$L_{min} = $" + str(L_max_lower), linestyle="", marker="x")
    ax.set_xlabel(r"$L_{max}$")
    ax.set_ylabel(r"$\nu$")
    configure_ax(fig, ax)
    plt.show()



def plt_inter(ax, fig, x, y_dic, res=1000, r=1):
    """
    plots the intersection of datasets in a dictionary
    :param res: number of datapoints to interpolate
    :param x: value which the data depends on (Temperature in our case)
    :param y_dic: dictionary with data for different cases (correlation length on Temperature for different system sizes)
    :return: None
    """
    # interpolate the x-values
    x_inter = np.linspace(x[0], x[-1], res)
    # create interpolations of the data for every size
    y_dic_inter = {}
    for key in y_dic.keys():
        print(x)
        print(y_dic[key])
        y_dic_inter[key] = np.interp(x_inter, x, y_dic[key])
    # now we find the intersection point?
    x_intersec, y_intersec, intersec_ind = det_intersection(x_inter, y_dic_inter)

    # plotting
    shown_inds = [int(i) for i in np.linspace(0, len(y_dic.keys()) - 1, r)]
    for i, key in enumerate(y_dic.keys()):
        if i in shown_inds:
            ax.plot(x, y_dic[key], ls="", marker="x", c=colors[i])         # discrete points
            ax.plot(x_inter, y_dic_inter[key], c=colors[i], label=key, linewidth=1)   # interpolated line
    mark_point(ax, x_intersec, y_intersec)

    return x_intersec, y_intersec

def plot_intersection(T, r, root, xix_dic, xiy_dic):
    figsize_factor = 1
    L_xi_dic = {}
    interpol_L_xi_dic = {}
    figsize = ( 1.2 * 6.4,4.8)
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    config = {}
    #    "ylabelsize": 20,
    #    "xlabelsize": 22,
    #    "titlesize": 25,
    #    "ticklength" : 9,
    #    "tickwidth" : 3,
    #}
    # Setze Tickmarken und Labels
    ax.tick_params(direction='in', which='both', length=9, width=3, labelsize=9)
    ax.tick_params(direction='in', which='minor', length=6, width=2, labelsize=9)
    x_span = np.max(T) - np.min(T)
    ax.xaxis.set_major_locator(ticker.MultipleLocator(base=np.round(x_span / 2, 2)))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(base=np.round(x_span / 2, 2) / 5))
    # TODO minor locator muss
    # Füge Gitterlinien hinzu
    ax.grid(which='major', linestyle='--', alpha=0.5)
    T_inter_arr = np.linspace(np.min(T), np.max(T), 300)
    line_nr = 0
    for i, size in enumerate(xix_dic.keys()):
        # only plot every n-th
        L_xi_dic[size] = int(size) / (0.5 * (xix_dic[size] + xiy_dic[size]))
        interpol_L_xi_dic[size] = np.interp(T_inter_arr, T, L_xi_dic[size])
        print(i)
        if i % r == 0:
            # print(size, cum_dic[size])
            # interpolation
            ax.errorbar(T, L_xi_dic[size], yerr=None, ls="",
                        marker="x", color=colors[2 * line_nr], ecolor="black", capsize=3)
            ax.plot(T_inter_arr, interpol_L_xi_dic[size], color=colors[2 * line_nr],
                    label=rf"L = {size}", linewidth=0.5 * figsize_factor)

            line_nr += 1
    # get intersection
    y_span = 0
    for key in L_xi_dic.keys():
        y_span = np.maximum(np.max(L_xi_dic[key]) - np.min(L_xi_dic[key]), y_span)

    T_xi_inter, xi_inter, xi_inter_ind = det_intersection(T_inter_arr, interpol_L_xi_dic)
    # mark_point(ax, T_xi_inter, xi_inter)

    plt.xticks(list(plt.xticks()[0]) + [T_xi_inter])
    ax.set_xlim(np.min(T) - 0.02 * np.mean(T), np.max(T) + 0.02 * np.mean(T))
    ax.set_xlabel(r"$\varepsilon$")
    ax.set_ylabel(r"$L/\xi$")
    ax.set_title(r" $L/\xi$ on $\varepsilon$")
    ax.set_ylim((0, ax.get_ylim()[1]))
    configure_ax(fig, ax, config)
    # for dummy plot ticks
    # x_span = np.max(T_inter_arr) - np.min(T_inter_arr)
    # tick_length = x_span / 5
    # x_ticks = np.linspace(T_xi_inter - 5 * tick_length, T_xi_inter + 5 * tick_length, 11, endpoint=True)
    # x_minor_ticks = np.linspace(T_xi_inter - 5 * tick_length, T_xi_inter + 5 * tick_length, 5 * (len(x_ticks) - 1) + 1, endpoint=True)
    # x_tick_labels = ["" for tick in x_ticks]
    # x_tick_labels[len(x_ticks) // 2] = "0"
    # ax.set_xticks(x_ticks, x_tick_labels, fontsize=20)
    # ax.set_xticks(x_minor_ticks, minor=True)
    # plt.tick_params(labelleft= False)
    # ax.get_legend().remove()
    # ax.set_xlim((np.min(T_inter_arr) - 0.05 * x_span, np.max(T_inter_arr) + 0.05 * x_span))

    fig.savefig(root + "/L_xi.png", format="png", dpi=300, transparent=True)
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