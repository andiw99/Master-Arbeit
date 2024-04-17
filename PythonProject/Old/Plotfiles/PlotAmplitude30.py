from Suite import amplitude_measurement
import matplotlib.pyplot as plt
from FunctionsAndClasses import *
def main():
    small_h()

    equil_cutoff = 0.1
    J_para = 3

    T_min = 0.298
    T_min_perp = 0.298

    #simpath_xix = "../../Generated content/Silicon/Subsystems/Suite/h/Large Jx/Jx=10-Lx_Ly=32/0.4161791450287818/Amplitude"
    #simpath_xix = "../../Generated content/Silicon/Subsystems/Suite/h/Large Jx/Jx=3-Lx_Ly=32/Amplitude"
    simpath_xix = "../../Generated content/Final/Amplitude/J_J=30/final/Amplitude"
    result_xix = amplitude_measurement.prep_sim_data(equil_cutoff, simpath_xix, "xix", T_min=T_min)
    #simpath_xiy = "../../Generated content/Silicon/Subsystems/Suite/h/Large Jx/Jx=10-Lx_Ly=1/0.4161791450287818/Amplitude"
    #simpath_xiy = "../../Generated content/Silicon/Subsystems/Suite/h/Large Jx/Jx=3-Lx_Ly=1/Amplitude"
    simpath_xiy = "../../Generated content/Final/Amplitude/J_J=30/final/Amplitude"
    result_xiy = amplitude_measurement.prep_sim_data(equil_cutoff, simpath_xiy, "xiy", T_min=T_min_perp)
    print(result_xix)
    print(result_xiy)
    fig, axx = plt.subplots(1, 1, figsize=(6.4 * 1.5, 4.8 * 1.5))
    axy = axx.twinx()
    axes = [axx, axy]

    results = [result_xix, result_xiy]
    directions = ["parallel", "perp"]
    size_factors = [1, 1 / 8]


    for j, result in enumerate(results):
        for i, size in enumerate(result):
            # We definitely plot the points
            T_xix, xix_arr, xix_err = result[size]
            T_xix /= J_para
            # I guess we also normalize the results?
            #result[size] = (result[size][0] / J_para, result[size][1])
            amplitude_measurement.plot_xi_points(axes[j], T_xix,
                                                 xix_arr, marker=markers[i], size=int(size * size_factors[j]),
                                                 color=colors[5 * j], direction=directions[j], yerr = None)


    print(result_xix)
    print(result_xiy)
    fit_sizes_x = list(result_xix.keys())
    T_xix, reg_x = amplitude_measurement.perform_fit_on_sizes(result_xix, fit_sizes_x, min_points=0, T_min=T_min)
    print("\nT_xix:")
    print(T_xix)
    amplitude_measurement.plot_xi_fit(axes[0], reg_x, T_xix, direction="parallel", Tc_label=False)

    fit_sizes_y = list(result_xiy.keys())
    T_xiy, reg_y = amplitude_measurement.perform_fit_on_sizes(result_xiy, fit_sizes_y, min_points=0, T_min=T_min_perp)
    amplitude_measurement.plot_xi_fit(axes[1], reg_y, T_xiy, color=colors[5], direction="perp", Tc_label=False)

    amplitude_measurement.configure_xi_plot(axes, fig)
    #handles, labels = plt.gca().get_legend_handles_labels()
    ## Rearrange handles and labels (e.g., swap them)
    #handles = [handles[1], handles[0], handles[2], handles[3], handles[4]]  # Swap the order of handles
    #labels = [labels[1], labels[0], labels[2], labels[2], labels[2]]  # Swap the order of labels
    #plt.legend(handles, labels)
    create_directory_if_not_exists(simpath_xix + "/plots")
    fig.savefig(simpath_xix + "/plots/xi-divergence", format="svg")
    fig.savefig(simpath_xix + "/plots/xi-divergence.png", format="png", dpi=300)
    plt.show()

    # Inverse plot
    fig, ax = plt.subplots(1, 1, figsize=(6.4 * 1.5, 4.8 * 1.5))

    for j, result in enumerate(results):
        for i, size in enumerate(result):
            # We definitely plot the points
            T_xi, xi_arr, xi_err = result[size]
            amplitude_measurement.plot_inverse_points(fig, ax, T_xi, xi_arr, marker=markers[i],
                                                      color=colors[5 * j], size=int(size * size_factors[j]), direction=directions[j])
            # Find out which size has the largest number of points and for that we do the fit


        # The thing is one could actually think about performing the fit on
        # combined measurement points, but that is for another time
        # one would have to replace smart the small size values for the large size

    Tc_x = amplitude_measurement.plot_xi_inv_fit(ax, T_xix, reg_x, color=colors[0], plot_until_pt=True, Tc_label=False, direction="parallel")
    Tc_y = amplitude_measurement.plot_xi_inv_fit(ax, T_xiy, reg_y, color=colors[5], plot_until_pt=True, Tc_label=False, direction="perp")
    xlims = ax.get_xlim()
    ax.hlines(0, *xlims, linestyles="--")
    ax.set_xlim(*xlims)

        # Now we still have to configure the inverse plot
    config = amplitude_measurement.configure_xi_inv_plot([ax], fig)
    ax.set_ylim(-0.01, ax.get_ylim()[1])
    configure_ax(fig, ax, config)
    plt.tight_layout()
    fig.savefig(simpath_xix + "/plots/xi-inv-divergence", format="svg")
    fig.savefig(simpath_xix + "/plots/xi-inv-divergence.png", format="png", dpi=300)
    plt.show()
    def construct_ratios(results_x, results_y):
        T_x = np.array([])
        xix = np.array([])
        T_y = np.array([])
        xiy = np.array([])
        for size in results_x:
            cur_T, cur_xi, cur_err = results_x[size]

            T_x = np.concatenate((T_x, cur_T))
            xix = np.concatenate((xix, cur_xi))
        for size in results_y:
            cur_T, cur_xi, cur_err = results_y[size]

            T_y = np.concatenate((T_y, cur_T))
            xiy = np.concatenate((xiy, cur_xi))

        xi_ratios = []
        Ts = []
        for ix in range(len(T_x)):
            T = T_x[ix]
            for iy in range(len(T_y)):
                Ty = T_y[iy]
                if T == Ty:
                    xi_ratios.append(xix[ix] / xiy[iy])
                    Ts.append(T)
        return Ts, xi_ratios

    same_Ts, xi_ratios = construct_ratios(result_xix, result_xiy)

    fig, ax = plt.subplots(1, 1)

    ax.plot(same_Ts, xi_ratios, **blue_square_kwargs, label=r"$\xi_\parallel / \xi_\perp$")

    config = {
        "labelrotation": 90,
    }

    configure_ax(fig, ax, config)
    ax.set_xlabel(r"T")
    ax.set_ylabel(r"$\xi_\parallel / \xi_\perp$")
    plt.show()

    # exit()
    # simpath = "../../Generated content/Final/Amplitude/J_J=30/final-small-h/Amplitude"
    # amplitude_measurement.plot_divergence(simpath, Tc=None, direction="both")

def small_h():
    equil_cutoff = 0.1
    J_para = 3
    T_min = 0.251
    T_min_perp = 0.251

    #simpath_xix = "../../Generated content/Silicon/Subsystems/Suite/h/Large Jx/Jx=10-Lx_Ly=32/0.4161791450287818/Amplitude"
    #simpath_xix = "../../Generated content/Silicon/Subsystems/Suite/h/Large Jx/Jx=3-Lx_Ly=32/Amplitude"
    simpath_xix = "../../Generated content/Final/Amplitude/J_J=30/final-small-h/Amplitude"
    result_xix = amplitude_measurement.prep_sim_data(equil_cutoff, simpath_xix, "xix", T_min=T_min)
    #simpath_xiy = "../../Generated content/Silicon/Subsystems/Suite/h/Large Jx/Jx=10-Lx_Ly=1/0.4161791450287818/Amplitude"
    #simpath_xiy = "../../Generated content/Silicon/Subsystems/Suite/h/Large Jx/Jx=3-Lx_Ly=1/Amplitude"
    simpath_xiy = "../../Generated content/Final/Amplitude/J_J=30/final-small-h/Amplitude"
    result_xiy = amplitude_measurement.prep_sim_data(equil_cutoff, simpath_xiy, "xiy", T_min=T_min_perp)
    print(result_xix)
    print(result_xiy)
    fig, axx = plt.subplots(1, 1, figsize=(6.4 * 1.5, 4.8 * 1.5))
    axy = axx.twinx()
    axes = [axx, axy]

    results = [result_xix, result_xiy]
    directions = ["parallel", "perp"]
    size_factors = [1, 1 / 8]


    for j, result in enumerate(results):
        for i, size in enumerate(result):
            # We definitely plot the points
            T_xix, xix_arr, xix_err = result[size]
            T_xix /= J_para
            # I guess we also normalize the results?
            #result[size] = (result[size][0] / J_para, result[size][1])
            amplitude_measurement.plot_xi_points(axes[j], T_xix,
                                                 xix_arr, marker=markers[i], size=int(size * size_factors[j]),
                                                 color=colors[5 * j], direction=directions[j], yerr = None)


    print(result_xix)
    print(result_xiy)
    fit_sizes_x = list(result_xix.keys())
    T_xix, reg_x = amplitude_measurement.perform_fit_on_sizes(result_xix, fit_sizes_x, min_points=0, T_min=T_min)
    print("\nT_xix:")
    print(T_xix)
    amplitude_measurement.plot_xi_fit(axes[0], reg_x, T_xix, direction="parallel", Tc_label=False)

    fit_sizes_y = list(result_xiy.keys())
    T_xiy, reg_y = amplitude_measurement.perform_fit_on_sizes(result_xiy, fit_sizes_y, min_points=0, T_min=T_min_perp)
    amplitude_measurement.plot_xi_fit(axes[1], reg_y, T_xiy, color=colors[5], direction="perp", Tc_label=False)

    amplitude_measurement.configure_xi_plot(axes, fig)
    #handles, labels = plt.gca().get_legend_handles_labels()
    ## Rearrange handles and labels (e.g., swap them)
    #handles = [handles[1], handles[0], handles[2], handles[3], handles[4]]  # Swap the order of handles
    #labels = [labels[1], labels[0], labels[2], labels[2], labels[2]]  # Swap the order of labels
    #plt.legend(handles, labels)
    create_directory_if_not_exists(simpath_xix + "/plots")
    fig.savefig(simpath_xix + "/plots/xi-divergence", format="svg")
    fig.savefig(simpath_xix + "/plots/xi-divergence.png", format="png", dpi=300)
    plt.show()

    # Inverse plot
    fig, ax = plt.subplots(1, 1, figsize=(6.4 * 1.5, 4.8 * 1.5))

    for j, result in enumerate(results):
        for i, size in enumerate(result):
            # We definitely plot the points
            T_xi, xi_arr, xi_err = result[size]
            amplitude_measurement.plot_inverse_points(fig, ax, T_xi, xi_arr, marker=markers[i],
                                                      color=colors[5 * j], size=int(size * size_factors[j]), direction=directions[j])
            # Find out which size has the largest number of points and for that we do the fit


        # The thing is one could actually think about performing the fit on
        # combined measurement points, but that is for another time
        # one would have to replace smart the small size values for the large size

    Tc_x = amplitude_measurement.plot_xi_inv_fit(ax, T_xix, reg_x, color=colors[0], plot_until_pt=True, Tc_label=False, direction="parallel")
    Tc_y = amplitude_measurement.plot_xi_inv_fit(ax, T_xiy, reg_y, color=colors[5], plot_until_pt=True, Tc_label=False, direction="perp")
    xlims = ax.get_xlim()
    ax.hlines(0, *xlims, linestyles="--")
    ax.set_xlim(*xlims)

        # Now we still have to configure the inverse plot
    config = amplitude_measurement.configure_xi_inv_plot([ax], fig)
    ax.set_ylim(-0.01, ax.get_ylim()[1])
    configure_ax(fig, ax, config)
    plt.tight_layout()
    fig.savefig(simpath_xix + "/plots/xi-inv-divergence", format="svg")
    fig.savefig(simpath_xix + "/plots/xi-inv-divergence.png", format="png", dpi=300)
    plt.show()
    def construct_ratios(results_x, results_y):
        T_x = np.array([])
        xix = np.array([])
        T_y = np.array([])
        xiy = np.array([])
        for size in results_x:
            cur_T, cur_xi, cur_err = results_x[size]

            T_x = np.concatenate((T_x, cur_T))
            xix = np.concatenate((xix, cur_xi))
        for size in results_y:
            cur_T, cur_xi, cur_err = results_y[size]

            T_y = np.concatenate((T_y, cur_T))
            xiy = np.concatenate((xiy, cur_xi))

        xi_ratios = []
        Ts = []
        for ix in range(len(T_x)):
            T = T_x[ix]
            for iy in range(len(T_y)):
                Ty = T_y[iy]
                if T == Ty:
                    xi_ratios.append(xix[ix] / xiy[iy])
                    Ts.append(T)
        return Ts, xi_ratios

    same_Ts, xi_ratios = construct_ratios(result_xix, result_xiy)

    fig, ax = plt.subplots(1, 1)

    ax.plot(same_Ts, xi_ratios, **blue_square_kwargs, label=r"$\xi_\parallel / \xi_\perp$")

    config = {
        "labelrotation": 90,
    }

    configure_ax(fig, ax, config)
    ax.set_xlabel(r"T")
    ax.set_ylabel(r"$\xi_\parallel / \xi_\perp$")
    plt.show()


if __name__ == "__main__":
    main()