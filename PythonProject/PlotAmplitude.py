from Suite import amplitude_measurement
import matplotlib.pyplot as plt
from FunctionsAndClasses import *
def main():

    equil_cutoff = 0.1

    #simpath_xix = "../../Generated content/Silicon/Subsystems/Suite/h/Large Jx/Jx=10-Lx_Ly=32/0.4161791450287818/Amplitude"
    #simpath_xix = "../../Generated content/Silicon/Subsystems/Suite/h/Large Jx/Jx=3-Lx_Ly=32/Amplitude"
    simpath_xix = "../../Generated content/Final/Amplitude/J_J=30/final/Amplitude"
    result_xix = amplitude_measurement.prep_sim_data(equil_cutoff, simpath_xix, "xix")
    #simpath_xiy = "../../Generated content/Silicon/Subsystems/Suite/h/Large Jx/Jx=10-Lx_Ly=1/0.4161791450287818/Amplitude"
    #simpath_xiy = "../../Generated content/Silicon/Subsystems/Suite/h/Large Jx/Jx=3-Lx_Ly=1/Amplitude"
    simpath_xiy = "../../Generated content/Final/Amplitude/J_J=30/final/Amplitude"
    result_xiy = amplitude_measurement.prep_sim_data(equil_cutoff, simpath_xiy, "xiy")
    print(result_xix)
    print(result_xiy)
    fig, axx = plt.subplots(1, 1, figsize=(6.4 * 1.5, 4.8 * 1.5))
    axy = axx.twinx()
    axes = [axx, axy]

    results = [result_xix, result_xiy]


    for j, result in enumerate(results):
        for i, size in enumerate(result):
            # We definitely plot the points
            T_xix, xix_arr = result[size]

            amplitude_measurement.plot_xi_points(axes[j], T_xix,
                                                 xix_arr, marker=markers[i], size=size, color=5 * j)



    fit_sizes_x = list(result_xix.keys())
    T_xix, reg_x = amplitude_measurement.perform_fit_on_sizes(result_xix, fit_sizes_x)
    amplitude_measurement.plot_xi_fit(axes[0], reg_x, T_xix)

    fit_sizes_y = list(result_xiy.keys())
    T_xiy, reg_y = amplitude_measurement.perform_fit_on_sizes(result_xiy, fit_sizes_y)
    amplitude_measurement.plot_xi_fit(axes[1], reg_y, T_xiy, color=5)

    amplitude_measurement.configure_xi_plot(axes, fig)
    plt.show()

    def construct_ratios(results_x, results_y):
        T_x = np.array([])
        xix = np.array([])
        T_y = np.array([])
        xiy = np.array([])
        for size in results_x:
            cur_T, cur_xi = results_x[size]

            T_x = np.concatenate((T_x, cur_T))
            xix = np.concatenate((xix, cur_xi))
        for size in results_y:
            cur_T, cur_xi = results_y[size]

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

    exit(0)
    simpath = "../../Generated content/Silicon/Subsystems/Suite/h/Large Jx/Jx=10-Lx_Ly=1/0.4161791450287818/Amplitude"

    amplitude_measurement.plot_divergence(simpath, Tc=1.7477, direction="xiy")


if __name__ == "__main__":
    main()