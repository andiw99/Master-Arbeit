from Suite import quench_measurement
import matplotlib.pyplot as plt
from FunctionsAndClasses import *
def main():
    plot_multiple_eta()
    plot_multiple_h()


def plot_multiple_h():
    simpaths = [
                "../../Generated content/Final/Quenches-old/0.5/Damping/Quench/1",
                "../../Generated content/Final/Quenches-old/1/Damping/Quench/1",
                "../../Generated content/Final/Quenches-old/5/Damping/Quench/1",
                "../../Generated content/Final/Quenches-old/large-h/10/Damping/Quench/1",
                ]
    names = [
        r"$, h = 0.5$",
        r"$, h = 1$",
        r"$, h = 5$",
        r"$, h = 10$",
    ]
    additional_ft_points = 5
    min_tau = 400
    cut_zero_impuls = True
    min_points = 4
    fitfunc = lorentz_offset
    plot_multiple_quench_scalings(additional_ft_points, cut_zero_impuls, fitfunc, min_tau, simpaths,
                                  min_points=min_points, names=names)

simpaths = [
    "../../Generated content/Paper content/Quenches/h=/1/Damping/Quench/1",
    "../../Generated content/Paper content/Quenches/h=/1/Damping/Quench/0.01",
    "../../Generated content/Paper content/Quenches/h=/1/Damping/Quench/10",
    #"../../Generated content/Final/Quenches-old/1/Damping/Quench/10",
    #"../../Generated content/Final/Quenches-old/1/Damping/Quench/50-short-equil",
    #"../../Generated content/Final/Quenches-old/1/Damping/Quench/100",
    #"../../Generated content/Final/Quenches-old/1/Damping/Quench/0.1",
    #"../../Generated content/Final/Quenches-old/end-equilibration-correct-range/1/Damping/Quench/0.01",
]
def plot_multiple_eta():

    names = [
        #" equil",
        #" no equil",
        r"$, \eta = 1$",
        r"$, \eta = 0.01$",
        r"$, \eta = 10$",
        r"$, \eta = 0.1$",
    ]
    additional_ft_points = 20
    min_tau = 600
    cut_zero_impuls = True
    min_points = 4
    fitfunc = lorentz_offset
    plot_multiple_quench_scalings_combined(additional_ft_points, cut_zero_impuls,
                                  fitfunc, min_tau, simpaths,
                                  min_points=min_points, names=names)

def plot_multiple_quench_scalings_combined(additional_ft_points, cut_zero_impuls, fitfunc, min_tau, simpaths, min_points=4, names=None):
    fig, ax = None, None
    config = {
        "increasefontsize" : 0.6,
        "labelverticalalignment": "bottom",
        "labelrotation": 90
    }
    for i, simpath in enumerate(simpaths):
        size_tau_xix_dic, size_tau_xiy_dic = (
            quench_measurement.get_size_quench_results(simpath,
                                                    cut_zero_impuls=cut_zero_impuls,
                                                    fitfunc=fitfunc,
                                                    additional_ft_points=additional_ft_points))
        tau_scaling, xix_scaling, reg_x, max_tau_ind_x, min_tau_ind_x = quench_measurement.fit_kzm(
            size_tau_xix_dic, min_tau=min_tau, min_points=min_points)
        tau_scaling, xiy_scaling, reg_y, max_tau_ind_y, min_tau_ind_y = quench_measurement.fit_kzm(
            size_tau_xiy_dic, min_tau=min_tau, min_points=min_points)
        if names:
            name = names[i]
        else:
            name = ""
        fig, ax = quench_measurement.plot_kzm_scaling(tau_scaling, size_tau_xix_dic, reg_x, max_tau_ind_x, min_tau_ind_x,
                                            direction="parallel", fig=fig, ax=ax, color=colors[2 * i], label=False,
                                                        marker=markers[0], last_tau=False, name=name, config=config)
        quench_measurement.plot_kzm_scaling(tau_scaling, size_tau_xiy_dic,
                                            reg_y, max_tau_ind_y, min_tau_ind_y,
                                            direction="perp", color=colors[2 * i], fig=fig, ax=ax, label=False, marker=markers[1],
                                            last_tau=False, name=name, config=config)
    ax.plot([], [], color=colors[0], marker="s", linestyle="", label=r"$\eta = 1$")
    ax.plot([], [], color=colors[2], marker="s", linestyle="", label=r"$\eta = 0.01$")
    ax.plot([], [], color=colors[4], marker="s", linestyle="", label=r"$\eta = 10$")
    ax.plot([], [], **get_point_kwargs_color("black"), marker=markers[0], label=r"$\delta = \parallel$")
    ax.plot([], [], **get_point_kwargs_color("black"), marker=markers[1], label=r"$\delta = \perp$")
    ax.legend()
    ax.set_xlim((ax.get_xlim()[0], 10 ** 4 * 4))
    leg_x = ax.get_legend()
    for line in leg_x.get_lines():
        line.set_linewidth(3.0)


    parent_path = os.path.dirname(simpath)
    create_directory_if_not_exists(parent_path + "/plots/")
    fig.savefig(parent_path + f"/plots/tau-xi-parallel.png", format="png")
    fig.savefig(parent_path + f"/plots/tau-xi-parallel.svg", format="svg")

    plt.show()


def plot_multiple_quench_scalings(additional_ft_points, cut_zero_impuls, fitfunc, min_tau, simpaths, min_points=4, names=None):
    figx, axx = None, None
    figy, axy = None, None
    config = {
        "increasefontsize" : 0.6,
        "labelverticalalignment": "bottom",
        "labelrotation": 90
    }
    for i, simpath in enumerate(simpaths):
        size_tau_xix_dic, size_tau_xiy_dic = (
            quench_measurement.get_size_quench_results(simpath,
                                                    cut_zero_impuls=cut_zero_impuls,
                                                    fitfunc=fitfunc,
                                                    additional_ft_points=additional_ft_points))
        tau_scaling, xix_scaling, reg_x, max_tau_ind_x, min_tau_ind_x = quench_measurement.fit_kzm(
            size_tau_xix_dic, min_tau=min_tau, min_points=min_points)
        tau_scaling, xiy_scaling, reg_y, max_tau_ind_y, min_tau_ind_y = quench_measurement.fit_kzm(
            size_tau_xiy_dic, min_tau=min_tau, min_points=min_points)
        if names:
            name = names[i]
        else:
            name = ""
        figx, axx = quench_measurement.plot_kzm_scaling(tau_scaling, size_tau_xix_dic, reg_x, max_tau_ind_x, min_tau_ind_x,
                                            direction="parallel", fig=figx, ax=axx, color=colors[2 * i], label=False,
                                                        marker=markers[i], last_tau=False, name=name, config=config)
        figy, axy = quench_measurement.plot_kzm_scaling(tau_scaling, size_tau_xiy_dic,
                                            reg_y, max_tau_ind_y, min_tau_ind_y,
                                            direction="perp",
                                            color=colors[3 + i], fig=figy, ax=axy, label=False, marker=markers[i],
                                                        last_tau=False, name=name, config=config)
    leg_x = axx.get_legend()
    leg_y = axy.get_legend()
    for line in leg_x.get_lines():
        line.set_linewidth(3.0)
    for line in leg_y.get_lines():
        line.set_linewidth(3.0)

    parent_path = os.path.dirname(simpath)
    create_directory_if_not_exists(parent_path + "/plots/")
    figx.savefig(parent_path + f"/plots/tau-xi-parallel.png", format="png")
    figx.savefig(parent_path + f"/plots/tau-xi-parallel.svg", format="svg")

    plt.show()

    figy.savefig(parent_path + f"/plots/tau-xi-perp.png", format="png")
    figy.savefig(parent_path + f"/plots/tau-xi-perp.svg", format="svg")
    plt.show()


if __name__ == "__main__":
    main()