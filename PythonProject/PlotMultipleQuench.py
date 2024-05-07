from Suite import quench_measurement
import matplotlib.pyplot as plt
from FunctionsAndClasses import *
def main():
    simpaths = [
                "../../Generated content/Final/Quenches-old/0.5/Damping/Quench/1",
                "../../Generated content/Final/Quenches-old/1/Damping/Quench/1",
                "../../Generated content/Final/Quenches-old/5/Damping/Quench/1",
                "../../Generated content/Final/Quenches-old/large-h/10/Damping/Quench/1",
                ]
    additional_ft_points = 5
    min_tau = 400
    cut_zero_impuls = True
    min_points = 4
    fitfunc = lorentz_offset
    plot_multiple_quench_scalings(additional_ft_points, cut_zero_impuls, fitfunc, min_tau, simpaths, min_points=min_points)


def plot_multiple_eta():
    simpaths = [
        "../../Generated content/Final/Quenches-old/1/Damping/Quench/1",
        "../../Generated content/Final/Quenches-old/1/Damping/Quench/10",
        "../../Generated content/Final/Quenches-old/1/Damping/Quench/0.1",
        "../../Generated content/Final/Quenches-old/1/Damping/Quench/0.01",
    ]
    additional_ft_points = 5
    min_tau = 400
    cut_zero_impuls = True
    min_points = 4
    fitfunc = lorentz_offset
    plot_multiple_quench_scalings(additional_ft_points, cut_zero_impuls,
                                  fitfunc, min_tau, simpaths,
                                  min_points=min_points)

def plot_multiple_quench_scalings(additional_ft_points, cut_zero_impuls, fitfunc, min_tau, simpaths, min_points=4):
    figx, axx = None, None
    figy, axy = None, None
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
        figx, axx = quench_measurement.plot_kzm_scaling(tau_scaling, size_tau_xix_dic, reg_x, max_tau_ind_x, min_tau_ind_x,
                                            direction="parallel", fig=figx, ax=axx, color=colors[i], label=False)
        figy, axy = quench_measurement.plot_kzm_scaling(tau_scaling, size_tau_xiy_dic,
                                            reg_y, max_tau_ind_y, min_tau_ind_y,
                                            direction="perp",
                                            color=colors[3 + i], fig=figy, ax=axy, label=False)
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