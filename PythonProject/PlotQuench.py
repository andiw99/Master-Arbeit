from Suite import quench_measurement
import matplotlib.pyplot as plt
from FunctionsAndClasses import *
def main():
    #simpath = "../../Generated content/Silicon/Subsystems/Suite/h/0.4161791450287818/Quench-eta/0.01"
    simpath = "../../Generated content/Silicon/Subsystems/Suite/h/1.7320508075688776/Quench/"
    #simpath = "../../Generated content/Final/Quenches-old/1/Damping/Quench/0.1"
    simpath = "../../Generated content/Final/Quenches-old/1/Damping/Quench/0.05"

    taus = [0.000002, 0.000488]
    taus = [4096, 1024, 32]
    # h = 0.4
    xi_ampl = 4.2
    xi_ampl_perp = 0.2
    Tc = 1.975
    xi_ampl = 1.2                   # h = 1.7
    xi_ampl_perp = 0.2
    Tc = 1.2
    additional_ft_points = 50
    min_tau = 50

    cut_from_equil = 0.2
    cut_zero_impuls = True
    fitfunc = lorentz_offset

    fig, axes = quench_measurement.plot_quench_process(simpath, taus, xi_ampl, Tc, cut_from_equil=cut_from_equil, direction="parallel")
    create_directory_if_not_exists(simpath + f"/plots/")
    plt.savefig(simpath + f"/plots/quench-process-tau-{taus}.png", format="png")
    plt.show()

    fig, axes = quench_measurement.plot_quench_process(simpath, taus, xi_ampl_perp, Tc, cut_from_equil=cut_from_equil,
                                                       direction="perp")
    plt.savefig(simpath + f"/plots/quench-process-perp-tau-{taus}.png", format="png")
    plt.show()

    size_tau_xix_dic, size_tau_xiy_dic = quench_measurement.get_size_quench_results(simpath,
                                                                                    cut_zero_impuls=cut_zero_impuls,
                                                                                    fitfunc=fitfunc,
                                                                                    additional_ft_points=additional_ft_points)
    tau_scaling, xix_scaling, reg_x, max_tau_ind_x, min_tau_ind_x = quench_measurement.fit_kzm(
        size_tau_xix_dic, min_tau=min_tau)

    tau_scaling, xiy_scaling, reg_y, max_tau_ind_y, min_tau_ind_y = quench_measurement.fit_kzm(
        size_tau_xiy_dic, min_tau=min_tau)

    quench_measurement.plot_kzm_scaling(tau_scaling, size_tau_xix_dic, reg_x, max_tau_ind_x, min_tau_ind_x, direction="parallel")
    plt.savefig(simpath + f"/plots/tau-xi-parallel.png", format="png")
    plt.show()

    quench_measurement.plot_kzm_scaling(tau_scaling, size_tau_xiy_dic, reg_y, max_tau_ind_y, min_tau_ind_y, direction="perp",
                                        color=colors[3])
    plt.savefig(simpath + f"/plots/tau-xi-perp.png", format="png")
    plt.show()

    quench_measurement.plot_ratio_after_quench(tau_scaling, xix_scaling, xiy_scaling, min_tau_ind_x)
    plt.savefig(simpath + "/plots/xix_xiy.png", format="png")
    plt.show()

    exit()
    fig, axes = quench_measurement.plot_quench_process(simpath, taus, xi_ampl, Tc)
    new_xlim = (0.3, 0.5)
    axes[1].set_xlim(new_xlim)
    axes[1].set_ylim(5, 20)
    configure_ax(fig, axes[1])
    #ax = zoom_plot(axes[1], new_xlim)
    ax = axes[1]
    ax.set_xlim(new_xlim)
    plt.show()

if __name__ == "__main__":
    main()