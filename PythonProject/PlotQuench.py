from Suite import quench_measurement
import matplotlib.pyplot as plt
from FunctionsAndClasses import *
def main():
    #plot_h_5_eta_1()
    ##exit()
    #plot_h_05_eta_1()
    #plot_h_1_eta_10()
    #plot_h_1_eta_001()
    #plot_h_1_eta_1()
    simpath = "../../Generated content/Final/Quenches-old/process-final1/Damping/Quench/1"
    simpath = "../../Generated content/Final/Quenches-old/large-h/10/Damping/Quench/1"
    taus = [2048, 256, 32]
    xi_ampl = 1.987                 # h = 1.7
    xi_ampl_perp = 0.194
    Tc = 1.972

    additional_ft_points = 50
    min_tau = 400
    cut_from_equil = 0
    cut_zero_impuls = True
    min_points = 4
    fitfunc = lorentz_offset
    yscale= "log"
    plot_quench(Tc, additional_ft_points, cut_from_equil, cut_zero_impuls, fitfunc, min_tau, simpath, taus, xi_ampl,
                xi_ampl_perp, min_points=min_points, yscale=yscale)

    #plot_h_05_eta_1()
    #plot_h_5_eta_1()

def plot_h_1_eta_10():
    simpath = "../../Generated content/Final/Quenches-old/1/Damping/Quench/10"

    taus = [2048, 256, 32]
    xi_ampl = 1.987                 # h = 1.7
    xi_ampl_perp = 0.194
    Tc = 1.972

    additional_ft_points = 50
    min_tau = 100
    cut_from_equil = 0.85
    cut_zero_impuls = True
    min_points = 4
    fitfunc = lorentz_offset
    yscale= "log"
    plot_quench(Tc, additional_ft_points, cut_from_equil, cut_zero_impuls, fitfunc, min_tau, simpath, taus, xi_ampl,
                xi_ampl_perp, min_points=min_points, yscale=yscale)

def plot_h_1_eta_001():
    simpath = "../../Generated content/Final/Quenches-old/1/Damping/Quench/0.01"

    taus = [2048, 256, 32]
    xi_ampl = 1.987                 # h = 1.7
    xi_ampl_perp = 0.194
    Tc = 1.972

    additional_ft_points = 50
    min_tau = 0
    cut_from_equil = 0.85
    cut_zero_impuls = True
    min_points = 4
    fitfunc = lorentz_offset
    yscale= "log"
    plot_quench(Tc, additional_ft_points, cut_from_equil, cut_zero_impuls, fitfunc, min_tau, simpath, taus, xi_ampl,
                xi_ampl_perp, min_points=min_points, yscale=yscale)

def plot_h_05_eta_1():
    simpath = "../../Generated content/Final/Quenches-old/0.5/Damping/Quench/1"

    taus = [2048, 1024, 32]
    xi_ampl = 1.2                   # h = 1.7
    xi_ampl_perp = 0.2
    Tc = 1.2

    additional_ft_points = 50
    min_tau = 20
    cut_from_equil = 0.2
    cut_zero_impuls = True
    min_points = 4
    fitfunc = lorentz_offset
    plot_quench(Tc, additional_ft_points, cut_from_equil, cut_zero_impuls, fitfunc, min_tau, simpath, taus, xi_ampl,
                xi_ampl_perp, min_points=min_points)
def plot_h_1_eta_1():
    simpath = "../../Generated content/Final/Quenches-old/1/Damping/Quench/1"

    taus = [2048, 256, 32]
    xi_ampl = 1.987                 # h = 1.7
    xi_ampl_perp = 0.194
    Tc = 1.972

    additional_ft_points = 50
    min_tau = 0
    cut_from_equil = 0.85
    cut_zero_impuls = True
    min_points = 4
    fitfunc = lorentz_offset
    yscale= "log"
    plot_quench(Tc, additional_ft_points, cut_from_equil, cut_zero_impuls, fitfunc, min_tau, simpath, taus, xi_ampl,
                xi_ampl_perp, min_points=min_points, yscale=yscale)
def plot_h_5_eta_1():
    taus = [ 1024, 32]
    xi_ampl = 1.2                   # h = 1.7
    xi_ampl_perp = 0.2
    Tc = 1.2
    additional_ft_points = 50
    min_tau = 50

    cut_from_equil = 0.2
    cut_zero_impuls = True
    fitfunc = lorentz_offset
    simpath = "../../Generated content/Final/Quenches-old/5/Damping/Quench/1"

    plot_quench(Tc, additional_ft_points, cut_from_equil, cut_zero_impuls, fitfunc, min_tau, simpath, taus, xi_ampl,
                xi_ampl_perp)

def plot_quench(Tc, additional_ft_points, cut_from_equil, cut_zero_impuls, fitfunc, min_tau, simpath, taus, xi_ampl,
                xi_ampl_perp, min_points=4, yscale="log"):
    fig, axes = quench_measurement.plot_quench_process(simpath, taus, xi_ampl, Tc, cut_from_equil=cut_from_equil,
                                                       direction="parallel", yscale=yscale)
    create_directory_if_not_exists(simpath + f"/plots/")
    plt.savefig(simpath + f"/plots/quench-process-tau-{taus}.png", format="png")
    plt.savefig(simpath + f"/plots/quench-process-tau-{taus}.svg", format="svg")
    plt.show()
    fig, axes = quench_measurement.plot_quench_process(simpath, taus, xi_ampl_perp, Tc, cut_from_equil=cut_from_equil,
                                                       direction="perp", yscale=yscale)
    plt.savefig(simpath + f"/plots/quench-process-perp-tau-{taus}.png", format="png")
    plt.savefig(simpath + f"/plots/quench-process-perp-tau-{taus}.svg", format="svg")
    plt.show()
    size_tau_xix_dic, size_tau_xiy_dic = quench_measurement.get_size_quench_results(simpath,
                                                                                    cut_zero_impuls=cut_zero_impuls,
                                                                                    fitfunc=fitfunc,
                                                                                    additional_ft_points=additional_ft_points)
    tau_scaling, xix_scaling, reg_x, max_tau_ind_x, min_tau_ind_x = quench_measurement.fit_kzm(
        size_tau_xix_dic, min_tau=min_tau, min_points=min_points)
    tau_scaling, xiy_scaling, reg_y, max_tau_ind_y, min_tau_ind_y = quench_measurement.fit_kzm(
        size_tau_xiy_dic, min_tau=min_tau, min_points=min_points)
    quench_measurement.plot_kzm_scaling(tau_scaling, size_tau_xix_dic, reg_x, max_tau_ind_x, min_tau_ind_x,
                                        direction="parallel")
    plt.savefig(simpath + f"/plots/tau-xi-parallel.png", format="png")
    plt.savefig(simpath + f"/plots/tau-xi-parallel.svg", format="svg")

    plt.show()
    quench_measurement.plot_kzm_scaling(tau_scaling, size_tau_xiy_dic, reg_y, max_tau_ind_y, min_tau_ind_y,
                                        direction="perp",
                                        color=colors[3])
    plt.savefig(simpath + f"/plots/tau-xi-perp.png", format="png")
    plt.savefig(simpath + f"/plots/tau-xi-perp.svg", format="svg")
    plt.show()
    quench_measurement.plot_ratio_after_quench(tau_scaling, xix_scaling, xiy_scaling, min_tau_ind_x)
    plt.savefig(simpath + "/plots/xix_xiy.png", format="png")
    plt.show()
    # exit()
    # fig, axes = quench_measurement.plot_quench_process(simpath, taus, xi_ampl, Tc)
    # new_xlim = (0.3, 0.5)
    # axes[1].set_xlim(new_xlim)
    # axes[1].set_ylim(5, 20)
    # configure_ax(fig, axes[1])
    # # ax = zoom_plot(axes[1], new_xlim)
    # ax = axes[1]
    # ax.set_xlim(new_xlim)
    # plt.show()


if __name__ == "__main__":
    main()