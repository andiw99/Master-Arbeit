import matplotlib.pyplot as plt

from Suite import *
import numpy as np

def main():
    # I want to do multiple Quenches of which I do not necissarily know the critical
    # temperature and I want to vary the h
    # If we choose our old values still, the h should go up to 30 which would be
    # the relation of J_parallel and h in the real system
    nr_gpus = 20
    # we somehow need the relevant parameters
    # The model defining parameters are J_perp J_para h eta
    # the simulation defining parameters are dt
    J_para = -10
    J_perp = -0.1
    h = 1
    eta = 1
    p = 2.5
    dt = 0.01

    project = "MinimalCudaProject"
    filepath = f"/home/weitze73/Documents/Master-Arbeit/Code/Master-Arbeit/{project}"
    filepath = f"/home/andi/Studium/Code/Master-Arbeit/{project}"
    simulation_path = "../../Generated content/Paper content/Binder intersection/nu scalar/"

    Tc_exec_file = "AutoCumulantOBC.cu"
    #runfile = "run_cuda_gpu_a100_low.sh"
    runfile = "run_cuda.sh"
    runfile = "run_cuda_gpu_a100_low_minimal.sh"
    # Tc parameters
    nr_Ts = 3
    # We use small equilibration errors since we want to have really accurate
    equil_error = 0.0004
    min_equil_error = 0.0004
    max_rel_intersection_error = 0.005       # is this to small or fine?
    equil_cutoff = 0.1
    # since we need quantitative exact values, we should know T_c beforehand
    #min_T = 0.94
    #max_T = 0.96
    # Values of J_J = 31
    # min_T = 0.98 * (0.853073)
    # max_T = 1.02 * (0.853073)

    min_T = 0.97 * (1.966)
    max_T = 1.03 * (1.966)
    nr_Ts = 5

    # we should start at another parameter file nr because yeah
    para_nr = 230
    file_ending = "mag"
    value_name = "U_L"
    val_write_density = 1 / 100
    min_mag_nr = 750
    process_file_func = process_new_mag_file_to_U_L
    equil_cutoff = 0.2



    Ls = [56, 84, 96, 112, 128]

    Ls = [84, 112]
    Ls = [32, 40, 48, 56, 64, 72, 80, 96, 112, 128]

    Ls = [32, 36, 40, 44, 48, 52, 56, 60, 64, 68, 72, 80]
    Ly_Lx = 1 / 2

    crit_temps = []
    crit_temp_errors = []
    Tc_sim = crit_temp_measurement(J_para, J_perp, h, eta, p, dt, filepath,
                                             simulation_path, Tc_exec_file, runfile=runfile, nr_GPUS=nr_gpus,
                                             sizes=Ls, equil_error=equil_error, Ly_Lx=Ly_Lx,
                                             T_min=min_T, T_max=max_T, nr_Ts=nr_Ts, min_equil_error=min_equil_error,
                                             intersection_error=max_rel_intersection_error, para_nr=para_nr,
                                             min_val_nr=min_mag_nr, file_ending=file_ending, value_name=value_name,
                                             process_file_func=process_file_func, val_write_density=val_write_density,
                                             equil_cutoff=equil_cutoff, project=project)
    T_c, T_c_error = Tc_sim.routine()
    crit_temps.append(T_c)
    crit_temp_errors.append(T_c_error)
    print("\n\n\n")
    print("Tcs", crit_temps)
    print("T_c_error", crit_temp_errors)
    print("\n\n\n")
    # plot the L dependence
    #all_Ls = np.concatenate((small_Ls, small_Ls_ef))
    all_Ls = Ls
    I = 1300        # moment of inertia

    crit_temps = np.array(crit_temps) / np.abs(J_para)
    crit_temp_errors = np.array(crit_temp_errors) / np.abs(J_para)
    fig, ax = plt.subplots(1, 1)

    config = {
        "labelrotation" : 90,
        "increasefontsize" : 0.3
    }

    ax.errorbar(all_Ls, crit_temps, yerr=crit_temp_errors, ecolor="black", elinewidth=1,  capsize=2, linestyle="None", marker="s", markerfacecolor="None", markeredgecolor="C0")
    low_lim = ax.get_xlim()[0]
    high_lim = ax.get_xlim()[1]
    ax.hlines(crit_temps[-1] * 0.9999, low_lim, high_lim, linestyles="dashed", color="C0", label=f"$T_c = {(crit_temps[-1] * 0.999):.4f}$")
    ax.set_xlim(low_lim, high_lim)
    #ax.set_title("Intersection of $L_{min}$ and $L_{max} = 2 L_{min}$ depending on $L_{min}$")
    ax.set_xlabel(r"$L_1$")
    ax.set_ylabel(r"$T_c~/~J_\parallel$")
    configure_ax(fig, ax, config)
    plt.savefig(f"{simulation_path}/plots/Tc-L.png", format="png", dpi=300)
    plt.show()

if __name__ == "__main__":
    main()