from Suite import *
import numpy as np

def main():
    h_arr = [10]
    nr_gpus = 20
    J_para = -10
    J_perp = -0.1

    p = 2.5
    eta_arr = [0.01]
    eta_arr = [50, 100]
    eta_arr = [1]
    dt = 0.01

    project = "MinimalCudaProject"
    filepath = f"/home/weitze73/Documents/Master-Arbeit/Code/Master-Arbeit/{project}"
    filepath = f"/home/andi/Studium/Code/Master-Arbeit/{project}"
    # simulation_path = "../../Generated content/h Quench/"
    simulation_path = "../../Generated content/Paper content/Quenches/h=/"
    equil_time_end = 0
    equil_time_start = 500
    gamma = 1
    base = np.sqrt(2)

    quench_exec_file = "AutoQuench.cu"
    runfile_quench = "run_cuda_gpu_a100_low_minimal.sh"
    runfile_quench = "run_cuda.sh"


    # Quench parameters
    max_size = 10000
    min_nr_sites = 1e6
    max_nr_quench_steps = 1e7
    nr_sites = 2e6
    max_tau = 10000
    min_nr_systems = 30
    Ly_Lx = 1 / 8

    #T_c = 21351     # maybe this T_c is to low?
    T_c = 0.296 * 3.11
    T_c = 0.3323 * 10
    #T_c = 1.15
    min_nr_corr_values = 100

    for h in h_arr:
        curr_sim_path = simulation_path + f"{h}/"
        for eta in eta_arr:
            input(f"h = {h}, eta = {eta}, Tc = {T_c} okay?")
            para_nr_quench = int(input("please just change the parameter nubmer :("))
            quench = quench_measurement(J_para, J_perp, h, eta, p, dt, filepath, curr_sim_path + f"Damping/Quench/{eta}",
                                        quench_exec_file, runfile_quench, T_c, nr_GPUS=nr_gpus, size_max=max_size,
                                        min_nr_sites=min_nr_sites, max_nr_steps=max_nr_quench_steps,
                                        para_nr=para_nr_quench, tau_max=max_tau, nr_sites=nr_sites,
                                        min_nr_systems=min_nr_systems, min_nr_corr_values=min_nr_corr_values,
                                        Ly_Lx=Ly_Lx, equil_time=equil_time_start, equil_time_end=equil_time_end,
                                        project=project, gamma=gamma, base=base)
            quench.run()

        # the good thing is, both of the simulation implement pickup capabilities so
        # I dont need to worry to much that my computer loses connection and stuff (which it will)


if __name__ == "__main__":
    main()