from Suite import *
import numpy as np

def main():
    h_arr = [1]
    nr_gpus = 20
    J_para = -10
    J_perp = -0.1

    p = 2.5
    eta_arr = [50, 100]
    eta_arr = [0.01]
    dt = 0.01

    filepath = "/home/weitze73/Documents/Master-Arbeit/Code/Master-Arbeit/CudaProject"
    filepath = "/home/andi/Studium/Code/Master-Arbeit/CudaProject"
    # simulation_path = "../../Generated content/h Quench/"
    simulation_path = "../../Generated content/Final/Quenches-old/end-equilibration-longer/"
    equil_time_end = 500
    equil_time_start = 50
    gamma = 2

    quench_exec_file = "AutoQuench.cu"
    runfile_quench = "run_cuda.sh"
    runfile_quench = "run_cuda_gpu_a100_low_minimal.sh"
    project = "MinimalCudaProject"


    # Quench parameters
    max_size = 10000
    min_nr_sites = 1e6
    max_nr_quench_steps = 1e7
    nr_sites = 2e6
    max_tau = 100000
    min_nr_systems = 10
    Ly_Lx = 1 / 8
    para_nr_quench = int(input("please just change the parameter nubmer :("))

    #T_c = 21351     # maybe this T_c is to low?
    T_c = 0.296 * 3.11
    #T_c = 1.15
    min_nr_corr_values = 100

    for h in h_arr:
        print(h)
        curr_sim_path = simulation_path + f"{h}/"
        for eta in eta_arr:
            quench = quench_measurement(J_para, J_perp, h, eta, p, dt, filepath, curr_sim_path + f"Damping/Quench/{eta}",
                                        quench_exec_file, runfile_quench, T_c, nr_GPUS=nr_gpus, size_max=max_size,
                                        min_nr_sites=min_nr_sites, max_nr_steps=max_nr_quench_steps,
                                        para_nr=para_nr_quench, tau_max=max_tau, nr_sites=nr_sites,
                                        min_nr_systems=min_nr_systems, min_nr_corr_values=min_nr_corr_values,
                                        Ly_Lx=Ly_Lx, equil_time=equil_time_start, equil_time_end=equil_time_end,
                                        project=project, gamma=gamma)
            quench.run()

        # the good thing is, both of the simulation implement pickup capabilities so
        # I dont need to worry to much that my computer loses connection and stuff (which it will)


if __name__ == "__main__":
    main()