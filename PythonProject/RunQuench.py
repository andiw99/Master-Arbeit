from Suite import *
import numpy as np

def main():
    h = 1
    nr_gpus = 10

    J_para = -10
    J_perp = -0.1

    p = 2.5
    eta = 1
    dt = 0.01

    filepath = "/home/andi/Studium/Code/Master-Arbeit/CudaProject"
    simulation_path = "../../Generated content/Final/Quenches-old-OBC2/"

    quench_exec_file = "AutoQuenchOBC.cu"
    runfile_quench = "run_cuda.sh"

    # Quench parameters
    max_size = 2100
    min_nr_sites = 4e6
    max_nr_quench_steps = 1e6
    nr_sites = 4e6
    max_tau = 10000
    min_nr_systems = 20
    para_nr_quench = int(input("please just change the parameter nubmer :("))

    #T_c = 21351     # maybe this T_c is to low?
    #T_c = 0.1975 * 10
    T_c = 2.325
    min_nr_corr_values = 100

    quench = quench_measurement(J_para, J_perp, h, eta, p, dt, filepath, simulation_path + f"/Quench/",
                                quench_exec_file, runfile_quench, T_c, nr_GPUS=nr_gpus, size_max=max_size,
                                min_nr_sites=min_nr_sites, max_nr_steps=max_nr_quench_steps,
                                para_nr=para_nr_quench, tau_max=max_tau, nr_sites=nr_sites,
                                min_nr_systems=min_nr_systems, min_nr_corr_values=min_nr_corr_values)
    quench.run()




if __name__ == "__main__":
    main()