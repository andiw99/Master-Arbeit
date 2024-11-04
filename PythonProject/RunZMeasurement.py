from Suite import *
import numpy as np

def main():

    h_arr = [1]
    nr_gpus = 10
    # we somehow need the relevant parameters
    # The model defining parameters are J_perp J_para h eta
    # the simulation defining parameters are dt
    J_para = -10
    J_perp = -0.1
    # J_para = -130000
    # J_perp = -1300

    p = 2.5
    eta_arr = [1]
    #eta_arr = [0.01, 0.05]
    #dt = 1e-5
    dt = 1e-4
    dt = 0.01

    project  = "MinimalCudaProject"
    filepath = f"/home/weitze73/Documents/Master-Arbeit/Code/Master-Arbeit/{project}"
    filepath = f"/home/andi/Studium/Code/Master-Arbeit/{project}"
    simulation_path = "../../Generated content/Paper content/z measurement/"

    z_exec_file = "AutoZ.cu"
    z_test_exec_file = "AutoCumulantOBC.cu"
    #runfile_z = "run_cuda_gpu_a100_low.sh"
    #runfile_z =  "run_cuda_casus_low.sh"
    runfile_z = "run_cuda.sh"
    runfile_z = "run_cuda_gpu_a100_low_minimal.sh"

    # z parameters
    para_nr_z = int(input("parameter number ..."))
    size_min_z = 48
    size_max_z = 144
    z_test_size = 24
    nr_sizes = 4
    z_min_nr_sites = 4e6
    z_min_nr_systems = 10000
    z_equil_error = 0.004
    fold=40

    # mag stuff
    file_ending = "mag"
    value_name = "m"

    test_min_val_nr = 100
    val_write_density = 1 / 10
    val_write_density_test = 1 / 10
    Ly_Lx = 1 / 8

    variation_error_rate = 0.002
    nr_sites = 10e6      # we use large systems because I think the cluster doesnt like it if we start very many runs
    T_c = 21700
    T_c = 0.1731 * 10
    T_c = 1.97

    for h in h_arr:
        print(h)
        curr_sim_path = simulation_path + f"{h}/"
        for eta in eta_arr:
            z_measure = z_measurement(J_para, J_perp, h, eta, p, dt, filepath, curr_sim_path + "z",
                                      z_exec_file, z_test_exec_file, runfile_z, T_c, nr_GPUS=nr_gpus, size_min=size_min_z, size_max=size_max_z,
                                      nr_sizes=nr_sizes, min_nr_sites=z_min_nr_sites, min_nr_systems=z_min_nr_systems,
                                      equil_error=z_equil_error, para_nr=para_nr_z, test_min_val_nr=test_min_val_nr,
                                      val_write_density=val_write_density, test_val_write_density=val_write_density_test,
                                      file_ending=file_ending, value_name=value_name, variation_error_rate=variation_error_rate,
                                      nr_sites=nr_sites, test_size=z_test_size, fold=fold, Ly_Lx=Ly_Lx, project=project)
            z_measure.run()

        # the good thing is, both of the simulation implement pickup capabilities so
        # I dont need to worry to much that my computer loses connection and stuff (which it will)


if __name__ == "__main__":
    main()