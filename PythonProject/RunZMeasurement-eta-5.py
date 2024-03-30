from Suite import *
import numpy as np

def main():
    h_arr = [1]
    nr_gpus = 10

    J_para = -10
    J_perp = -0.1


    p = 2.5
    eta_arr = [5]

    dt = 0.01

    filepath = "/home/weitze73/Documents/Master-Arbeit/Code/Master-Arbeit/CudaProject"
    filepath = "/home/andi/Studium/Code/Master-Arbeit/CudaProject"
    #simulation_path = "../../Generated content/Final/z-measurement-small/eta=5/"
    simulation_path = "../../../Generated content without sync/Final/z-measurement-small/eta=5/"
    z_exec_file = "AutoZ.cu"
    z_test_exec_file = "AutoCumulantOBC.cu"
    runfile_z = "run_cuda_gpu_a100_low.sh"
    #runfile_z =  "run_cuda_casus_low.sh"

    # z parameters
    para_nr_z = int(input("parameter number ..."))
    size_min_z = 48
    size_max_z = 144
    z_test_size = 24
    nr_sizes = 4
    z_min_nr_sites = 10e6
    z_min_nr_systems = 50000
    z_equil_error = 0.004
    fold=40

    # mag stuff
    file_ending = "mag"
    value_name = "m"
    process_file_func = process_new_mag_file_to_U_L
    # We use relatively large equilibration errors since for the quenches we only need a
    # rough estimate of the transition temperature
    # for future use we could extend the pickup of the Tc measurement to work with
    # any previous measurements, not only the the ones the coincide with the current one
    # test_min_val_nr = 1000
    # val_write_density = 1 / 1000
    # val_write_density_test = 1 / 1000
    test_min_val_nr = 200
    val_write_density = 1 / 25
    val_write_density_test = 1 / 25
    Ly_Lx = 1 / 12

    variation_error_rate = 0.01
    nr_sites = 10e6      # we use large systems because I think the cluster doesnt like it if we start very many runs
    T_c = 21700
    T_c = 0.1965 * 10
    notest=True

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
                                      nr_sites=nr_sites, test_size=z_test_size, fold=fold, Ly_Lx=Ly_Lx, notest=notest)
            z_measure.run()

        # the good thing is, both of the simulation implement pickup capabilities so
        # I dont need to worry to much that my computer loses connection and stuff (which it will)


if __name__ == "__main__":
    main()