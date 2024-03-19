from Suite import *
import numpy as np

def main():
    # I want to do multiple Quenches of which I do not necissarily know the critical
    # temperature and I want to vary the h
    # If we choose our old values still, the h should go up to 30 which would be
    # the relation of J_parallel and h in the real system
    #h_arr = np.logspace(0.857840941039747, np.log10(30), 2)     # maybe logarithmic?
    h_arr = [5200] #, 1000, 10000, 20000]
    nr_gpus = 10
    # we somehow need the relevant parameters
    # The model defining parameters are J_perp J_para h eta
    # the simulation defining parameters are dt
    J_para = -130000
    J_perp = -1300
    p = 2.5
    eta_arr = [1]
    #eta_arr = [0.01, 0.05]
    dt = 1e-5

    filepath = "/home/weitze73/Documents/Master-Arbeit/Code/Master-Arbeit/CudaProject"
    filepath = "/home/andi/Studium/Code/Master-Arbeit/CudaProject"
    simulation_path = "../../Generated content/Final/z-measurement/"

    z_exec_file = "AutoZ.cu"
    z_test_exec_file = "AutoCumulant.cu"

    runfile_z = "run_cuda.sh"

    # z parameters
    para_nr_z = int(input("parameter number ..."))
    size_min_z = 64
    size_max_z = 256
    nr_sizes = 3
    z_min_nr_sites = 1e6
    z_min_nr_systems = 500
    z_equil_error = 0.002

    # mag stuff
    file_ending = "mag"
    value_name = "m"
    process_file_func = process_new_mag_file_to_U_L
    # We use relatively large equilibration errors since for the quenches we only need a
    # rough estimate of the transition temperature
    # for future use we could extend the pickup of the Tc measurement to work with
    # any previous measurements, not only the the ones the coincide with the current one
    test_min_val_nr = 200
    val_write_density = 1 / 100
    val_write_density_test = 1 / 10000
    moving_factor = 0.02
    min_equil_error = 0.01
    equil_cutoff_Tc = 0.5
    max_rel_intersection_error = 0.03

    T_c = 21700

    for h in h_arr:
        print(h)
        curr_sim_path = simulation_path + f"{h}/"
        for eta in eta_arr:
            z_measure = z_measurement(J_para, J_perp, h, eta, p, dt, filepath, curr_sim_path + "z",
                                      z_exec_file, z_test_exec_file, runfile_z, T_c, nr_GPUS=nr_gpus, size_min=size_min_z, size_max=size_max_z,
                                      nr_sizes=nr_sizes, min_nr_sites=z_min_nr_sites, min_nr_systems=z_min_nr_systems,
                                      equil_error=z_equil_error, para_nr=para_nr_z, test_min_val_nr=test_min_val_nr,
                                      val_write_density=val_write_density, test_val_write_density=val_write_density_test,
                                      file_ending=file_ending, value_name=value_name)
            z_measure.run()

        # the good thing is, both of the simulation implement pickup capabilities so
        # I dont need to worry to much that my computer loses connection and stuff (which it will)


if __name__ == "__main__":
    main()