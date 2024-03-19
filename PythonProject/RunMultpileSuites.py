from Suite import *
import numpy as np

def main():
    # I want to do multiple Quenches of which I do not necissarily know the critical
    # temperature and I want to vary the h
    # If we choose our old values still, the h should go up to 30 which would be
    # the relation of J_parallel and h in the real system
    #h_arr = np.logspace(0.857840941039747, np.log10(30), 2)     # maybe logarithmic?
    h_arr = np.array([0.4161791450287818])
    #h_arr = [7.208434242404265]
    #h_arr = [0.2, 3]
    #h_arr = np.array([1.7320508075688776])
    h_arr = np.array([10])
    h_arr = [0.1, 0.2, 0.4161791450287818, 1.7320508075688776, 3, 7.208434242404265, 10]
    h_arr = [2500, 10000, 20000]
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
    simulation_path = "../../Generated content/Final/Quenches/"

    Tc_exec_file = "AutoCumulant.cu"
    quench_exec_file = "AutoQuench.cu"
    runfile_Tc = "run_cuda_gpu_a100_low.sh"
    runfile_quench = "run_cuda_gpu_a100_low.sh"

    # Tc parameters
    max_size_Tc = 96
    min_size_Tc = 48
    nr_sizes_Tc = 2
    nr_Ts = 3
    para_nr_Tc = int(input("para nr, please take seriously:"))
    file_ending = "mag"
    value_name = "U_L"
    val_write_density = 1 / 5000
    min_mag_nr = 2500
    process_file_func = process_new_mag_file_to_U_L
    equil_error = 0.003
    equil_cutoff = 0.5
    # We add the moving factor because I think that the point at 1.2 is not equilibrated
    moving_factor = 0.005
    min_equil_error = 0.005
    max_rel_intersection_error = 0.025

    # Quench parameters
    max_size = 8200
    min_nr_sites = 1e6
    max_nr_quench_steps = 4e7
    max_tau = 17000
    para_nr_quench = 230

    for h in h_arr:
        curr_sim_path = simulation_path + f"{h}/"

        # Run Tc Sim:
        eta = np.mean(eta_arr)          # Which eta is the fastest?
        #Tc_sim = efficient_crit_temp_measurement(J_para, J_perp, h, eta, p, dt, filepath,
        #                                         curr_sim_path + "Tc", Tc_exec_file, nr_GPUS=nr_gpus,
        #                                         size_min=min_size_Tc, size_max=max_size_Tc, equil_error=equil_error,
        #                                         min_equil_error=min_equil_error,
        #                                         intersection_error=max_rel_intersection_error,
        #                                         max_moving_factor=moving_factor,
        #                                         para_nr=para_nr_Tc, min_val_nr=min_mag_nr, value_name="U_L")
        Tc_sim = efficient_crit_temp_measurement(J_para, J_perp, h, eta, p, dt, filepath,
                                                 curr_sim_path + "Tc", Tc_exec_file, runfile=runfile_Tc, nr_GPUS=nr_gpus,
                                                 size_min=min_size_Tc, size_max=max_size_Tc, equil_error=equil_error,
                                                min_equil_error=min_equil_error,
                                                 intersection_error=max_rel_intersection_error, para_nr=para_nr_Tc,
                                                 min_val_nr=min_mag_nr, file_ending=file_ending, value_name=value_name,
                                                 process_file_func=process_file_func, val_write_density=val_write_density,
                                                 equil_cutoff=equil_cutoff)
        #T_c, T_c_error = Tc_sim.routine()
        T_c, T_c_error = Tc_sim.routine()
        # We could in principle run the quenches in parallel, but that would
        # require some work on my end

        # Run Quench
        # Quench can be influenced by eta
        for eta in eta_arr:
            quench = quench_measurement(J_para, J_perp, h, eta, p, dt, filepath, curr_sim_path + f"Quench/{eta}",
                                        quench_exec_file, runfile_quench, T_c, nr_GPUS=nr_gpus, size_max=max_size,
                                        min_nr_sites=min_nr_sites, max_nr_steps=max_nr_quench_steps,
                                        para_nr=para_nr_quench, tau_max=max_tau)
            quench.run()

        # the good thing is, both of the simulation implement pickup capabilities so
        # I dont need to worry to much that my computer loses connection and stuff (which it will)


if __name__ == "__main__":
    main()