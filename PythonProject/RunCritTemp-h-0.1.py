from Suite import *

def main():
    # okay what is the first thing we need to do?
    # we need parameters like the number of gpus we are able to use
    nr_gpus = 20
    # we somehow need the relevant parameters
    # The model defining parameters are J_perp J_para h eta
    # the simulation defining parameters are dt
    #J_para = -130000
    J_para = -10
    #J_perp = -1300
    J_perp = -0.1
    h = 0.2
    #h = 0.5
    eta = 1
    p = 2.5
    dt = 0.01
    max_size_Tc = 64
    min_size_Tc = 128
    nr_sizes_Tc = 2
    nr_Ts = 10
    cum_error = 0.0015
    equil_cutoff_Tc = 0.1
    value_name = "U_L"
    file_ending = "mag"
    process_file_func = recalculate_mag_file_to_U_L
    value_write_density = 0.01
    nr_sites = 4e6
    Ly_Lx = 1 / 8


    random_init = 0.0
    #filepath = "/home/weitze73/Documents/Master-Arbeit/Code/Master-Arbeit/CudaProject"
    filepath = "/home/andi/Studium/Code/Master-Arbeit/CudaProject"
    simulation_path = "../../Generated content/Final/CriticalTemperature/J_J=100-old-OBC-h-0.2/"

    Tc_exec_file = "AutoCumulantOBC.cu"
    quench_exec_file = "AutoQuench.cu"
    amplitude_exec_file = "AutoAmplitude.cu"
    z_exec_file = "AutoZ.cu"        # what dow we need here? Maybe different files depending if we are doing the test measurement or the real one?
    z_test_exec_file = "AutoCumulant.cu"
    # for the real measurements we have fixed end times and we extract the cumulant a fixed number of times
    # for the test measurement we extract a density of cumulants, calculate the error and run until we are equilibrated.
    # In both cases we start in a high temperature phase. For the testmeasurement we can actually just use the amplitude file?
    # for the other ones I think we need a new file. Okay we can maybe use the amplitude file, but it observes the correlation length
    # We want to observe the binder cumulant. But for the equilibration it should not make to much difference. But tbh i also
    # want to work with the new error
    runfile = "run_cuda_gpu_a100_low.sh"
    runfile = "run_cuda.sh"

    # T- parameters?
    max_rel_intersection_error = 0.01
    min_cum_nr = 5000
    moving_factor = 0.001
    T_min = None
    T_max = None



    para_nr = int(input("parameter number.."))
    sim = efficient_crit_temp_measurement(J_para, J_perp, h, eta, p, dt, filepath, simulation_path + "Tc", Tc_exec_file,
                                            runfile, nr_GPUS=nr_gpus, T_min=T_min, T_max=T_max, Ly_Lx=Ly_Lx,
                                            size_min=min_size_Tc, size_max=max_size_Tc, equil_error=cum_error,
                                            min_equil_error=cum_error, intersection_error=max_rel_intersection_error,
                                            max_moving_factor=moving_factor, para_nr=para_nr, random_init=random_init,
                                            min_val_nr=min_cum_nr, file_ending=file_ending, value_name=value_name,
                                            process_file_func=process_file_func, val_write_density=value_write_density,
                                            equil_cutoff=equil_cutoff_Tc, nr_sites=nr_sites)
    T_c, T_c_error = sim.routine()


if __name__ == '__main__':
    main()