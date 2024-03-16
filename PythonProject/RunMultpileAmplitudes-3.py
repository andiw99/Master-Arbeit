from Suite import *
import numpy as np

def main():
    # I want to do multiple Quenches of which I do not necissarily know the critical
    # temperature and I want to vary the h
    # If we choose our old values still, the h should go up to 30 which would be
    # the relation of J_parallel and h in the real system
    #    h_arr = np.logspace(-1, np.log10(30), 5)     # maybe logarithmic?
    h_arr = np.array([0.4161791450287818])
    nr_gpus = 6
    # we somehow need the relevant parameters
    # The model defining parameters are J_perp J_para h eta
    # the simulation defining parameters are dt
    #J_para = -120000
    J_para = -3
    #J_perp = -2000
    J_perp = -0.1
    #Ly_Lx = 1 / 16
    Ly_Lx = 1
    p = 2.54
    eta = 1.5
    dt = 0.01

    #filepath = "/home/weitze73/Documents/Master-Arbeit/Code/Master-Arbeit/CudaProject"
    filepath = "/home/andi/Studium/Code/Master-Arbeit/CudaProject"
    simulation_path = "../../Generated content/Silicon/Subsystems/Suite/h/Large Jx/Jx=3-Lx_Ly=1/"

    Tc_exec_file = "AutoCumulant.cu"
    amplitude_exec_file = "AutoAmplitude.cu"
    runfile = "run_cuda_gpu_a100_low.sh"

    # Tc parameters
    max_size_Tc = 192
    min_size_Tc = 64
    file_ending = "mag"
    value_name = "m"
    process_file_func = process_new_mag_file_to_U_L
    nr_sizes_Tc = 2
    nr_Ts = 3
    para_nr_Tc = int(input("para nr, please take seriously:"))
    # We use relatively large equilibration errors since for the quenches we only need a
    # rough estimate of the transition temperature
    # for future use we could extend the pickup of the Tc measurement to work with
    # any previous measurements, not only the the ones the coincide with the current one
    min_val_nr = 2000
    equil_error = 0.04
    val_write_density = 1 / 1000            # otherwise the files become to large?
    moving_factor = 0.02
    min_equil_error = 0.01
    max_rel_intersection_error = 0.01

    # Amplitude parameters
    amplitude_size = 4048
    equil_error_amplitude = 0.05
    equil_cutoff = 0.01
    para_nr_ampl = int(input("para nr amplitude, please take seriously:"))
    observed_direction = int(input("observed direction :"))
    #T_min_fraction = 0.0025
    T_min_fraction = -0.0025
    T_range_fraction = 0.03
    nr_Ts = 4
    T_c = 0.903
    #T_c = 1.7268

    amplitude_sizes = [512] #[2048, 1024]
    #amplitude_sizes = [512, 256, 128]
    T_ranges = [0.00]#, 0.02, 0.03]
    nr_Ts_per_range = 2
    last_T = None

    for size, T_up in zip(amplitude_sizes, T_ranges):
        h = h_arr[0]
        if last_T:
            T_min = T_min_fraction + last_T
        else:
            T_min = T_min_fraction
        ampl = amplitude_measurement(J_para, J_perp, h, eta, p, dt, filepath, simulation_path + f"/Amplitude",
                                     amplitude_exec_file, runfile, T_c, nr_GPUS=nr_gpus, size=size,
                                     equil_error=equil_error_amplitude, equil_cutoff=equil_cutoff, para_nr=para_nr_ampl,
                                     T_min_fraction=T_min, T_range_fraction=T_up, nr_Ts=nr_Ts_per_range,
                                     Ly_Lx=Ly_Lx, observed_direction=observed_direction)
        ampl.run()
        last_T = T_up

    for h in h_arr:
        curr_sim_path = simulation_path + f"{h}/"

        # Run Tc Sim:
        # Tc_sim = efficient_crit_temp_measurement(J_para, J_perp, h, eta, p, dt, filepath, curr_sim_path + "Tc",
        #                                          Tc_exec_file, runfile, nr_GPUS=nr_gpus,
        #                             size_min=min_size_Tc, size_max=max_size_Tc, equil_error=equil_error,
        #                                          min_equil_error=min_equil_error, intersection_error=max_rel_intersection_error,
        #                                          max_moving_factor=moving_factor, para_nr=para_nr_Tc, Ly_Lx=Ly_Lx,
        #                                          min_val_nr=min_val_nr, file_ending=file_ending, value_name=value_name,
        #                                          process_file_func=process_file_func, val_write_density=val_write_density)
        # T_c, T_c_error = Tc_sim.routine()
        # We could in principle run the quenches in parallel, but that would
        # require some work on my end

        ampl = amplitude_measurement(J_para, J_perp, h, eta, p, dt, filepath, curr_sim_path + "Amplitude",
                                     amplitude_exec_file, runfile, T_c, nr_GPUS=nr_gpus, size=amplitude_size,
                                     equil_error=equil_error_amplitude, equil_cutoff=equil_cutoff, para_nr=para_nr_ampl,
                                     T_min_fraction=T_min_fraction, T_range_fraction=T_range_fraction, nr_Ts=nr_Ts,
                                     Ly_Lx=Ly_Lx, observed_direction=observed_direction)
        ampl.run()

        # the good thing is, both of the simulation implement pickup capabilities so
        # I dont need to worry to much that my computer loses connection and stuff (which it will)


if __name__ == "__main__":
    main()