from Suite import *
import numpy as np

def main():
    # I want to do multiple Quenches of which I do not necissarily know the critical
    # temperature and I want to vary the h
    # If we choose our old values still, the h should go up to 30 which would be
    # the relation of J_parallel and h in the real system
    #    h_arr = np.logspace(-1, np.log10(30), 5)     # maybe logarithmic?
    h_arr = np.array([3300])
    nr_gpus = 6
    # we somehow need the relevant parameters
    # The model defining parameters are J_perp J_para h eta
    # the simulation defining parameters are dt
    J_para = -77500
    J_perp = -2500
    Ly_Lx = 1/8
    p = 2.54        # forgott all the time to overgive p so we will use this real quck
    #p = 2.33
    eta = 1
    dt = 1e-5

    #filepath = "/home/weitze73/Documents/Master-Arbeit/Code/Master-Arbeit/CudaProject"
    filepath = "/home/andi/Studium/Code/Master-Arbeit/CudaProject"
    simulation_path = "../../Generated content/Silicon/Subsystems/Suite/Exp/h=3300/Jx_Jy=31/"

    Tc_exec_file = "AutoCumulant.cu"
    amplitude_exec_file = "AutoAmplitude.cu"
    runfile = "run_cuda_gpu_a100_low.sh"

    # Tc parameters
    max_size_Tc = 80
    min_size_Tc = 48
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
    min_val_nr = 200
    equil_error = 0.005
    val_write_density = 1 / 10000
    moving_factor = 0.02
    min_equil_error = 0.01
    max_rel_intersection_error = 0.01

    # Amplitude parameters
    amplitude_size = 2048
    equil_error_amplitude = 0.03
    equil_cutoff = 0.01
    para_nr_ampl = 170
    T_min_fraction = 0.0075
    T_range_fraction = 0.0125
    nr_Ts = 2

    for h in h_arr:
        curr_sim_path = simulation_path + f"{h}/"

        # Run Tc Sim:
        Tc_sim = efficient_crit_temp_measurement(J_para, J_perp, h, eta, p, dt, filepath,
                                                 curr_sim_path + "Tc", Tc_exec_file, runfile, nr_GPUS=nr_gpus,
                                                 size_min=min_size_Tc, size_max=max_size_Tc, equil_error=equil_error,
                                                 min_equil_error=min_equil_error, intersection_error=max_rel_intersection_error,
                                                 max_moving_factor=moving_factor, para_nr=para_nr_Tc, Ly_Lx=Ly_Lx,
                                                 min_val_nr=min_val_nr, file_ending=file_ending, value_name=value_name,
                                                 process_file_func=process_file_func, val_write_density=val_write_density)
        T_c, T_c_error = Tc_sim.routine()
        # We could in principle run the quenches in parallel, but that would
        # require some work on my end

        ampl = amplitude_measurement(J_para, J_perp, h, eta, p, dt, filepath, curr_sim_path + "Amplitude",
                                     amplitude_exec_file, runfile, T_c, nr_GPUS=nr_gpus, size=amplitude_size,
                                     equil_error=equil_error_amplitude, equil_cutoff=equil_cutoff, para_nr=para_nr_ampl,
                                     T_min_fraction=T_min_fraction, T_range_fraction=T_range_fraction, nr_Ts=nr_Ts, Ly_Lx=Ly_Lx)
        ampl.run()

        # the good thing is, both of the simulation implement pickup capabilities so
        # I dont need to worry to much that my computer loses connection and stuff (which it will)


if __name__ == "__main__":
    main()