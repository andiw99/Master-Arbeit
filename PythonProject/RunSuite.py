from Suite import *

def main():
    # okay what is the first thing we need to do?
    # we need parameters like the number of gpus we are able to use
    nr_gpus = 1
    # we somehow need the relevant parameters
    # The model defining parameters are J_perp J_para h eta
    # the simulation defining parameters are dt
    #J_para = -130000
    J_para = -10
    #J_perp = -1300
    J_perp = -0.1
    # h = 0.282727
    h = 1
    #h = 0.5
    eta = 1
    p = 2.5
    dt = 0.00001
    dt = 0.01
    max_size_Tc = 32
    min_size_Tc = 128
    nr_sizes_Tc = 2
    nr_Ts = 5
    cum_error = 0.001
    equil_cutoff_Tc = 0.1
    value_name = "U_L"
    file_ending = "mag"
    process_file_func = recalculate_mag_file_to_U_L
    value_write_density = 0.01
    nr_sites = 4e6
    Ly_Lx = 1 / 8


    random_init = 0.0
    project = "MinimalCudaProject"
    filepath = f"/home/andi/Studium/Code/Master-Arbeit/{project}"
    filepath = f"/home/weitze73/Documents/Master-Arbeit/Code/Master-Arbeit/{project}"
    simulation_path = "../../Generated content/Paper content/Binder intersection/m vectorial"


    Tc_exec_file = "AutoCumulantVectorialOBC.cu"
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
    runfile = "run_cuda.sh"
    runfile = "run_cuda_gpu_a100_low_minimal.sh"

    # T- parameters?
    max_rel_intersection_error = 0.5
    min_cum_nr = 750
    moving_factor = 0.001
    # T_min = 29071.961123
    # T_max = 31494.624550
    # for the rough J_J = 100, the minimum should be at T / J = 0.15
    # Tc is at T / J = 0.197 so almost at 0.2?
    # J_J = 31 detailed settings:
    # T_min = 0.267 * np.abs(J_para)
    # T_max = 0.281 * np.abs(J_para)
    # rough measurement J_J = 100
    # T_min = 0.155 * np.abs(J_para)
    # T_max = 0.245 * np.abs(J_para)
    T_min = 0.17 * np.abs(J_para)
    T_max = 0.22 * np.abs(J_para)


    #T_min = 0.83601154
    #T_max = 0.8701344599999999

    # Quench parameters
    max_size = 1024
    min_nr_sites = 1e6


    # Amplitude parameters
    amplitude_size = 1024
    equil_error_ampl = 0.05
    equil_cutoff = 0.01

    # z parameters
    size_min_z = 64
    size_max_z = 256
    nr_sizes = 3
    z_min_nr_sites = 1e6
    z_min_nr_systems = 500
    z_equil_error = 0.005

    # Enter which calculations are supposed to run here
    measurements = {
        "Tc": True,
        "efficient Tc": False,
        "Quench": False,
        "Amplitude": False,
        "z": False,
    }

    # I honestly have no idea on how to account h, that is really a problem
    # the Scanned interval
    if measurements["Tc"]:
        para_nr = int(input("parameter number.."))
        sim = crit_temp_measurement(J_para, J_perp, h, eta, p, dt, filepath, simulation_path + "Tc", Tc_exec_file, nr_GPUS=nr_gpus,
                                    size_min=min_size_Tc, size_max=max_size_Tc, nr_sizes=nr_sizes_Tc, nr_Ts=nr_Ts,
                                    intersection_error=max_rel_intersection_error, max_moving_factor=moving_factor,
                                    min_val_nr=min_cum_nr, equil_error=cum_error, equil_cutoff=equil_cutoff_Tc,
                                    file_ending=file_ending, value_name=value_name, process_file_func=process_file_func,
                                    val_write_density=value_write_density, runfile=runfile,
                                    T_min=T_min, T_max=T_max, para_nr=para_nr, project=project)
        T_c, T_c_error = sim.routine()
    elif measurements["efficient Tc"]:
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
    else:
        T_c = float(input("Enter critical temperature:"))
        T_c_error = 0
    if measurements["Quench"]:
        quench = quench_measurement(J_para, J_perp, h, eta, p, dt, filepath, simulation_path + "Quench", quench_exec_file, T_c, nr_GPUS=nr_gpus, size_max=max_size, min_nr_sites=min_nr_sites )
        quench.run()
    if measurements["Amplitude"]:
        ampl = amplitude_measurement(J_para, J_perp, h, eta, p, dt, filepath, simulation_path + "Amplitude",
                                     amplitude_exec_file, T_c, nr_GPUS=nr_gpus, size=amplitude_size,
                                     equil_error=equil_error_ampl, equil_cutoff=equil_cutoff)
        ampl.run()
    if measurements["z"]:
        z_measure = z_measurement(J_para, J_perp, h, eta, p, dt, filepath, simulation_path + "z",
                                        z_exec_file, z_test_exec_file, T_c, nr_GPUS=nr_gpus, size_min_z=size_min_z, size_max_z=size_max_z,
                                        nr_sizes=nr_sizes, min_nr_sites=z_min_nr_sites, min_nr_systems=z_min_nr_systems,
                                     equil_error=z_equil_error)
        z_measure.run()

if __name__ == '__main__':
    main()