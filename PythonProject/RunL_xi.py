from Suite import *
import numpy as np

def main():
    # I want to do multiple Quenches of which I do not necissarily know the critical
    # temperature and I want to vary the h
    # If we choose our old values still, the h should go up to 30 which would be
    # the relation of J_parallel and h in the real system
    #h_arr = np.logspace(0.857840941039747, np.log10(30), 2)     # maybe logarithmic?
    #h_arr = [7.208434242404265]
    #h_arr = [0.2, 3]
    #h_arr = np.array([1.7320508075688776])
    h_arr = np.array([10])
    h_arr = [0.1, 0.2, 0.4161791450287818, 1.7320508075688776, 3, 7.208434242404265, 10]
    h_arr = [0.001]
    h_arr = np.array([0.4161791450287818])
    nr_gpus = 10
    # we somehow need the relevant parameters
    # The model defining parameters are J_perp J_para h eta
    # the simulation defining parameters are dt
    J_para = -3.11
    J_perp = -0.1
    p = 2.54
    eta_arr = [1]
    #eta_arr = [0.01, 0.05]
    dt = 0.01

    filepath = "/home/weitze73/Documents/Master-Arbeit/Code/Master-Arbeit/CudaProject"
    filepath = "/home/andi/Studium/Code/Master-Arbeit/CudaProject"
    simulation_path = "../../Generated content/Silicon/Subsystems/Suite/L_xi/Check-OBC/"

    Tc_exec_file = "AutoAmplitude.cu"
    quench_exec_file = "AutoQuench.cu"

    # Tc parameters
    max_size_Tc = 160
    min_size_Tc = 120
    nr_sizes_Tc = 2
    nr_Ts = 3
    para_nr_Tc = int(input("para nr, please take seriously:"))
    min_val_nr = 2000
    # We use relatively large equilibration errors since for the quenches we only need a
    # rough estimate of the transition temperature
    # for future use we could extend the pickup of the Tc measurement to work with
    # any previous measurements, not only the the ones the coincide with the current one
    equil_error = 0.025
    # We add the moving factor because I think that the point at 1.2 is not equilibrated
    moving_factor = 0.005
    min_equil_error = 0.005
    max_rel_intersection_error = 0.025

    # sample parameters to get a fee
    T_min = 0.8
    T_max = 1.1
    nr_Ts = 1

    for h in h_arr:
        curr_sim_path = simulation_path + f"{h}/"

        # Run Tc Sim:
        eta = np.mean(eta_arr)          # Which eta is the fastest?
        # Tc_sim = efficient_crit_temp_measurement_corr(J_para, J_perp, h, eta, p, dt, filepath,
        #                                          curr_sim_path + "Tc", Tc_exec_file, nr_GPUS=nr_gpus,
        #                                          size_min=min_size_Tc, size_max=max_size_Tc, equil_error=equil_error,
        #                                          min_equil_error=min_equil_error,
        #                                          intersection_error=max_rel_intersection_error,
        #                                          max_moving_factor=moving_factor,
        #                                          para_nr=para_nr_Tc, min_val_nr=min_cum_nr, value_name="xix", random_init=1)
        Tc_sim = crit_temp_measurement_corr(J_para, J_perp, h, eta, p, dt, filepath,
                                                 curr_sim_path + "Tc", Tc_exec_file, nr_GPUS=nr_gpus,
                                                 size_min=min_size_Tc, size_max=max_size_Tc, nr_sizes=2, equil_error=equil_error,
                                                 min_equil_error=min_equil_error, intersection_error=max_rel_intersection_error,
                                                 max_moving_factor=moving_factor, para_nr=para_nr_Tc, min_val_nr=min_val_nr,
                                                 T_min=T_min, T_max=T_max, nr_Ts=nr_Ts,
                                                 value_name="xix", random_init=0)
        T_c, T_c_error = Tc_sim.routine()
        # the good thing is, both of the simulation implement pickup capabilities so
        # I dont need to worry to much that my computer loses connection and stuff (which it will)


if __name__ == "__main__":
    main()