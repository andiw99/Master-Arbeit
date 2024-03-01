from Suite import *
import numpy as np

def main():
    # I want to do multiple Quenches of which I do not necissarily know the critical
    # temperature and I want to vary the h
    # If we choose our old values still, the h should go up to 30 which would be
    # the relation of J_parallel and h in the real system
    #h_arr = np.logspace(0.857840941039747, np.log10(30), 2)     # maybe logarithmic?
    h_arr = np.array([29.999999999999996])
    nr_gpus = 6
    # we somehow need the relevant parameters
    # The model defining parameters are J_perp J_para h eta
    # the simulation defining parameters are dt
    J_para = -3.11
    J_perp = -0.1
    p = 2.54
    eta = 1.5
    dt = 0.01

    #filepath = "/home/weitze73/Documents/Master-Arbeit/Code/Master-Arbeit/CudaProject"
    filepath = "/home/andi/Studium/Code/Master-Arbeit/CudaProject"
    simulation_path = "../../Generated content/Silicon/Subsystems/Suite/h/"

    Tc_exec_file = "AutoCumulant.cu"
    quench_exec_file = "AutoQuench.cu"

    # Tc parameters
    max_size_Tc = 80
    min_size_Tc = 48
    nr_sizes_Tc = 2
    nr_Ts = 3
    # We use relatively large equilibration errors since for the quenches we only need a
    # rough estimate of the transition temperature
    # for future use we could extend the pickup of the Tc measurement to work with
    # any previous measurements, not only the the ones the coincide with the current one
    equil_error = 0.05
    # We add the moving factor because I think that the point at 1.2 is not equilibrated
    moving_factor = 0.2
    min_equil_error = 0.01
    max_rel_intersection_error = 0.05

    # Quench parameters
    max_size = 1024
    min_nr_sites = 1e6
    max_nr_quench_steps = 4e7

    for h in h_arr:
        curr_sim_path = simulation_path + f"{h}/"

        # Run Tc Sim:
        Tc_sim = efficient_crit_temp_measurement(J_para, J_perp, h, eta, p, dt, filepath,
                                                 curr_sim_path + "Tc", Tc_exec_file, nr_GPUS=nr_gpus,
                                                 size_min=min_size_Tc, size_max=max_size_Tc, equil_error=equil_error,
                                                 min_equil_error=min_equil_error,
                                                 intersection_error=max_rel_intersection_error,
                                                 max_moving_factor=moving_factor)
        T_c, T_c_error = Tc_sim.routine()
        # We could in principle run the quenches in parallel, but that would
        # require some work on my end

        # Run Quench
        quench = quench_measurement(J_para, J_perp, h, eta, p, dt, filepath, curr_sim_path + "Quench",
                                    quench_exec_file, T_c, nr_GPUS=nr_gpus, size_max=max_size,
                                    min_nr_sites=min_nr_sites, max_nr_steps=max_nr_quench_steps)
        quench.run()

        # the good thing is, both of the simulation implement pickup capabilities so
        # I dont need to worry to much that my computer loses connection and stuff (which it will)


if __name__ == "__main__":
    main()