from Suite import *
import numpy as np

def main():
    # I want to do multiple Quenches of which I do not necissarily know the critical
    # temperature and I want to vary the h
    # If we choose our old values still, the h should go up to 30 which would be
    # the relation of J_parallel and h in the real system
    nr_gpus = 6
    # we somehow need the relevant parameters
    # The model defining parameters are J_perp J_para h eta
    # the simulation defining parameters are dt
    J_para = -3.11
    J_perp = -0.1
    h = 0.5
    eta = 1.5
    dt = 0.01

    # filepath = "/home/weitze73/Documents/Master-Arbeit/Code/Master-Arbeit/CudaProject"
    filepath = "/home/andi/Studium/Code/Master-Arbeit/CudaProject"
    simulation_path = "../../Generated content/Silicon/Subsystems/Suite/L/"

    Tc_exec_file = "AutoCumulant.cu"

    # Tc parameters
    max_size_Tc = 80
    min_size_Tc = 48
    nr_sizes_Tc = 2
    nr_Ts = 3
    # We use small equilibration errors since we want to have really accurate
    equil_error = 0.002
    min_equil_error = 0.0005
    max_rel_intersection_error = 0.01       # is this to small or fine?
    # since we need quantitative exact values, we should know T_c beforehand
    min_T = 0.94
    max_T = 0.96
    # we should start at another parameter file nr because yeah
    para_nr = 120
    # what L-pairs do we want to check? we always just use a pair for now?
    # We could think about not using Ly / Lx = 1 / 8, then it would be easier
    # to use intermediate values and we could start at 4 or so.
    # With OBC we actually dont need to use powers of 2 to calculate U_L
    # We also dont need to use twice the size for the larger partner, but
    # I feel like this would work better and be mor precise.
    small_Ls = [8, 12, 16, 20, 24, 28, 32]
    large_Ls = [24, 28, 32, 36, 40, 44, 48]

    # I think this is good? almost every size is reused and we always have twice the size

    small_Ls = [8,  12, 16, 24, 32, 48, 64]
    large_Ls = [16, 24, 32, 48, 64, 96, 128]
    Ly_Lx = 1 / 2

    for (L_small, L_large) in zip(small_Ls, large_Ls):
        # we dont even need another path? all in one path
        #curr_sim_path = simulation_path + f"{h}/"

        # Run Tc Sim:
        Tc_sim = crit_temp_measurement(J_para, J_perp, h, eta, dt, filepath, simulation_path, Tc_exec_file, nr_GPUS=nr_gpus,
                                    size_min=L_small, size_max=L_large, Ly_Lx=Ly_Lx, nr_sizes=nr_sizes_Tc, nr_Ts=nr_Ts, T_min=min_T, T_max=max_T,
                                    equil_error=equil_error, min_equil_error=min_equil_error, intersection_error=max_rel_intersection_error, para_nr=para_nr)
        T_c, T_c_error = Tc_sim.routine()
        # to make use of the calculations we did for the larger sizes which become the smaller sizes, we need to implement a better
        # pickup capability of the Tc calculation


if __name__ == "__main__":
    main()