import numpy as np
from scipy.optimize import fsolve
from itertools import product

def T_c_XY(T, J_parallel, J_perp):
    return 2 * T / J_parallel * np.log(2 * T / J_perp) - 1

def T_c_est(J_para, J_perp, h):
    # we need to solve the transcendet equation of the XY model
    # scipy f_solve needs a guess of the critical temperature which is difficult but lets try with the mean of the Js
    T_c_est = (J_para + J_perp) / 2
    T_c_est = fsolve(T_c_XY, T_c_est, args=(J_para, J_perp))
    return T_c_est


def start_cum_measurement(J_para, J_perp, h, eta, dt, filepath, measurement_path, nr_GPUS=6, nr_Ts=5, size_min=48,
                          size_max=80, nr_sizes=3, max_steps=1e9, nr_sites=5e6, Ly_Lx = 1/8):
    # this somehow needs the parameters, where do we put them? In a file? On the moon? User input?
    T_min = T_c_est(J_para, J_perp, h)
    # TODO really crude approximation of the maximum T, I really should think about something better
    T_max = T_min + (h/J_perp) * T_min
    # We use nr_Ts datapoints
    T_arr = np.linspace(T_min, T_max, nr_Ts)
    sizes = np.linspace(size_min, size_max, nr_sizes, endpoint=True, dtype=np.int32)
    max_time = dt * max_steps

    for i, (T, size) in enumerate(product(T_arr, sizes)):
        # We now need to construct the parameterfile with the appropriate temperature and size
        # Is it okay if we construct all files in the beginning and deal with the threading of the gpus later?
        # to construct the para set we need to know how many subsystems we should initialize
        nr_subsystems = int(nr_sites / (size ** 2 * Ly_Lx))
        with open(filepath + "/para_set_100" + str(i) + '.txt', 'w') as f:
            f.write(measurement_path)
            f.write(f"end_time, {max_time} \n"
                    f"dt, {dt} \n"
                    f"J, {J_para} \n"
                    f"Jy, {J_perp} \n"
                    f"alpha, {h} \n"
                    f"eta, {eta} \n"
                    f"nr_saves, 4 \n"
                    f"nr_repeat, 0 \n"
                    f"min_temp, {T} \n"
                    f"max_temp, {T} \n"
                    f"nr_runs, 0.0 \n"
                    f"random_init, 0.0 \n"
                    f"curand_random, 1 \n"
                    f"subsystem_min_Lx, {size} \n"
                    f"subsystem_max_Lx, {size} \n"
                    f"nr_subsystem_sizes, 0 \n"
                    f"")


def main():
    # okay what is the first thing we need to do?
    # we need parameters like the number of gpus we are able to use
    nr_gpus = 6
    # we somehow need the relevant parameters
    # The model defining parameters are J_perp J_para h eta
    # the simulation defining parameters are dt
    J_para = 110000
    J_perp = 3500
    h = 500
    # I honestly have no idea on how to account h, that is really a problem
    # the Scanned interval

    T_c_est(J_para, J_perp, h)



if __name__ == '__main__':
    main()