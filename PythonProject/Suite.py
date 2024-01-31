import matplotlib.pyplot as plt
import numpy as np
from itertools import product
import re
from fabric import Connection
import time
from FunctionsAndClasses import *
import subprocess
import pathlib
from scipy.interpolate import CubicSpline, PchipInterpolator


def extract_numbers_after_newline(input_string):
    # Define the regular expression pattern
    pattern = r'\n\s*(\d+)'
    # Use re.findall to find all matches of the pattern
    matches = re.findall(pattern, input_string)
    # Convert the matched strings to integers
    numbers = [int(match) for match in matches]
    return numbers

def extract_numbers_before_newline(input_string):
    # Define the regular expression pattern
    pattern = r'(\d+)\n'
    # Use re.findall to find all matches of the pattern
    matches = re.findall(pattern, input_string)
    # Convert the matched strings to integers
    numbers = [int(match) for match in matches]
    return numbers

def check_completed_status(number, input_string):
    # Split the input string into lines
    lines = input_string.split('\n')

    # Iterate through each line
    for line in lines:
        # Check if the number is present in the line
        if str(number) in line:
            # Check if 'COMPLETED' exists in the rest of the line
            if 'COMPLETED' in line.split(str(number), 1)[1]:
                return True  # 'COMPLETED' found
    return False  # 'COMPLETED' not found

def check_directory_structure(sizes, temperatures, directory_path):
    # Check if the directory exists
    if not os.path.exists(directory_path) or not os.path.isdir(directory_path):
        return False

    # Iterate over sizes
    for size in sizes:
        size_path = os.path.join(directory_path, str(size))

        # Check if the size folder exists
        if not os.path.exists(size_path) or not os.path.isdir(size_path):
            print("No size folders available")
            return False

        # Iterate over temperatures
        for temp in temperatures:
            temp_path = os.path.join(size_path, "{:.6f}".format(temp))      # we only write 6

            # Check if the temperature folder exists
            if not os.path.exists(temp_path) or not os.path.isdir(temp_path):
                print("No temperature folders available")
                return False

            # Check if there is at least one "*.cum" file in the temperature folder
            cum_files = [f for f in os.listdir(temp_path) if f.endswith('.cum')]
            if not cum_files:
                print("No cumulant files available")
                return False

    # If all checks pass, return True
    return True

def get_avail_simulations(sizes, temperatures, directory_path, check_function):
    """
    Checks if the directory structure is valid based on a custom check function.

    Parameters:
    - sizes: List of sizes.
    - temperatures: List of temperatures.
    - directory_path: Path to the main directory.
    - check_function: A custom function to check the validity of each temperature folder.

    Returns:
    - List of (size, temp) pairs for which a valid folder exists.
    """
    valid_folders = []

    # Iterate over sizes
    for size in sizes:
        size_path = os.path.join(directory_path, str(size))

        # Check if the size folder exists
        if os.path.exists(size_path) and os.path.isdir(size_path):
            # Iterate over temperatures
            for temp in temperatures:
                temp_path = os.path.join(size_path, str(temp))

                # Check if the temperature folder exists and is valid based on the custom function
                if os.path.exists(temp_path) and os.path.isdir(temp_path) and check_function(temp_path):
                    valid_folders.append((size, temp))

    return valid_folders

class crit_temp_measurement():
    def __init__(self, J_para, J_perp, h, eta, dt, filepath, simulation_path, nr_GPUS=6, nr_Ts=5, size_min=48,
                          size_max=80, nr_sizes=3, max_steps=1e9, nr_sites=5e5, Ly_Lx = 1/8, equil_error=0.003, intersection_error=0.002):
        self.J_para = J_para
        self.J_perp = J_perp
        self.h = h
        self.eta = eta
        self.dt = dt
        self.filepath = filepath
        self.simulation_path = simulation_path
        self.nr_GPUS = nr_GPUS
        self.nr_Ts = nr_Ts
        self.size_min = size_min
        self.size_max = size_max
        self.nr_sizes = nr_sizes
        self.max_steps = max_steps
        self.nr_sites = nr_sites
        self.Ly_Lx = Ly_Lx

        self.T_arr = np.array([])
        self.sizes = np.array([])
        self.max_time = 0
        self.total_runs = 0
        self.running_jobs = set()
        self.completed_jobs = set()
        self.connection = None

        self.all_T_arr = np.array([])       # Bookkeeping for all the temperatures we have simulated in this setting
        self.max_rel_intersec_error = intersection_error  # standard maximum error of 2%
        self.equil_error = equil_error           # standard equilibration error for the U_L runs
        self.maximum_iterations = 4
        self.iteration_nr = 0
        self.repeat = False             # variable that is set to true if we have to repeat a simulation
        self.min_cum_nr = 100
    def init(self):
        # this somehow needs the parameters, where do we put them? In a file? On the moon? User input?
        T_min = T_c_est(np.abs(self.J_para), np.abs(self.J_perp), self.h)[0]
        # TODO really crude approximation of the maximum T, I really should think about something better
        print(f"T_min = {T_min}\n"
              f"J_perp = {self.J_perp}\n"
              f"h={self.h}")
        T_max = T_min + (self.h / (5 * np.abs(self.J_perp))) * T_min

        # We use nr_Ts datapoints
        self.T_arr = np.linspace(T_min, T_max, self.nr_Ts)
        self.sizes = np.linspace(self.size_min, self.size_max, self.nr_sizes, endpoint=True,
                            dtype=np.int32)
        print(f"Initializing Simulation with T_min = {T_min} and T_max = {T_max}")
        self.max_time = self.dt * self.max_steps
        self.total_runs = self.nr_sizes * self.nr_Ts  # every temp size combo will be a seperate simulation
    def routine(self, walltime='24:00:00', file='SubsystemRelaxation.cu', folder='simulations', user='weitze73', wait=40):
        """
        the outer routine for the T_c calculation
        :param walltime: walltime per job+3,
        :param file: file to be compiled, we may want to have a seperate file for the autonomous exection
        :param folder: folder for the executables
        :param wait: sleep for wait seconds before checking again if a job is finished        
        :return: -
        """
        self.init()                 # initializing determines the T range
        T_c, T_c_error = self.iteration(file, folder, user, wait, walltime)  # recursive fuction that returns T_c
        self.conclude()
        return T_c, T_c_error

    def conclude(self):
        # This function should plot stuff etc. keep it simple at this point
        # copied just the cumulanttimeaverage script!
        # has to be rewritten...
        threshold = 0.1
        results = {}
        for size_folder in os.listdir(self.simulation_path):
            if (size_folder[0] != ".") & (size_folder != "plots"):
                size_folder_path = os.path.join(self.simulation_path, size_folder)
                if os.path.isdir(size_folder_path):
                    size_result = process_size_folder(size_folder_path,
                                                      threshold, selected_temperatures=self.T_arr)
                    results[int(size_folder)] = size_result

        # interpolate and minimize is deprecated, we use the technique we also use in iteration
        intersections = []
        intersections_y = []
        for i in range(len(self.sizes)):
            U_L_1 = results[self.sizes[i]]["U_L"]
            U_L_2 = results[self.sizes[(i + 1) % len(self.sizes)]]["U_L"]
            # print("U_L_1 = ", U_L_1)
            T_arr = results[self.sizes[i]]["T"]
            intersection, U_L_intersection = find_intersection(T_arr, U_L_1, U_L_2)
            intersections.append(intersection)
            intersections_y.append(U_L_intersection)

        T_c = np.mean(intersections)
        U_L_intersection = np.mean(intersections_y)
        T_c_error = np.ptp(intersections)

        fig, ax = plt.subplots(1, 1)

        y_upper_lim = 0
        y_lower_lim = np.infty
        shown_inds = np.linspace(0, len(self.sizes), len(self.sizes) + 1, endpoint=True,
                                 dtype=np.int64)
        ind = 0
        max_T = np.max(self.T_arr)
        min_T = np.min(self.T_arr)
        for i, size in enumerate(sorted(results.keys())):
            if i in shown_inds:
                T = np.array(results[size]["T"])
                U_L = np.array(results[size]["U_L"])
                ax.plot(T, U_L, linestyle="-", marker="x", color=colors[ind])
                ind += 1
                if max_T:
                    y_upper_lim = np.maximum(
                        np.max(U_L[(min_T < T) & (T < max_T)]), y_upper_lim)
                    y_lower_lim = np.minimum(
                        np.min(U_L[(min_T < T) & (T < max_T)]), y_lower_lim)

        y_span = y_upper_lim - y_lower_lim

        ax.set_xlabel("T")
        ax.set_ylabel(r"$U_L$")
        ax.set_title("Binder Cumulant on T")
        if min_T:
            ax.set_xlim(min_T, ax.get_xlim()[1])
            ax.set_ylim(y_lower_lim - 0.2 * y_span, y_upper_lim + 0.2 * y_span)
        if max_T:
            ax.set_xlim(ax.get_xlim()[0], max_T)
            ax.set_ylim(y_lower_lim - 0.2 * y_span, y_upper_lim + 0.2 * y_span)
        mark_point(ax, T_c, U_L_intersection,
                   label=rf"$T_c = {T_c:.4f}$")
        configure_ax(fig, ax)
        fig.savefig(self.simulation_path + "/cum_time_avg.png", format="png",
                    dpi=300, transparent=False)
        plt.show()

        # constructing cum dic
        cum_dic = {}
        for size in results:
            cum_dic[size] = results[size]["U_L"]

        diff_arr, size_arr = calc_diff_at(T_c,
                                          list(results.values())[0]["T"],
                                          cum_dic)

        popt, _ = curve_fit(linear_fit, np.log(size_arr),
                            np.log(diff_arr))
        nu = 1 / popt[0]

        print("FITTING RESULTS:")
        print("nu = ", nu)

        fig, ax = plt.subplots(1, 1)
        L_fit = np.linspace(0, np.max(size_arr) + 0.2 * np.max(size_arr), 101)
        ax.plot(L_fit, poly(L_fit, 1 / nu, np.exp(popt[1])),
                label=rf"$\nu = {nu:.2f}$", color=colors[0])
        ax.plot(size_arr, diff_arr, linestyle="", marker="x", color=colors[0])
        ax.set_xlabel("L")
        ax.set_ylabel(r"$\frac{d U_L}{d \varepsilon}$")
        ax.legend()
        ax.set_title(
            r"$\frac{d U_L}{d \varepsilon}$ for different System sizes $L$")
        configure_ax(fig, ax)
        # save_plot(root, "/critical_exponent.pdf", format="pdf")
        fig.savefig(self.simulation_path + "/critical_exponent_time_avg.png",
                    format="png", dpi=250, transparent=False)
        plt.show()

    def iteration(self, file, folder, user, wait, walltime):
        self.iteration_nr += 1
        # Here I want to have something that checks whether there is already a measurement
        simulation_available = check_directory_structure(self.sizes, self.T_arr, self.simulation_path)
        print("simulation_available", simulation_available)
        if not simulation_available or self.repeat:
            self.repeat = False                 # we set repeat to false again
            para_nr = self.write_para_files()  # setting up the parameter files for every simulation
            self.run_jobs(file, folder, para_nr, user, wait, walltime)
        else:
            print("Found valid simulation, evaluating")
        # after running the jobs, we need to
        # calculate the binder cumulant for every run
        return self.evaluate_simulation(file, folder, user, wait, walltime)

    def evaluate_simulation(self, file, folder, user, wait, walltime):
        self.all_T_arr = np.concatenate((self.all_T_arr, self.T_arr))   # we do it here, before evaluating
        threshold = 0.1  # in the simulation we calculated the mean from the last 90% of the values and achieved a small error
        all_results = self.construct_results(threshold, self.all_T_arr)     # for all results we use the self.all_T_arr? Since we want to reproduce always the measurements that we already did so that we can resume it where we left it of
        # current_results = self.construct_results(threshold, self.T_arr)
        # we switch to all results again but we set the min/max limites then also to min and max of all_T
        # We only enlarge the limits of all_T if we dont have an intersection anywhere
        # how does results look again precisely?
        # I have for every size a dictionary with {'T' : [T_min, ..., T_max], 'U_L': [U_L(T_min), ..., U_L(T_max)]}
        # once we have the U_L values the first thing we should do is check whether we have an intersection
        # something should be recursively if we need a second or third iteration...
        # The U_Ls should be the ones that we calculated this run..., so if we calculated a bad interval that does not contain Tc
        # that we can corr
        U_L_min_T_min = np.min(all_results[np.min(self.sizes)]['U_L'])
        U_L_max_T_min = np.min(all_results[np.max(self.sizes)]['U_L'])
        U_L_min_T_max = np.max(all_results[np.min(self.sizes)]['U_L'])
        U_L_max_T_max = np.max(all_results[np.max(self.sizes)]['U_L'])
        # we say we have an intersection if U_L_min_T_min > U_L_max_T_min
        # and U_L_min_T_max < U_L_max_T_max
        intersection = (U_L_min_T_min > U_L_max_T_min) & (
                U_L_min_T_max < U_L_max_T_max)
        if intersection:
            # good sign
            # determine Tc AND its error
            T_range, U_L_intersection, T_intersection, U_L_interpolated = interpolate_and_minimize(
                all_results)
            # I think T_range is the common Ts for the sizes, I think here every size should
            # definetely have the same temperatures
            # U_L_intersection is the U_L_value at the intersection
            # T_intersection is the temperature at the intersection
            # U_L_interpolated is a list with the numerical interpolated values at the stützstellen with resolution of
            # 1000, the resolution can be adjusted as a parameter
            # how do i find the size?
            # I think we can just assume that it is sorted
            # we now want to find the intersections of pairs of lines and want to
            # estimate an error out of this
            intersections = []
            #print("T_range")
            #print(T_range)
            #for i in range(len(U_L_interpolated)):
            #    nr_sizes = len(U_L_interpolated)
            #    U_L_1 = U_L_interpolated[i]
            #    U_L_2 = U_L_interpolated[(i + 1) % nr_sizes]
            #
            #    print("U_L_1 = ", U_L_1)
            #
            #    intersection = find_intersection(T_range, U_L_1, U_L_2)
            #    intersections.append(intersection)
            # TODO this works only for 3 different sizes
            for i in range(len(self.sizes)):
                U_L_1 = all_results[self.sizes[i]]["U_L"]
                U_L_2 = all_results[self.sizes[(i + 1) % len(self.sizes)]]["U_L"]
                #print("U_L_1 = ", U_L_1)
                T_arr = all_results[self.sizes[i]]["T"]
                intersection, _ = find_intersection(T_arr, U_L_1, U_L_2)
                intersections.append(intersection)
            # more is it net?
            # simple error would be to be just max_intersection - min_intersection?
            T_c = np.mean(intersections)
            print(f"Found an intersection at T_c = {T_c}")
            print("intersections: ", intersections)
            T_c_error = np.ptp(intersections)
            print(f"T_c_error = {T_c_error}")
            # relative error
            rel_intersec_error = T_c_error / T_c
            print(f"rel_intersec_error = {rel_intersec_error}")
            if rel_intersec_error < self.max_rel_intersec_error:
                # best case, now we are done?
                print(f"Determined crit. Temp T_c = {T_c} +- {rel_intersec_error}")
                return T_c, T_c_error
            else:
                # we check how many iterations we did so that we are not in an endless loop
                if self.iteration_nr > self.maximum_iterations:
                    print(f"Doing to many iterations, returning Temp T_c = {T_c} +- {T_c_error}")
                    return T_c, T_c_error
                # case that the error is to large, meaning we are moving closer to the critical temperature
                T_min = max(T_c - 5 * T_c_error, np.min(self.all_T_arr))  # standard smaller interval of 4*T_c_error?
                T_max = min(T_c + 5 * T_c_error, np.max(self.all_T_arr))
                if (T_min == np.min(self.T_arr)) & (T_max == np.max(self.T_arr)):
                    self.repeat = True
                # what if the new interval is larger than the last one because the errors were so large?
                # therefore we added max and min. If thats the case we want to reduce the equil error and repeat the
                # simulation
                # we need to set some variable so that the check folder structure function does not return that
                # the simulation is already done

                self.T_arr = np.linspace(T_min, T_max, self.nr_Ts)
                print(f"Error was too large: Temp T_c = {T_c} +- {T_c_error} \n"
                      f"Starting new run with T_min = {T_min}, T_max = {T_max}")
                # If we are moving closer to the critical point we should decrease the allowed error
                self.equil_error /= 2  # standard devide by two or is that not enough?
                return self.iteration(file, folder, user, wait, walltime)
        else:
            # bad sign, we do not seem to have an intersection
            # the question is now whether we are below or above
            print("We do not see a intersection")
            if U_L_min_T_max > U_L_max_T_max:
                print("The maximum temperature is too low")
                # this means we are below the critical temperature.
                # TODO can we somehow approximate how far away we are?
                # for now we just double the range i would say?
                T_range = np.ptp(self.all_T_arr)
                T_min = np.max(self.T_arr) + (self.T_arr[1] - self.T_arr[
                    0])  # the next T_min is the current maximum plus the current stepsize
                T_max = np.max(self.T_arr) + T_range
                self.T_arr = np.linspace(T_min, T_max, self.nr_Ts)
                # Okay we updated the temperature array... now we just run everything again?
            elif U_L_min_T_min < U_L_max_T_min:
                print("The minimum temperature is too high")
                T_range = np.ptp(self.all_T_arr)
                T_min = np.maximum(np.min(self.T_arr) - T_range, 0.0)  # We should not consider negative temperatures
                T_max = np.min(self.T_arr) - (self.T_arr[1] - self.T_arr[0])
                self.T_arr = np.linspace(T_min, T_max, self.nr_Ts)
            return self.iteration(file, folder, user, wait, walltime)

    def construct_results(self, threshold, selected_temps=None):
        results = {}
        for size_folder in os.listdir(self.simulation_path):
            if (size_folder[0] != ".") & (size_folder != "plots"):
                size_folder_path = os.path.join(self.simulation_path,
                                                size_folder)
                if os.path.isdir(size_folder_path):
                    size_result = process_size_folder(size_folder_path,
                                                      threshold, selected_temperatures=selected_temps)
                    results[int(size_folder)] = size_result
        return results

    def run_jobs(self, file, folder, para_nr, user, wait, walltime):
        # how do we make this routine? First we can make one cycle and submit nr gpus jobs?
        # or just the routine that will be waiting for the jobs to finish instantly?
        # establish the connection to hemera
        self.connection = Connection('hemera')
        print("connecting to hemera...")
        next_job = 0
        self.completed_jobs = set()     # we have to reset the completed jobs otherwise the program thinks we already compleated all the jobs
        self.nr_GPUS = min(self.nr_GPUS, self.total_runs)   # we dont need more GPUS than we have jobs
        while (len(self.completed_jobs)) < self.total_runs:
            # after this the set of running jobs is guaranteed to be empty
            # now we should check wheter some jobs are completed
            # just by getting which jobs are pending etc?
            queue_command = f"squeue -u {user}"
            queue_feedback = self.connection.run(queue_command)
            jobs_on_hemera = extract_numbers_after_newline(
                queue_feedback.stdout)

            # now we just check if all ids in running jobs are still in jobs_on_hemera
            # if not, we check if the id is completed
            # if thats true, we move the job id to the completed jobs
            print("jobs on hemera: ", jobs_on_hemera)
            print("running jobs: ", self.running_jobs)
            just_completed_jobs = set()
            for job_id in self.running_jobs:
                # check if it is still in jobs_on_hemera
                if job_id not in jobs_on_hemera:
                    # if not, we check what the sacc thingy command says
                    job_status_command = f'sacct -j {job_id} -o jobid,submit,start,end,state'
                    status_feedback = self.connection.run(job_status_command)
                    completed = check_completed_status(job_id,
                                                       status_feedback.stdout)
                    if completed:
                        # if its completed we double checked it
                        # the job id can be removed from running jobs and added to completed_jobs
                        just_completed_jobs.add(job_id)
                        self.completed_jobs.add(job_id)
                        # we rsync the new files
                        subprocess.call("./rsync.sh", cwd=pathlib.Path.home())
                    else:
                        # if it is not completed but not in jobs_on_hemera anymore,
                        # we have a problem
                        print(
                            "Oups! The Job vanished. Please check what happend")
            # remove the completed jobs from the running jobs list
            for job_id in just_completed_jobs:
                self.running_jobs.remove(job_id)

            # now we need to determine how many jobs are currently running
            # and we need to know how to submit jobs
            # while the number running jobs is smaller than the number of GPUs
            # to use, we submit new jobs
            # indeed it would be optimal if the long running jobs would start first
            # but that is hard to do? For now we just submit the jobs in arbitrary
            # order
            # we know which job is next through the next job variable,
            # but we still need to know the number at which to start
            while (len(self.running_jobs) < self.nr_GPUS) & (len(self.completed_jobs) + len(self.running_jobs) < self.total_runs):          # I think we forgot that we should not submit more jobs than we were expecting to? Is a job always either running or completed?
                # if this is true we are definetly going to submit a new job
                # so we can construct the para set string and advance next job
                para_set_nr = str(para_nr) + str(next_job)
                next_job += 1
                # now... we just submit it?
                # the command is:
                submit_job = f'sbatch --time {walltime} --output logs/%j.log --error' \
                             f' logs/errors/%j.err run_cuda.sh {file} {folder} {para_set_nr}'
                submit_feedback = self.connection.run(submit_job)
                job_id = extract_numbers_before_newline(submit_feedback.stdout)[
                    0]  # extract numbers before newline returns a list
                print(f"with para nr {para_set_nr}")
                # I think we can be sure that the job is running if we just commited it
                # better double check?
                self.running_jobs.add(job_id)
            # now we just wait some time before we do the next check?
            time.sleep(wait)
        # if we are here that means that all runs are done
        # we add the currently simulated temperatures to the bookkeeping variable

    def write_para_files(self, para_nr=100):
        # you ..., you know that you have to construct the parameter file at hemera?
        # and you need to do rsync after the jobs are finished!
        print("Writing the parameter files...")
        for i, (T, size) in enumerate(product(self.T_arr, self.sizes)):
            # We now need to construct the parameterfile with the appropriate temperature and size
            # Is it okay if we construct all files in the beginning and deal with the threading of the gpus later?
            # to construct the para set we need to know how many subsystems we should initialize
            nr_subsystems = int(self.nr_sites / (size ** 2 * self.Ly_Lx))
            with open(self.filepath + "/parameters/para_set_" + str(para_nr) + str(i) + '.txt', 'w') as f:
                f.write(self.simulation_path)
                f.write(f"\nend_time, {self.max_time} \n"
                        f"dt, {self.dt} \n"
                        f"J, {self.J_para} \n"
                        f"Jy, {self.J_perp} \n"
                        f"alpha, {self.h} \n"
                        f"eta, {self.eta} \n"
                        f"nr_saves, 4 \n"
                        f"nr_repeat, 0 \n"
                        f"min_temp, {T} \n"
                        f"max_temp, {T} \n"
                        f"nr_runs, 0.0 \n"
                        f"random_init, 0.0 \n"
                        f"curand_random, 1 \n"
                        f"subsystem_min_Lx, {size} \n"
                        f"subsystem_max_Lx, {size} \n"
                        f"nr_subsystem_sizes, 0  \n"
                        f"nr_subsystems, {nr_subsystems} \n"
                        f"x_y_factor, {self.Ly_Lx} \n"
                        f"nr_corr_values, 0 \n"
                        f"nr_ft_values, 0 \n"
                        f"equil_error, {self.equil_error}")
        # we need to copy the files to hemera
        rsync_command = ["rsync", "-auv", "--rsh", "ssh",
                         f"{self.filepath}/parameters/",
                         "hemera:~/Code/Master-Arbeit/CudaProject/parameters/"]
        subprocess.run(rsync_command, cwd=pathlib.Path.home())
        return para_nr

class quench_measurement():
    def __init__(self, J_para, J_perp, h, eta, dt, filepath, simulation_path, Tc, nr_GPUS=6, size_min=64,
                 size_max=4096, nr_sites=5e5, Ly_Lx=1/8, min_quench_steps=100, min_nr_sites=1e6,
                 min_nr_systems=10, host="hemera", user="weitze73"):
        self.J_para = J_para
        self.J_perp = J_perp
        self.h = h
        self.eta = eta
        self.dt = dt
        self.filepath = filepath
        self.simulation_path = simulation_path
        self.nr_GPUS = nr_GPUS
        self.size_min = size_min        # The starting size at which we do go higher
        self.size_max = size_max        # maximum size, if xi = size_max / 10 we stop the simulation
        self.nr_sites = nr_sites        # nr of sites for one run
        self.Ly_Lx = Ly_Lx              # the usual ratio between Lx and Ly
        self.Tc = Tc                    # The critical temperature that was calculated in the measurement before
        self.min_quench_steps = min_quench_steps    # the minimum number of steps done during a quench, influences tau_min
        self.min_nr_sites = min_nr_sites            # the minimum nr of sites i want to simulate for every tau
        self.min_nr_systems = min_nr_systems        # the minimum nr of systems to go through the qunech, guaranteeing an approximately precice xi value
        self.host = host                            # adress of the cluster
        self.user = user                            # user on the cluster


        self.T_start = Tc
        self.T_end = Tc
        self.tau_min = 1
        self.tau_factor = 1                    # current tau
        self.tau_list = []                 # empty list to keep track of the taus that we used
        self.size = size_min            # current size, starts at size_min
        self.equil_error = 0.02        # the equilibration error, so the error of U_L for which we assume that we are approximately equilibrated, doesnt need to be as small as for the T_c measurement
        self.para_nr = 100
        self.walltime = "24:00:00"
        self.file = "AutoQuench.cu"
        self.folder = "simulations"
        self.wait = 20
        self.cut_zero_impuls = True     # will probably always be true since we are quenching
        self.fitfunc = lorentz_offset      # We usually use this atm?
        self.nr_measured_values = 300   # standard number of measured values during the quench
        self.min_tau_scaling_fit = 10

        self.connection = None
        self.completed_jobs = set()     # bookkeeping of completed jobs
        self.running_jobs = set()       # bookkeeping of running jobs
    def setup(self):
        # We decide the start and end temperature
        self.T_start = 3 / 2 * self.Tc
        self.T_end = 1/2 * self.Tc
        # the minimum tau is supposed to be 1 / Tc, but depends on the stepsize
        min_pow = np.ceil(np.log2(self.min_quench_steps * self.dt / self.Tc))
        self.tau_factor = min_pow
        self.tau_min = 2 ** min_pow
        print("Tau min: ", self.tau_min)
        print("Tau factor: ", self.tau_factor)
        # We also add the tau to the list
        self.tau_list.append(self.tau_min)

    def tau(self):
        return 2 ** self.tau_factor

    def systems_per_job(self):
        return int(self.nr_sites / (self.size ** 2 * self.Ly_Lx))

    def write_para_file(self):
        # I think we only construct one file at a time and ggf run it multiple times.
        # because we have to wait for the result auf the last tau to decide whether to increase the system size for the next tau
        print("Writing the parameter file...")
        nr_subsystems = self.systems_per_job()
        with open(self.filepath + "/parameters/para_quench_set_" + str(self.para_nr) + '.txt', 'w') as f:
            f.write(self.simulation_path)
            f.write(f"\n"
                    f"dt, {self.dt} \n"
                    f"J, {self.J_para} \n"
                    f"Jy, {self.J_perp} \n"
                    f"alpha, {self.h} \n"
                    f"eta, {self.eta} \n"
                    f"nr_saves, 2 \n"
                    f"nr_repeat, 0 \n"
                    f"starting_temp, {self.T_start} \n"
                    f"end_temp, {self.T_end} \n"
                    f"nr_runs, 0.0 \n"
                    f"min_tau_factor, {self.tau_factor} \n"
                    f"max_tau_factor, {self.tau_factor} \n"
                    f"random_init, 1.0 \n"
                    f"curand_random, 1 \n"
                    f"subsystem_min_Lx, {self.size} \n"
                    f"subsystem_max_Lx, {self.size} \n"
                    f"nr_subsystem_sizes, 0  \n"
                    f"nr_subsystems, {nr_subsystems} \n"
                    f"x_y_factor, {self.Ly_Lx} \n"
                    f"nr_corr_values, {self.nr_measured_values} \n"
                    f"nr_ft_values, {self.nr_measured_values} \n"           # TODO probably deprecated, have to see later!
                    f"equil_error, {self.equil_error} \n"
                    f"min_cum_nr, 50")
        # we need to copy the files to hemera
        rsync_command = ["rsync", "-auv", "--rsh", "ssh",
                         f"{self.filepath}/parameters/",
                         "hemera:~/Code/Master-Arbeit/CudaProject/parameters/"]
        subprocess.run(rsync_command, cwd=pathlib.Path.home())

    def iteration(self):
        # Check directory structure again
        sim_available = check_directory_structure([self.size], [self.tau()], self.simulation_path) # the check directory structure takes lists in
        # one of the first things we do in the iteration is writing the parameter file
        if not sim_available:
            # If no simulation is available, we have to run the simulation
            self.write_para_file()
            # then we just run the jobs
            self.run_jobs()
        # now the evaluation of the job and what we do accordningly
        return self.evaluate()      # Is this proper design actually? letting iteration return evalutate while evaluate calls iteration recursively?

    def evaluate(self):
        # do we have to write this new or can we use stuff from previous scripts?
        # Will we have a new kind of quench without equilibration in the end? Probably yes
        # meaning we should probably write new stuff to satisfy our needs...
        # we now the current measurement path
        cur_path = os.path.join(self.simulation_path, str(self.size), f"{self.tau():.6f}")
        # we need the fourier transforms of this measurement and we need to fit xi
        ft_k, ft_l = average_lastline_ft(cur_path)
        p_k = get_frequencies_fftw_order(len(ft_k))
        p_l = get_frequencies_fftw_order(len(ft_l))

        if self.cut_zero_impuls:
            p_k, ft_k = cut_zero_imp(p_k, ft_k)
            p_l, ft_l = cut_zero_imp(p_l, ft_l)

        popt_x, perr_x = fit_lorentz(p_k, ft_k, fitfunc=self.fitfunc)
        popt_y, perr_y = fit_lorentz(p_l, ft_l, fitfunc=self.fitfunc)
        xix = np.abs(popt_x[0])     # I dont think that we need the minimum here anymore since if xi is larger than L/10 we repeat the measurement?
        xiy = np.abs(popt_y[0])

        # Okay so we have now the xi value of the current simulation.
        # We now act accordingly.
        if (xix > self.size / 10) or (xiy > self.size * self.Ly_Lx / 10):
            # if in one of the two directions there is a correlation length that is too large, we repeat with twice
            # the system size
            print(f"The Correlation length is too large, xi_x = {xix},"
                  f"  xi_y = {xiy} with Lx = {self.size}, Ly = {self.size * self.Ly_Lx}")
            self.size *= 2

            if self.size > self.size_max:
                print("We exceed the maximum size, we stop the measurement here")
                return
            else:
                print(f"Repeating the simulation with the size of Lx = {self.size}")
                return self.iteration()
        else:
            # If it is not to large, we just start the next measurement with the larger tau?
            # Should we preventively enlarge the size if for example xix > self.size / 15 or something like that?
            # not for now i would say, keep it simple at this point
            self.tau_factor += 1       # increase tau
            self.tau_list.append(self.tau())    # and add it to the bookkeeping list
            return self.iteration()

    def run_jobs(self):
        # this method is responsible for running the jobs for one tau, one system size before moving to the next tau
        # first we need to find out how many jobs this will be
        # it doenst matter if this part of the code is optimized so just naive:
        # find out whether min_nr_systems or min_nr_sites is limiting:
        sys_per_job = int(np.ceil(self.nr_sites / (self.size ** 2 * self.Ly_Lx)))
        min_jobs = int(self.min_nr_sites / self.nr_sites)            # the minimum number of jobs is the minimum_nr of total sites divided by the number of sites per job
        if sys_per_job >= self.min_nr_systems / min_jobs:
            # if the number of systems per job * the minimumb number of jobs given by the minimum number of system sites
            # divided by the nr of sites per job is larger than the required amount of systems, wo only do the
            # minimum number of jobs
            nr_jobs = min_jobs
        else:
            # else we do as many jobs as we need to exceed the min_nr_systems
            nr_jobs = int(np.ceil(self.min_nr_systems / sys_per_job))
        # now we do basically the same stuff as last time, only that we onyl submit the same job
        self.connection = Connection('hemera')
        self.completed_jobs = set()
        while (len(self.completed_jobs)) < nr_jobs:
            # after this the set of running jobs is guaranteed to be empty
            # now we should check wheter some jobs are completed
            # just by getting which jobs are pending etc?
            queue_command = f"squeue -u {self.user}"
            queue_feedback = self.connection.run(queue_command)
            jobs_on_hemera = extract_numbers_after_newline(
                queue_feedback.stdout)

            # now we just check if all ids in running jobs are still in jobs_on_hemera
            # if not, we check if the id is completed
            # if thats true, we move the job id to the completed jobs
            print("jobs on hemera: ", jobs_on_hemera)
            print("running jobs: ", self.running_jobs)
            just_completed_jobs = set()
            for job_id in self.running_jobs:
                # check if it is still in jobs_on_hemera
                if job_id not in jobs_on_hemera:
                    # if not, we check what the sacc thingy command says
                    job_status_command = f'sacct -j {job_id} -o jobid,submit,start,end,state'
                    status_feedback = self.connection.run(job_status_command)
                    completed = check_completed_status(job_id,
                                                       status_feedback.stdout)
                    if completed:
                        # if its completed we double checked it
                        # the job id can be removed from running jobs and added to completed_jobs
                        just_completed_jobs.add(job_id)
                        self.completed_jobs.add(job_id)
                        # we rsync the new files
                        subprocess.call("./rsync.sh", cwd=pathlib.Path.home())
                    else:
                        # if it is not completed but not in jobs_on_hemera anymore,
                        # we have a problem
                        print(
                            "Oups! The Job vanished. Please check what happend")
            # remove the completed jobs from the running jobs list
            for job_id in just_completed_jobs:
                self.running_jobs.remove(job_id)

            # now we need to determine how many jobs are currently running
            # and we need to know how to submit jobs
            # while the number running jobs is smaller than the number of GPUs
            # to use, we submit new jobs
            # indeed it would be optimal if the long running jobs would start first
            # but that is hard to do? For now we just submit the jobs in arbitrary
            # order
            # we know which job is next through the next job variable,
            # but we still need to know the number at which to start
            while (len(self.running_jobs) < self.nr_GPUS) and (len(self.running_jobs) + len(self.completed_jobs) < nr_jobs):
                submit_job = f'sbatch --time {self.walltime} --output logs/%j.log --error' \
                             f' logs/errors/%j.err run_cuda.sh {self.file} {self.folder} {self.para_nr}'
                submit_feedback = self.connection.run(submit_job)
                job_id = extract_numbers_before_newline(submit_feedback.stdout)[
                    0]  # extract numbers before newline returns a list
                print(f"with para nr {self.para_nr}")
                # I think we can be sure that the job is running if we just commited it
                # better double check?
                self.running_jobs.add(job_id)
            # now we just wait some time before we do the next check?
            time.sleep(self.wait)

    def conclude(self):
        # like last time, I think for now I will just paste the avargeFTsOverTime script
        size_x_dic = {}
        size_y_dic = {}

        # for size in os.listdir(self.simulation_path):
        #     t_xix = {}
        #     t_xiy = {}
        #     if (size != "plots") & (size[0] != "."):
        #         sizepath = os.path.join(self.simulation_path, size)
        #         if os.path.isdir(sizepath):
        #             for setting in os.listdir(sizepath):
        #                 if (setting != "plots") & (setting[0] != "."):
        #                     settingpath = os.path.join(sizepath, setting)
        #                     print(settingpath)
        #                     setting = float(setting)
        #                     parapath = find_first_txt_file(self.simulation_path)
        #                     parameters = read_parameters_txt(parapath)
        #
        #                     Lx = parameters["subsystem_Lx"]
        #                     Ly = parameters["subsystem_Ly"]
        #
        #                     if os.path.isdir(settingpath):
        #                         t_xix[setting] = {}
        #                         t_xiy[setting] = {}
        #                         ft_k, ft_l = average_ft_unequal_times(settingpath)
        #                         for t in ft_k:
        #                             p_k = get_frequencies_fftw_order(len(ft_k[t]))
        #                             if self.cut_zero_impuls:
        #                                 p_k, ft_k[t] = cut_zero_imp(p_k, ft_k[t])
        #                             popt_x, perr_x = fit_lorentz(p_k, ft_k[t],
        #                                                          fitfunc=self.fitfunc)
        #                             xix = np.minimum(np.abs(popt_x[0]), Lx)
        #                             t_xix[setting][t] = xix
        #                         for t in ft_l:
        #                             p_l = get_frequencies_fftw_order(len(ft_l[t]))
        #                             if self.cut_zero_impuls:
        #                                 p_l, ft_l[t] = cut_zero_imp(p_l, ft_l[t])
        #                             popt_y, perr_y = fit_lorentz(p_l, ft_l[t],
        #                                                          fitfunc=self.fitfunc)
        #                             xiy = np.minimum(np.abs(popt_y[0]), Ly)
        #                             t_xiy[setting][t] = xiy
        #                 size_x_dic[int(size)] = t_xix.copy()
        #                 size_y_dic[int(size)] = t_xiy.copy()
        # okay now sadly that wont just work like that. I have to write the own funcitons.
        # What functions do I need? I have for every size just the dicionary with
        # size_x_dic = {
        #                   64: {   tau_1: {t11: xix11, t12: xix12, ...},
        #                           tau_2 :  {t11: xix11, t12: xix12, ...}, ...
        #                        },
        #                   128 {
        #                           tau_2 :  {t21: xix21, t12: xix22, ...},
        #                           tau_3 :  {t22: xix22, t22: xix22, ...}, ..
        #                       },...
        #               }
        # If we plot the process, it doesnt really matter which size we use since we made sure that the correlation length
        # is small enough. Though it would be better to use the larger systemsize, if it has one
        # We should specify a number of lines that should be shown
        # self.nr_process_curves = 5
        # # Now we have to select the taus that we want, We use the bookkeeping list therefore
        # # the tau list will be sorted, right?
        # plot_inds = np.linspace(0, len(self.tau_list), self.nr_process_curves, dtype=np.int32, endpoint=False)
        # tau_plot = np.array(sorted(self.tau_list))[plot_inds][::-1]         # reverse order, start with largest tau
        # # now we iterate through the sorted taus and look for the largest system size that can deliver it
        #
        # figx, axx = plt.subplots(1, 1)
        # figy, axy = plt.subplots(1, 1)
        # for tau in tau_plot:
        #     # focus! write this stuff down and then you are allowed to clean your keyboard
        #     # search the largest system size containing tau
        #     # If we start from behind with the largest tau, we can be sure that only the largest system size has this tau
        #     # we then continue to use this size until the tau we are looking for is not there anymore for this size
        #     # then we reduce the system size by one
        #     t_xix, t_xiy, largest_size = self.get_xi_of_tau(size_x_dic, size_y_dic, tau)      # returns the fitting dics aswell as the size that we used
        #     tx = list(t_xix.keys())
        #     xix = list(t_xix.values())
        #     ty = list(t_xiy.keys())
        #     xiy = list(t_xiy.values())
        #
        #     # Plotting
        #     # We sadly dont know the equilibration time atm
        #     # Oh damn I just  thought of the point that we atm dont know when we average runs If they were even
        #     # In the same state of simulation at the same time. The one system might haven taken some longer time
        #     # to equilibrate
        #     # We could plot just one run, this one would then look a bit shakey probably
        #     # we could just omit the plotting of the process and just plot the endresult,
        #     # although some kind of process would be nice?
        #     # We could determine the equilibration time for every run on its own and have its own dictionary and afterwards
        #     # take the mean by folding?
        #     # We have to deal with the equilibration time anyways, otherwise it would look weird?
        #     # So we need for every run its own dicitonary and we need to know what the equilibration time was and somehow
        #     # we have to combine all that into a huge dictionary...
        #     # Another problem is that for large systems we want to average the ft before fitting as the ft of a single system might
        #     # have large uncertainties. If we have different timepoints we can not really do that
        #     # we would have to select and map the timepoints... the first point after equilibration in all runs can be averaged.
        #     # The points before can not really be averaged. They cound be folded afterwards
        #     # damn this is all more complicated than expected, i will continue with first just looking at the quench exponent
        #     # so the dictionary now has the ebenen:
        #     # size_dic = { size: {tau: { equil_time: {t:
        #     axx.plot(tx/tau, xix, label=rf"$\tau = $ {tau}  $L_x = ${largest_size}", linestyle="", marker=".")
        #     axy.plot(ty/tau, xiy, label=rf"$\tau = $ {tau}  $L_y = ${largest_size * self.Ly_Lx}", linestyle="", marker=".")
        #
        # configure_ax(figx, axx)
        # configure_ax(figy, axy)
        # plt.show()

        # for the plotting of the exponent we only need the value at the end of the quench
        # just plot it for every size and every tau? worry about the fitting later
        size_tau_xix_dic, size_tau_xiy_dic = self.get_size_quench_results()

        # Now we want to do some fitting
        # I have the feeling that it will be faster if i just rewrite it
        # Okay so we will use from the size_tau_xi_dics for every size
        # all but the last value. Then we should have a continous tau-xi curve
        # About the minimum tau, I don't know how to deal with that atm,
        # we will just set it to 10 or let it be a parameter

        tau_scaling = []
        xix_scaling = []
        xiy_scaling = []

        for size in sorted(list(size_tau_xix_dic.keys())):
            tau_xix_dic = size_tau_xix_dic[size]
            tau_xiy_dic = size_tau_xiy_dic[size]
            # from those dictionaries I want to extract the largest tau. Are the taus floats here?
            for tau in np.sort(list(tau_xix_dic.keys()))[:-1]:
                # By the -1 we exclude hopefully the largest tau
                tau_scaling.append(tau)
                xix_scaling.append(tau_xix_dic[tau])
                xiy_scaling.append(tau_xiy_dic[tau])
        # We dont need a max tau in this case as the largest tau value is excluded and all other are basically guaranteed to be accurate
        xix_scaling = np.array(xix_scaling)[
            np.array(tau_scaling) > self.min_tau_scaling_fit]
        xiy_scaling = np.array(xiy_scaling)[
            np.array(tau_scaling) > self.min_tau_scaling_fit]
        tau_scaling = np.array(tau_scaling)[
            np.array(tau_scaling) > self.min_tau_scaling_fit]

        # Do the fitting
        popt_x, _ = curve_fit(linear_fit, np.log(tau_scaling),
                              np.log(xix_scaling))
        popt_y, _ = curve_fit(linear_fit, np.log(tau_scaling),
                              np.log(xiy_scaling))


        # now we have them dictionaries, we can plot the points
        figx, axx = plt.subplots(1, 1)
        # set the axes to be logarithmic
        axx.set_yscale("log")
        axx.set_xscale("log")
        axx.set_xlabel(r"$\tau$")
        axx.set_ylabel(r"$\xi_x$")
        for i, size in enumerate(size_tau_xix_dic):
            # Plotting, every size should get its own color and/or symbol?
            # construct the lists to plot
            tau = list(size_tau_xix_dic[size].keys())
            xix = list(size_tau_xix_dic[size].values())

            axx.plot(tau, xix, marker=markers[i], linestyle="None", label=rf"$L_x = ${size}", color="C0")
        # Plot the fit
        axx.plot(tau_scaling, poly(tau_scaling, popt_x[0], np.exp(popt_x[1])),
                 color="black", alpha=0.5, linestyle="dashed",
                 label=r"$\frac{\nu}{1 + \nu z} =$" + f"{popt_x[0]:.2f}")
        configure_ax(figx, axx)
        create_directory_if_not_exists(self.simulation_path+ "/plots")
        plt.savefig(self.simulation_path + "/plots/tau-xix.png", format="png")
        plt.show()

        figy, axy = plt.subplots(1, 1)
        axy.set_yscale("log")
        axy.set_xscale("log")
        axy.set_xlabel(r"$\tau$")
        axy.set_ylabel(r"$\xi_y$")
        for i, size in enumerate(size_tau_xiy_dic):
            # Plotting, every size should get its own color and/or symbol?
            # construct the lists to plot
            tau = list(size_tau_xiy_dic[size].keys())
            xiy = list(size_tau_xiy_dic[size].values())

            axy.plot(tau, xiy, marker=markers[i], linestyle="None", label=rf"$L_y = ${size * self.Ly_Lx}", color="C1")
        axy.plot(tau_scaling, poly(tau_scaling, popt_y[0], np.exp(popt_y[1])),
                 color="black", alpha=0.5, linestyle="dashed",
                 label=r"$\frac{\nu}{1 + \nu z} =$" + f"{popt_y[0]:.2f}")
        configure_ax(figy, axy)
        plt.savefig(self.simulation_path + "/plots/tau-xiy.png", format="png")
        plt.show()









    def get_size_quench_results(self):
        size_tau_xix_dic = {}
        size_tau_xiy_dic = {}
        for size in os.listdir(self.simulation_path):
            tau_xix = {}
            tau_xiy = {}
            if (size != "plots") & (size[0] != "."):
                sizepath = os.path.join(self.simulation_path, size)
                if os.path.isdir(sizepath):
                    for tau in os.listdir(sizepath):
                        if (tau != "plots") & (tau[0] != "."):
                            taupath = os.path.join(sizepath, tau)
                            if os.path.isdir(taupath):

                                ft_k, ft_l = average_lastline_ft(taupath)
                                p_k = get_frequencies_fftw_order(len(ft_k))
                                p_l = get_frequencies_fftw_order(len(ft_l))

                                if self.cut_zero_impuls:
                                    p_k, ft_k = cut_zero_imp(p_k, ft_k)
                                    p_l, ft_l = cut_zero_imp(p_l, ft_l)

                                popt_x, perr_x = fit_lorentz(p_k, ft_k, fitfunc=self.fitfunc)
                                popt_y, perr_y = fit_lorentz(p_l, ft_l, fitfunc=self.fitfunc)
                                xix = np.abs(popt_x[
                                                 0])  # Those are the xi values for a certain size and tau at the end of the quench
                                xiy = np.abs(popt_y[0])

                                # We add them to the tau_xi dictionaries
                                tau_xix[float(tau)] = xix
                                tau_xiy[float(tau)] = xiy
                # Now we add them to the outer dictionary
                size_tau_xix_dic[int(size)] = tau_xix
                size_tau_xiy_dic[int(size)] = tau_xiy
        return size_tau_xix_dic, size_tau_xiy_dic

    def get_xi_of_tau(self, size_x_dic, size_y_dic, tau, cur_size=None):
        if cur_size == None:
            # this means that we call the function the first time and are supposed to use the largest size available
            cur_size = np.max(list(size_x_dic.keys()))
        tau_t_xix = size_x_dic[cur_size]
        tau_t_xiy = size_y_dic[cur_size]
        # this should now contain the tau
        try:
            t_xix = tau_t_xix[tau]
            t_xiy = tau_t_xiy[tau]
            return t_xix, t_xiy, cur_size
        except KeyError as e:
            # Keyerror means that the tau we want is not available for the largest size
            # we Have to reduce the current size
            avail_sizes = sorted(list(size_x_dic.keys()))
            # get the index of the current size
            cur_index = avail_sizes.index(cur_size)
            # we reduce it by one
            cur_index -= 1
            if cur_index < 0:
                # should not happen, but who knows?
                print("Index is below zero, means we didnt find some tau?")
                return
            cur_size = avail_sizes[cur_index]
            return self.get_xi_of_tau(size_x_dic, size_y_dic, tau, cur_size)

    def run(self):
        self.setup()
        self.iteration()
        self.conclude()

class amplitude_measurement():
    def __init__(self, J_para, J_perp, h, eta, dt, filepath, simulation_path, Tc, nr_GPUS=6, nr_Ts=5, size=1024,
                 max_steps=1e9, Ly_Lx = 1/8, equil_error=0.003, T_range_fraction=0.05):
        self.J_para = J_para
        self.J_perp = J_perp
        self.h = h
        self.eta = eta
        self.dt = dt
        self.filepath = filepath                # Path where the parameter files are stored
        self.simulation_path = simulation_path  # path where the simulation data is stored
        self.nr_GPUS = nr_GPUS
        self.nr_Ts = nr_Ts                      # nr of temperatures used to fit
        self.T_range_fraction = T_range_fraction    # This is the fraction of Tc that is used to determine the interval of Ts in [Tc, (1+T_range_raction) * Tc]
        self.Tc = Tc                            # We need to know the critical temperature that we determined in previous measurements
        # I think we only use one size, being the largest possible
        # Is there some way of guessing how long a simulation will take?
        # 1000 is a save bet but I guess we would like to go up to even larger sizes?
        self.size = size
        self.max_steps = max_steps
        self.Ly_Lx = Ly_Lx

        self.T_arr = np.array([])
        self.total_runs = 0                 # just a parameter to keep track of how many jobs we want to run in this iteration
        self.running_jobs = set()
        self.completed_jobs = set()
        self.connection = None

        self.all_T_arr = np.array([])       # Bookkeeping for all the temperatures we have simulated in this setting
        self.equil_error = equil_error           # standard equilibration error for the xi runs
        self.maximum_iterations = 2         # we first look at the interval [Tc, 1.05Tc] and If this doesnt work we inrease to [Tc, 1.1Tc] and If this doesnt work we abort
        self.iteration_nr = 0
        self.min_corr_nr = 500
        self.equil_cutoff = 0.3             # This is the values that we cut off because we think we are still equilibrating. Since we definitely want the values in equilibration we use a relatively large cutoff here


    def setup(self):
        # This function will determine the initial T range
        # we want nr_Ts temperatures between
        Tmax = self.Tc * (1 + self.T_range_fraction)
        self.T_arr = np.linspace(self.Tc, Tmax, num=self.nr_Ts)
        # Is that all, what else do we need to setup?
        # For now the total number of runs is just the nr of temps
        self.total_runs = self.nr_Ts
    def run(self):
        # runs the complete simulation. Some initialization, a routine and a finish, plotting and fitting
        self.setup()

    def iteration(self):
        # Okay again this is the function that runs the jobs, evaluates them and runs additional jobs if necessary. Also this
        # this function should be able to pick up on simulations that were not finished
        # increase the iteration number
        self.iteration_nr += 1
        # Check if we have a simulation available
        # I think we should improve the pickup capability, if we run for example half the jobs of simulation we should be able to use them

def main():
    # okay what is the first thing we need to do?
    # we need parameters like the number of gpus we are able to use
    nr_gpus = 15
    # we somehow need the relevant parameters
    # The model defining parameters are J_perp J_para h eta
    # the simulation defining parameters are dt
    J_para = -3.11
    J_perp = -0.1
    h = 0.5
    eta = 1.5
    dt = 0.01
    max_size_Tc = 80
    min_size_Tc = 48
    nr_sizes_Tc = 3
    #filepath = "/home/andi/Studium/Code/Master-Arbeit/CudaProject"
    filepath = "/home/weitze73/Documents/Master-Arbeit/Code/Master-Arbeit/CudaProject"
    simulation_path = "../../Generated content/Silicon/Subsystems/Suite/Test4/"

    max_rel_intersection_error = 0.001

    # Quench parameters
    max_size = 1024
    min_nr_sites = 1e6

    # Enter which calculations are supposed to run here
    measurements = {
        "Tc": False ,
        "Quench": True
    }

    # I honestly have no idea on how to account h, that is really a problem
    # the Scanned interval
    if measurements["Tc"]:
        sim = crit_temp_measurement(J_para, J_perp, h, eta, dt, filepath, simulation_path + "Tc", nr_GPUS=nr_gpus,
                                    size_min=min_size_Tc, size_max=max_size_Tc, nr_sizes=nr_sizes_Tc,
                                    intersection_error=max_rel_intersection_error)
        T_c, T_c_error = sim.routine()
    else:
        T_c = float(input("Enter critical temperature:"))
        T_c_error = 0
    if measurements["Quench"]:
        quench = quench_measurement(J_para, J_perp, h, eta, dt, filepath, simulation_path + "Quench", T_c, nr_GPUS=nr_gpus, size_max=max_size, min_nr_sites=min_nr_sites )
        quench.run()


if __name__ == '__main__':
    main()