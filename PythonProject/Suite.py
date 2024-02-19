import os

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
from glob import glob

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
        size_path = os.path.join(directory_path, str(int(size)))

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

def get_avail_simulations(sizes, temperatures, directory_path, check_function, check_function_args):
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
                temp_path = os.path.join(size_path, f"{temp:.6f}")

                # Check if the temperature folder exists and is valid based on the custom function
                if (os.path.exists(temp_path) and os.path.isdir(temp_path)):
                        sim_valid = check_function(temp_path, *check_function_args)
                        if sim_valid:
                            valid_folders.append((size, temp))

    return valid_folders

def check_corr_valid(folderpath, equil_error, equil_cutoff):
    # This function is supposed to check whether the .corr file in the folder (There should only be one)
    # has a low enough error, the error is specified by equil error in the class
    # We also need the threshold to calculate the accurate values with the accurate errors
    xix_avg, xix_error = process_temp_folder(folderpath, equil_cutoff, value="xix", file_ending="corr")
    xiy_avg, xiy_error = process_temp_folder(folderpath, equil_cutoff, value="xiy", file_ending="corr")

    # We have the error, so we can check if it is small enough
    # How do we check in the observer again? is not written actually...
    # We need the relative errors here
    xix_rel_error = xix_error / xix_avg
    xiy_rel_error = xiy_error / xiy_avg
    avg_error = 1/ 2 * (xix_rel_error + xiy_rel_error)

    if avg_error < equil_error:
        return True
    else:
        return False

def check_cum_valid(folderpath, variation_error_rate, num_per_var_val, ds):
    """
    this function checks if a simulation in folderpath has .cum file with
    a variation error that is small than given
    :param folderpath: path to the folder containing the simulation
    :param num_per_var_val: the number of values we use to calculate one dif_std value
    :param ds: time between the measurements
    :return: bool whether the simulation is valid or not
    """
    # I dont think we need anything else do do this?
    # first we need to read in every .cum file and average it. We assume that they have
    # the same timepoints aswell as the same number of values?
    cum, times = get_folder_average(folderpath)
    # I am still not sure if we should do the thing with the nr of points per stddev value
    # you know what I think we should do it, if would be like linearly approximating every part of the curve
    # if the curve would vary very smoothly, the variation error that we want should approach zero
    # if we do not subdivide the cumulant curve, the stddev will somehow be capped at a minimum value
    # the number of intervals is integer division of the nr of cum values devided
    # by the number of values in one interval. If this goes not perfectly
    # up, we add one to have one more interval
    dif_var_mean, squared_dif_var_mean = get_mean_dif_var(cum, num_per_var_val)
    # so the thing is that this is dependent on how large the stepsize is, do we have to correlate them again
    # or we consider something more like the dif_var_mean relative to the duration of the intervals
    # if the duration of the intevals is really small, the dif_var_mean will be small
    dif_var_mean_rate = dif_var_mean / (num_per_var_val * ds)
    squared_dif_var_mean_rate = squared_dif_var_mean / (num_per_var_val * ds)
    # If this mean is small than the validation error that we want to use, we accept the measurement
    if (dif_var_mean_rate < variation_error_rate) and (squared_dif_var_mean_rate < variation_error_rate ** 2):
        return True
    else:
        return False


def get_mean_dif_var(cum,num_per_var_val):
    nr_intervals = len(cum) // num_per_var_val + (len(cum) % num_per_var_val != 0)
    # we change the behavior, strong flucutations at the beginning are averaged out for the large sizes in the long run
    # num_per_var_val = len(cum) // nr_intervals
    dif_var_arr = []
    for i in range(nr_intervals):
        # We have to be careful if we are at the last i?
        cum_interval = cum[
                       i * num_per_var_val:min((i + 1) * num_per_var_val, len(cum))]  # the cum values in this interval
        # from those values we want to get the stddev of the differences
        dif_var = get_difference_std(cum_interval)  # and now we collect them in an array?
        dif_var_arr.append(dif_var)
    # If we are finished we just take the mean of those var values?
    # we choose the 5 largest variances
    dif_var_arr = np.array(sorted(dif_var_arr))[-5:]
    dif_var_mean = np.mean(dif_var_arr)
    squared_dif_var_mean = np.mean(np.array(dif_var_arr) ** 2)    # this way large variances are punished harder
    return dif_var_mean, squared_dif_var_mean


def get_folder_average(folderpath, file_ending="cum", value="U_L"):
    """
    This is supposed to look at a file that has observables like corr length and average them if we have multiple files
    this one only works for same times
    :param folderpath: path to the simulation
    :param file_ending: file ending of the file that holds the observables
    :param value: name of the value in the file
    :return: averaged value array, should it maybe also return the time?
    """
    cum_files = glob(f"{folderpath}/*.{file_ending}")
    cum = np.array([])
    times = np.array([])
    for i, cum_path in enumerate(cum_files):
        df = pd.read_csv(cum_path, delimiter=",", index_col=False)
        this_cum = np.array(df[value])
        # I think I dont need the times, I just want to calculate the difference variation
        if i == 0:
            cum = this_cum
            times = df['t']
        else:
            cum += this_cum
    # averaging
    cum = cum / len(cum_files)
    times = np.array(times)
    return cum, times


def get_difference_std(cum_arr):
    """
    Calculates the stddev of the difference of consecutive cumulant errors
    :param cum_arr: array with the cumulant values
    :return: stddevs of the differences
    """
    # first we have to calculate the array with the differences, maybe numpy
    # actually has something for this
    differences = np.ediff1d(cum_arr)
    # and now we just want the stddev of this, or maybe the relative stddev?
    dif_var = np.std(differences)
    # I dont think we want to use the relative difference variance, because the difference variance is already
    # a measure of how much the difference flucutates, if the difference is large, but doesnt fluctuate, the variance
    # will already be small
    # rel_dif_var = dif_var / (np.mean(np.abs(differences)))
    return dif_var

class autonomous_measurement():
    def __init__(self, J_para, J_perp, h, eta, dt, filepath, simulation_path, exec_file, nr_GPUS=6, Ly_Lx = 1/8,
                 host="hemera", user="weitze73"):
        # This class is supposed to encapsulate some of the functionality that the following classes share
        # For every simulation I need the basic simulation parameters
        self.J_para = J_para
        self.J_perp = J_perp
        self.h = h
        self.eta = eta
        self.dt = dt
        self.filepath = filepath
        self.simulation_path = simulation_path
        self.nr_GPUS = nr_GPUS
        self.Ly_Lx = Ly_Lx

        # also some parameters for the cluster
        self.host = host                            # adress of the cluster
        self.user = user                            # user on the cluster
        self.walltime = "12:00:00"
        self.file = exec_file
        self.folder = "simulations"
        self.wait = 60

        # Besides the external simulation parameters that I have to provide there are some other attributes that
        # every autonomous suite needs
        # The bookkeeping variables for commiting the jobs
        self.running_jobs = set()
        self.completed_jobs = set()
        self.connection = None
        self.total_runs = 0
        self.para_nr = 100                            # every mearsurement has at least one current parameter number


    def check_completed_jobs(self):
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

    def submit_jobs(self, file=None):
        if file == None:
            file = self.file
        while (len(self.running_jobs) < self.nr_GPUS) & (len(self.completed_jobs) + len(
                self.running_jobs) < self.total_runs):  # I think we forgot that we should not submit more jobs than we were expecting to? Is a job always either running or completed?
            # if this is true we are definetly going to submit a new job
            # so we can construct the para set string and advance next job
            # TODO If i overwrite this get para function in a child class and call THIS method, submit jobs, the function should
            # be overwritten here right? I think you make yourselfs to many thoughts, In python this should be way easier than in c++
            para_nr = self.get_para_nr()
            # now... we just submit it?
            # the command is:
            submit_job = f'sbatch --time {self.walltime} --output logs/%j.log --error' \
                         f' logs/errors/%j.err run_cuda.sh {file} {self.folder} {para_nr}'
            submit_feedback = self.connection.run(submit_job)
            job_id = extract_numbers_before_newline(submit_feedback.stdout)[
                0]  # extract numbers before newline returns a list
            print(f"with para nr {para_nr}")
            # I think we can be sure that the job is running if we just commited it
            # better double check?
            self.running_jobs.add(job_id)

    def get_para_nr(self):
        return self.para_nr
    def run_jobs(self, file=None):
        # how do we make this routine? First we can make one cycle and submit nr gpus jobs?
        # or just the routine that will be waiting for the jobs to finish instantly?
        # establish the connection to hemera
        self.connection = Connection('hemera')
        print("connecting to hemera...")
        self.completed_jobs = set()     # we have to reset the completed jobs otherwise the program thinks we already compleated all the jobs
        self.nr_GPUS = min(self.nr_GPUS, self.total_runs)   # we dont need more GPUS than we have jobs
        while (len(self.completed_jobs)) < self.total_runs:
            # after this the set of running jobs is guaranteed to be empty
            # now we should check wheter some jobs are completed
            # just by getting which jobs are pending etc?
            self.check_completed_jobs()

            # now we need to determine how many jobs are currently running
            # and we need to know how to submit jobs
            # while the number running jobs is smaller than the number of GPUs
            # to use, we submit new jobs
            # indeed it would be optimal if the long running jobs would start first
            # but that is hard to do? For now we just submit the jobs in arbitrary
            # order
            # we know which job is next through the next job variable,
            # but we still need to know the number at which to start
            self.submit_jobs(file)
            # now we just wait some time before we do the next check?
            time.sleep(self.wait)
        # if we are here that means that all runs are done
        # we add the currently simulated temperatures to the bookkeeping variable

class crit_temp_measurement(autonomous_measurement):
    def __init__(self, J_para, J_perp, h, eta, dt, filepath, simulation_path, exec_file, nr_GPUS=6, nr_Ts=5, size_min=48,
                          size_max=80, nr_sizes=3, max_steps=1e9, nr_sites=5e5, Ly_Lx = 1/8, equil_error=0.004,
                 intersection_error=0.02, T_min=None, T_max=None):
        # call the constructor of the parent classe
        super().__init__(J_para, J_perp, h, eta, dt, filepath, simulation_path, exec_file,  nr_GPUS=nr_GPUS, Ly_Lx=Ly_Lx)
        self.nr_Ts = nr_Ts
        self.size_min = size_min
        self.size_max = size_max
        self.nr_sizes = nr_sizes
        self.max_steps = max_steps
        self.nr_sites = nr_sites

        # We can also specify the temperature interval ourselves If we want that
        self.T_min = T_min
        self.T_max = T_max

        self.T_arr = np.array([])
        self.max_T_step = 0.1               # fraction of the critical temperature that is the maximum stepsize for accepted measurements
        self.sizes = np.array([])
        self.max_time = 0
        self.total_runs = 0

        self.all_T_arr = np.array([])       # Bookkeeping for all the temperatures we have simulated in this setting
        self.all_T_dic = {}                 # Bookkeeping dictionary for the different simulation spheres
        self.max_rel_intersec_error = intersection_error  # standard maximum error of 2%
        self.equil_error = equil_error           # standard equilibration error for the U_L runs
        self.maximum_iterations = 4
        self.iteration_nr = 0
        self.repeat = False             # variable that is set to true if we have to repeat a simulation
        self.min_cum_nr = 100000
        self.cum_write_density = 1 / 5
        self.discard_threshold = 0.1     # discards 10% of the U_L values when calculating the mean U_L

        self.cur_para_nr = 0
    def init(self):
        # this somehow needs the parameters, where do we put them? In a file? On the moon? User input?
        T_min = T_c_est(np.abs(self.J_para), np.abs(self.J_perp), self.h)[0]
        # TODO really crude approximation of the maximum T, I really should think about something better
        print(f"T_min = {T_min}\n"
              f"J_perp = {self.J_perp}\n"
              f"h={self.h}")
        # TODO this is still all a bit fishy but...
        T_max = T_min + self.h
        T_max = 10 * T_min
        # or just 10 %?


        if self.T_min is None:
            # If we do not specify the critical temperature, we use the critical temperature estimation
            self.T_min = T_min
        if self.T_max is None:
            self.T_max = T_max
            # If T_max is smaller, equal or slitghly larger than our specified T_min, what do we do then?
            if self.T_max < self.T_min * 1.05:
                # will this happen and how should I deal with it?
                self.T_max = 2 * self.T_min
        # We use nr_Ts datapoints
        self.T_arr = np.linspace(self.T_min, self.T_max, self.nr_Ts)
        # the T_array has to be added to the hierarchy
        self.all_T_dic[0] = self.T_arr  # This first array is on level 0

        self.sizes = np.linspace(self.size_min, self.size_max, self.nr_sizes, endpoint=True,
                            dtype=np.int32)
        print(f"Initializing Simulation with T_min = {T_min} and T_max = {T_max}")
        self.max_time = self.dt * self.max_steps
        self.total_runs = self.nr_sizes * self.nr_Ts  # every temp size combo will be a seperate simulation
    def routine(self):
        """
        the outer routine for the T_c calculation
        :param walltime: walltime per job+3,
        :param file: file to be compiled, we may want to have a seperate file for the autonomous exection
        :param folder: folder for the executables
        :param wait: sleep for wait seconds before checking again if a job is finished        
        :return: -
        """
        self.init()                 # initializing determines the T range
        T_c, T_c_error = self.iteration()  # recursive fuction that returns T_c
        self.conclude()
        return T_c, T_c_error

    def conclude(self):
        # In the new conclude only the simulation with the highest hierarchy is used
        T_arr = self.all_T_dic[np.max(list(self.all_T_dic.keys()))]

        results = self.construct_results(self.discard_threshold, T_arr)

        # interpolate and minimize is deprecated, we use the technique we also use in iteration
        intersections, intersections_y = get_intersections(results)

        T_c = np.mean(intersections)
        U_L_intersection = np.mean(intersections_y)
        T_c_error = np.ptp(intersections)

        self.plotBinderCurve(T_c, U_L_intersection, results)

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

    def plotBinderCurve(self, T_c, U_L_intersection, results):
        fig, ax = plt.subplots(1, 1)
        y_upper_lim = 0
        y_lower_lim = np.infty
        shown_inds = np.linspace(0, len(self.sizes), len(self.sizes) + 1, endpoint=True,
                                 dtype=np.int64)
        ind = 0
        max_T = np.max(self.T_arr) * 1.01
        min_T = np.min(self.T_arr) * 0.99
        for i, size in enumerate(sorted(results.keys())):
            if i in shown_inds:
                T = np.array(results[size]["T"])
                U_L = np.array(results[size]["U_L"])
                ax.plot(T, U_L, linestyle="-", marker="x", color=colors[ind])
                ind += 1
                if max_T:
                    y_upper_lim = np.maximum(
                        np.max(U_L[(min_T <= T) & (T <= max_T)]), y_upper_lim)
                    y_lower_lim = np.minimum(
                        np.min(U_L[(min_T <= T) & (T <= max_T)]), y_lower_lim)
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

    def conclude_old(self):
        # This function should plot stuff etc. keep it simple at this point
        # copied just the cumulanttimeaverage script!
        # has to be rewritten...
        results = {}
        for size_folder in os.listdir(self.simulation_path):
            if (size_folder[0] != ".") & (size_folder != "plots"):
                size_folder_path = os.path.join(self.simulation_path, size_folder)
                if os.path.isdir(size_folder_path):
                    size_result = process_size_folder(size_folder_path,
                                                      self.discard_threshold, selected_temperatures=self.T_arr)
                    results[int(size_folder)] = size_result

        # interpolate and minimize is deprecated, we use the technique we also use in iteration
        intersections = []
        intersections_y = []
        intersections, intersections_y = get_intersections(results)

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

    def iteration(self):
        self.iteration_nr += 1
        # Here I want to have something that checks whether there is already a measurement
        simulation_available = check_directory_structure(self.sizes, self.T_arr, self.simulation_path)
        print("simulation_available", simulation_available)
        if not simulation_available or self.repeat:
            self.repeat = False                 # we set repeat to false again
            self.cur_para_nr = 0                # reset the parameter number
            self.write_para_files()  # setting up the parameter files for every simulation
            self.cur_para_nr = 0                # reset the parameter number
            self.run_jobs()
        else:
            print("Found valid simulation, evaluating")
        # after running the jobs, we need to
        # calculate the binder cumulant for every run
        return self.evaluate_simulation()

    def evaluate_simulation(self):
        # The simulation we are at right now is the one with the highest key
        sim_hierarchy_nr = np.max(list(self.all_T_dic.keys()))
        T_arr = self.all_T_dic[sim_hierarchy_nr]
        # for this one we want to check if it has an intersection
        results = self.construct_results(self.discard_threshold, T_arr)

        U_L_min_T_min = np.min(results[np.min(self.sizes)]['U_L'])
        U_L_max_T_min = np.min(results[np.max(self.sizes)]['U_L'])
        U_L_min_T_max = np.max(results[np.min(self.sizes)]['U_L'])
        U_L_max_T_max = np.max(results[np.max(self.sizes)]['U_L'])
        # we say we have an intersection if U_L_min_T_min > U_L_max_T_min
        # and U_L_min_T_max < U_L_max_T_max
        # This is we have an intersection at all, but we dont know if we have an intersection in the current T_simulation.
        intersection = ((U_L_min_T_min >= U_L_max_T_min) & (
                U_L_min_T_max <= U_L_max_T_max)) or (U_L_max_T_max / U_L_max_T_min > 2.95)       # TODO this is a bit fishy but will probably work in 99% of times, catches the case that the maximum temperature is far in the high temperaturef region and therefore inflicated with strong fluctiations
        dT = T_arr[1] - T_arr[0]
        if intersection:
            # now the usual stuff, estimate the postion of the intersection
            intersections, intersections_y = get_intersections(results)

            T_c = np.mean(intersections)

            print(f"Found an intersection at T_c = {T_c}")

            print("intersections: ", intersections)

            T_c_error = np.ptp(intersections)
            # the error will now be the minimum of this and a fourth of the stepsize
            T_c_error = max(T_c_error, dT / 5)
            print(f"T_c_error = {T_c_error}")
            rel_intersec_error = T_c_error / T_c
            print(f"rel_intersec_error = {rel_intersec_error}")
            if rel_intersec_error < self.max_rel_intersec_error:
                # In this case we are done since we covered the case of large stepsizes, small errors with the dT / 4
                print(f"Determined crit. Temp T_c = {T_c} +- {rel_intersec_error}")
                return T_c, T_c_error
            else:
                # We still want to plot the stuff to see where we at
                U_L_intersection = np.mean(intersections_y)
                self.plotBinderCurve(T_c, U_L_intersection, results)
                # If the error is too large we do a child measurement with smaller error
                self.equil_error /= 2
                # We wanted to have our edgecases of within 5 or 20 percent of the interval edges...
                T_interval_low = np.max(T_arr[T_arr <= T_c])
                T_interval_up = np.min(T_arr[T_arr >= T_c]) # this should hopefully guaranteed to be always dT larger than lower interval?
                if T_c < (T_interval_low + 0.05 * dT):
                    # in this case we want to half the new interval
                    T_interval_up = T_interval_low + dT / 2
                    T_interval_low -= 0.02 * dT
                elif T_c > (T_interval_up - 0.2 * dT):
                    T_interval_low = T_interval_up - dT / 2
                    T_interval_up += 0.02 * dT
                # else we just take the whole interval
                # We actually dont want to have exactly the same temperature as in the previous run
                # we think we are fairly sure to include the critical point if we subtract 0.01dT from the interval edges
                T_min = T_interval_low + 0.01 * dT
                T_max = T_interval_up - 0.01 * dT
                self.T_arr = np.linspace(T_min, T_max, self.nr_Ts)
                # If we are here this directly means that we started a new child a a new level in the hierarchy
                self.all_T_dic[sim_hierarchy_nr + 1] = self.T_arr

                print(f"Error was too large: Temp T_c = {T_c} +- {T_c_error} \n"
                      f"Starting new run with T_min = {T_min}, T_max = {T_max}")
                return self.iteration()
        else:
            # This means we missed the intersection so we want to restart a simulation with the same error
            print("We do not see a intersection")
            if U_L_min_T_max > U_L_max_T_max:
                print("The maximum temperature is too low")
                # this means we are below the critical temperature.
                # TODO can we somehow approximate how far away we are?
                # for now we just double the range i would say?
                T_range = np.ptp(self.T_arr)
                T_min = np.max(self.T_arr) + (dT)  # the next T_min is the current maximum plus the current stepsize
                T_max = np.max(self.T_arr) + T_range
                self.T_arr = np.linspace(T_min, T_max, self.nr_Ts)
                # Okay we updated the temperature array... now we just run everything again?
            elif U_L_min_T_min < U_L_max_T_min:
                print("The minimum temperature is too high")
                T_range = np.ptp(self.all_T_arr)
                T_min = np.maximum(np.min(self.T_arr) - T_range, 0.0)  # We should not consider negative temperatures
                T_max = np.min(self.T_arr) - (self.T_arr[1] - self.T_arr[0])
                self.T_arr = np.linspace(T_min, T_max, self.nr_Ts)
            # we add the new stuff to the dictionary, but we already overwrote the T_arr, but the T_arr should still be
            # fine
            self.all_T_dic[sim_hierarchy_nr] = np.concatenate((T_arr, self.T_arr))

            return self.iteration()

    def evaluate_simulation_old(self):
        self.all_T_arr = np.concatenate((self.all_T_arr, self.T_arr))   # we do it here, before evaluating
        all_results = self.construct_results(self.discard_threshold, self.all_T_arr)     # for all results we use the self.all_T_arr? Since we want to reproduce always the measurements that we already did so that we can resume it where we left it of
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
        # This is we have an intersection at all, but we dont know if we have an intersection in the current T_simulation.
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
            intersections, intersections_y = get_intersections(all_results)
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
                # If the T stepsize is to large, lets say larger than 10 % of the critical temperature, we repeat the
                # simulation anyway
                T_stepsize = self.T_arr[1] - self.T_arr[0]
                if T_stepsize > 0.1 * T_c:
                    # so now again a measurment +- 5% of the critical temperature?
                    self.T_arr = np.linspace(0.95 * T_c, 1.05 * T_c, self.nr_Ts)
                    self.equil_error /= 2
                    return self.iteration()
                else:
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
                return self.iteration()
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
            return self.iteration()

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

    def get_para_nr(self):
        # this tactic inreases the parameter number everytime get_para_nr is called so that we do not submit any job twice
        self.cur_para_nr += 1
        return self.para_nr + self.cur_para_nr - 1

    def write_para_files(self):
        # you ..., you know that you have to construct the parameter file at hemera?
        # and you need to do rsync after the jobs are finished!
        print("Writing the parameter files...")
        for i, (T, size) in enumerate(product(self.T_arr, self.sizes)):
            # We now need to construct the parameterfile with the appropriate temperature and size
            # Is it okay if we construct all files in the beginning and deal with the threading of the gpus later?
            # to construct the para set we need to know how many subsystems we should initialize
            nr_subsystems = int(self.nr_sites / (size ** 2 * self.Ly_Lx))
            with open(self.filepath + "/parameters/para_set_" + str(self.para_nr + i) + '.txt', 'w') as f:
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
                        f"equil_error, {self.equil_error}\n"
                        f"cum_write_density, {self.cum_write_density}\n"
                        f"min_cum_nr, {self.min_cum_nr}")
        # we need to copy the files to hemera
        rsync_command = ["rsync", "-auv", "--rsh", "ssh",
                         f"{self.filepath}/parameters/",
                         "hemera:~/Code/Master-Arbeit/CudaProject/parameters/"]
        subprocess.run(rsync_command, cwd=pathlib.Path.home())

class quench_measurement(autonomous_measurement):
    def __init__(self, J_para, J_perp, h, eta, dt, filepath, simulation_path, exec_file, Tc, nr_GPUS=6, size_min=64,
                 size_max=4096, nr_sites=5e5, Ly_Lx=1/8, min_quench_steps=100, min_nr_sites=1e6,
                 min_nr_systems=10, host="hemera", user="weitze73"):
        super().__init__(J_para, J_perp, h, eta, dt, filepath, simulation_path, exec_file,  nr_GPUS=nr_GPUS, Ly_Lx=Ly_Lx)
        self.size_min = size_min        # The starting size at which we do go higher
        self.size_max = size_max        # maximum size, if xi = size_max / 10 we stop the simulation
        self.nr_sites = nr_sites        # nr of sites for one run
        self.Tc = Tc                    # The critical temperature that was calculated in the measurement before
        self.min_quench_steps = min_quench_steps    # the minimum number of steps done during a quench, influences tau_min
        self.min_nr_sites = min_nr_sites            # the minimum nr of sites i want to simulate for every tau
        self.min_nr_systems = min_nr_systems        # the minimum nr of systems to go through the qunech, guaranteeing an approximately precice xi value

        self.T_start = Tc
        self.T_end = Tc
        self.tau_min = 1
        self.tau_factor = 1                    # current tau
        self.tau_list = []                 # empty list to keep track of the taus that we used
        self.size = size_min            # current size, starts at size_min
        self.equil_error = 0.02        # the equilibration error, so the error of U_L for which we assume that we are approximately equilibrated, doesnt need to be as small as for the T_c measurement
        self.para_nr = 100

        self.cut_zero_impuls = True     # will probably always be true since we are quenching
        self.fitfunc = lorentz_offset      # We usually use this atm?
        self.nr_measured_values = 300   # standard number of measured values during the quench
        self.min_tau_scaling_fit = 10

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
        return int(np.ceil(self.nr_sites / (self.size ** 2 * self.Ly_Lx)))

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
        sim_available = check_directory_structure([self.size], [self.tau()], self.simulation_path)      # the check directory structure takes lists in
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
        self.total_runs = nr_jobs
        # now we do basically the same stuff as last time, only that we onyl submit the same job
        super().run_jobs()

    def get_para_nr(self):
        return self.para_nr
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

class amplitude_measurement(autonomous_measurement):
    def __init__(self, J_para, J_perp, h, eta, dt, filepath, simulation_path, exec_file, Tc, nr_GPUS=6, nr_Ts=6, size=1024,
                 max_steps=1e9, Ly_Lx = 1/8, equil_error=0.001, equil_cutoff=0.1, T_range_fraction=0.05):
        super().__init__(J_para, J_perp, h, eta, dt, filepath, simulation_path, exec_file,  nr_GPUS=nr_GPUS, Ly_Lx=Ly_Lx)

        self.nr_Ts = nr_Ts                      # nr of temperatures used to fit
        self.T_range_fraction = T_range_fraction    # This is the fraction of Tc that is used to determine the interval of Ts in [Tc, (1+T_range_raction) * Tc]
        self.Tc = Tc                            # We need to know the critical temperature that we determined in previous measurements
        # I think we only use one size, being the largest possible
        # Is there some way of guessing how long a simulation will take?
        # 1000 is a save bet but I guess we would like to go up to even larger sizes?
        self.size = size
        self.max_steps = max_steps

        self.T_arr = np.array([])
        self.total_runs = 0                 # just a parameter to keep track of how many jobs we want to run in this iteration


        self.all_T_arr = np.array([])       # Bookkeeping for all the temperatures we have simulated in this setting
        self.equil_error = equil_error           # standard equilibration error for the xi runs
        self.maximum_iterations = 2         # we first look at the interval [Tc, 1.05Tc] and If this doesnt work we inrease to [Tc, 1.1Tc] and If this doesnt work we abort
        self.iteration_nr = 0
        self.min_corr_nr = 5000
        self.corr_write_density = 1 / 10
        self.equil_cutoff = equil_cutoff             # This is the values that we cut off because we think we are still equilibrating. Since we definitely want the values in equilibration we use a relatively large cutoff here
        self.max_time = 0
        self.Tc_fit_tolerance = 0.025        # 5% tolerance for the Tc obtained from the linear regression around the critical point. If its further away, we do not accept the fit
        self.min_r_sqaured = 0.98           # The r²-value of the linear regression should be fairly high so that we can be sure that the interval that we fit is really linear

        self.para_nr = 110                 # A different para nr for this simulation?
        self.cur_para_nr = 0                # same method as in Tc?
    def setup(self):
        # This function will determine the initial T range
        # we want nr_Ts temperatures between
        Tmax = self.Tc * (1 + self.T_range_fraction)
        self.T_arr = np.linspace(self.Tc, Tmax, num=self.nr_Ts)
        # Is that all, what else do we need to setup?
        # For now the total number of runs is just the nr of temps
        self.total_runs = self.nr_Ts
        self.max_time = self.dt * self.max_steps
    def run(self):
        # runs the complete simulation. Some initialization, a routine and a finish, plotting and fitting
        self.setup()
        return self.iteration()     # In this suite  the conclude is called directly in iteration
    def iteration(self):
        # Okay again this is the function that runs the jobs, evaluates them and runs additional jobs if necessary. Also this
        # this function should be able to pick up on simulations that were not finished
        # increase the iteration number
        self.iteration_nr += 1
        # We add the temperatures from now to all temperatures
        self.all_T_arr = np.concatenate((self.all_T_arr, self.T_arr))
        # Check if we have a simulation available
        # I think we should improve the pickup capability, if we run for example half the jobs of simulation we should be able to use them
        valid_simulations = get_avail_simulations([self.size], self.T_arr, self.simulation_path, check_corr_valid,
                                                  check_function_args=(self.equil_error, self.equil_cutoff))
        # For every valid simulation we do not have to do this simulation in the following
        for valid_simulation in valid_simulations:
            valid_temp = valid_simulation[1]        # valid simulations is tuple of (size, temp), we only need the temp here
            Ts = list(self.T_arr)                   # change to list as we then can easier use remove
            Ts.remove(valid_temp)                   # remove the valid T from the Ts we still have to do
            self.T_arr = np.array(Ts)               # overwrite the original array with the new Ts

        if self.T_arr.size != 0:
            # if the array is not empty this means that there are still simulations to do
            # write the parameter files
            # I also need to reset this current parameter number here...
            self.cur_para_nr = 0
            self.write_para_files()
            # submit the jobs, should be handlede by the super method
            self.cur_para_nr = 0
            self.run_jobs()
            # If the array is empty this means that we have all measurements done already and we can call evaluate
        # We want to evaluate anyways if we have a valid simulation or if we just ran one
        return self.evaluate()
    def write_para_files(self):
        print("Writing the parameter files...")
        for i, T in enumerate(self.T_arr):
            # We now need to construct the parameterfile with the appropriate temperature
            # Is it okay if we construct all files in the beginning and deal with the threading of the gpus later?
            with open(self.filepath + "/parameters/para_set_" + str(self.get_para_nr()) + '.txt', 'w') as f:
                f.write(self.simulation_path)
                f.write(f"\nend_time, {self.max_time} \n"
                        f"dt, {self.dt} \n"
                        f"J, {self.J_para} \n"
                        f"Jy, {self.J_perp} \n"
                        f"alpha, {self.h} \n"
                        f"eta, {self.eta} \n"
                        f"nr_saves, 4 \n"           # We dont know how long the simulation will go so we could either use a density observer or weeeee just dont care
                        f"nr_repeat, 0 \n"
                        f"min_temp, {T} \n"
                        f"max_temp, {T} \n"
                        f"nr_runs, 0.0 \n"
                        f"random_init, 1.0 \n"      # for the amplitude we want to initialize randomly
                        f"curand_random, 1 \n"
                        f"subsystem_min_Lx, {self.size} \n"
                        f"subsystem_max_Lx, {self.size} \n"
                        f"nr_subsystem_sizes, 0  \n"
                        f"nr_subsystems, {1} \n"    # The number of subsystems will be one, we use large systems that will run long to eliminate the statistical deivations
                        f"x_y_factor, {self.Ly_Lx} \n"
                        f"nr_corr_values, 0 \n"     # We need a new corr observer that just observes with density and doesnt switch after quench     
                        f"nr_ft_values, 0 \n"       # Ah we still wanted to check whether the values of the ft and fit and python or direct fit in c++ are the same, but they should be fairly similar
                        f"equil_error, {self.equil_error}\n"
                        f"equil_cutoff, {self.equil_cutoff}\n"
                        f"min_corr_nr, {self.min_corr_nr}\n"
                        f"corr_write_density, {self.corr_write_density}\n")
        # we need to copy the files to hemera
        rsync_command = ["rsync", "-auv", "--rsh", "ssh",
                         f"{self.filepath}/parameters/",
                         "hemera:~/Code/Master-Arbeit/CudaProject/parameters/"]
        subprocess.run(rsync_command, cwd=pathlib.Path.home())

    def evaluate(self):
        # Okay the evaluation logic, will be done after lunch, you will just extract the xi, fitting and see if its okay
        # What do we need to to, get right into it, just turn on some music and focus
        # we know the simulation path and the temperatures that we just simulated, but we actually want to evaluate all
        # temperatures that we have
        # Here we can call the process size folder method
        size_path = f"{self.simulation_path}/{self.size}"
        xix_dic = process_size_folder(size_path, threshold=self.equil_cutoff, key="T", value="xix",
                                  file_ending="corr")      # the selected temperatures are just all temperatures
        xiy_dic = process_size_folder(size_path, threshold=self.equil_cutoff, key="T", value="xiy",
                                  file_ending="corr")
        # Those dictionaries map the temperature to the according xi, ahh not quite those are dictionaries
        # with xix_dic = {"T" : [T1, T2, ...], "xix" : [xix1, xix2, ...]}
        T_xix = xix_dic["T"]
        xix_arr = xix_dic["xix"]
        T_xiy = xiy_dic["T"]
        xiy_arr = xiy_dic["xiy"]

        # They should be sorted and we are actually interested in the inverse correlation lengths
        xix_arr = xix_arr[np.argsort(T_xix)]
        xiy_arr = xiy_arr[np.argsort(T_xiy)]
        xix_inv = 1 / xix_arr
        xiy_inv = 1 / xiy_arr
        T_xix = T_xix[np.argsort(T_xix)]
        T_xiy = T_xiy[np.argsort(T_xiy)]
        # Now we have the inverse arrays for the two directions.
        # For the linear regression we need the critical temperature which we can access by self.Tc
        # so we can call the best_fit_inv function that we wrote for the evaluation outside of the suite
        reg_x, T_include_start_x, T_include_end_x = best_fit_inv(T_xix, xix_inv, self.Tc, self.Tc_fit_tolerance, self.min_r_sqaured)
        reg_y, T_include_start_y, T_include_end_y = best_fit_inv(T_xiy, xiy_inv, self.Tc, self.Tc_fit_tolerance,
                                                                 self.min_r_sqaured)
        # The best fit_inv function returns None if we didnt find a fitting segment
        if (reg_x is None) or (reg_y is None):
            print("No fitting fit.")
            # if we are doing to many iterations we have to return here
            if self.iteration_nr > self.maximum_iterations:
                print("Maximum iterations reached. Aborting")
                return
            else:
                print("We have to repeat the simulation")
            # so we extend the temperatures that we investigate
            stepsize = self.T_arr[1] - self.T_arr[0]
            self.T_arr = np.linspace(np.max(self.T_arr) + stepsize, (1 + 2 * self.T_range_fraction) * self.Tc)
            # Just iteration now?
            return self.iteration()
        else:
            # means borh are not none so we found fits
            # If we found something we can call conclude directly from here where we have the fits and the start and endpoints in scope?
            print("Fitting worked, concluding...")
            # We have to hand over pretty much data seems like
            x_data = (T_xix, xix_arr, xix_inv)
            y_data = (T_xiy, xiy_arr, xiy_inv)
            x_result = (reg_x, T_include_start_x, T_include_end_x)
            y_result = (reg_y, T_include_start_y, T_include_end_y)
            self.conclude(x_data, y_data, x_result, y_result)
    def conclude(self, x_data, y_data, x_result, y_result):
        # This function is supposed to plot the fits that we did alongside with the data and present the amplitude that
        # we calculated
        # So we use the fits and all the data to do some plots
        # Do we want to plot into the same plot or in different plots?
        # If we want to plot in the same plot we need to use logarithmic axes. Otherwise one direction will not be easily recognizable
        # If we use logarithmic axes, we would not have the usual critical divergence form or respectively the linear scaling in 1 / xi
        # Or i mean it will be linear but every other polynomial would also be linear
        # what will it be... decide! I see them in a single log plot tbh
        # first we will calculate the critical temperatures of the two cases aswell as the critical amplitudes
        reg_x = x_result[0]
        reg_y = y_result[0]
        T_x = x_data[0]
        T_y = y_data[0]
        xix_ampl = 1 / reg_x.slope
        xiy_ampl = 1 / reg_y.slope
        Tc_x = - reg_x.intercept * xix_ampl
        Tc_y = - reg_y.intercept * xiy_ampl

        # With this we can just plot the stuff? first the inverse stuff
        fig, ax = plt.subplots(1, 1)
        # First plot the data points
        ax.plot(T_x, x_data[2], label=rf"$1 / \xi_x$", linestyle="", marker="x")
        ax.plot(T_y, y_data[2], label=rf"$1 / \xi_y$", linestyle="", marker="x")
        # Now plot the fits
        ax.plot(T_x, reg_x.intercept + reg_x.slope * T_x,
                 label=rf"$\xi_x^+ = {xix_ampl:.2f}, T_c = {Tc_x:.3f}$")
        ax.plot(T_y, reg_y.intercept + reg_y.slope * T_y,
                label=rf"$\xi_y^+ = {xiy_ampl:.2f}, T_c = {Tc_y:.3f}$")
        # We want to show the ratio on the plot
        ax.plot([], [], label=rf"$\xi_x / \xi_y  = {xix_ampl / xiy_ampl}$", linestyle="", marker="")
        # like I said if we want them to be in one plot we need to scale the y axix to be logarthimic, you cant do that, the fits will cross zero and so wont look good logarithmic
        ax.set_xlabel("T")
        ax.set_ylabel(r"$1 / \xi$")
        # We set the lower limit of the y axis to be zero since negative xi are not sensible but the fit can become negative
        ax.set_ylim(0, ax.get_ylim()[1])
        configure_ax(fig, ax)
        # saving the plot
        create_directory_if_not_exists(self.simulation_path + "/plots/")
        plt.savefig(self.simulation_path + "/plots/xi-inv.png", format="png")
        plt.show()

        # We also want to show the divergence and include the actual xi values
        fig, ax = plt.subplots(1, 1)
        # First plot the data points
        ax.plot(T_x, x_data[1], label=rf"$1 / \xi_x$", linestyle="", marker="x")
        ax.plot(T_y, y_data[1], label=rf"$1 / \xi_y$", linestyle="", marker="x")
        # We also want to plot some kind of fit but for this we need the eps arrays
        # We need to use the critical temperature of the fit
        # I think we want some more points than 8 to plot the critical amplitude plot
        T_x_plot = np.linspace(np.min(T_x), np.max(T_x), 200)
        T_y_plot = np.linspace(np.min(T_y), np.max(T_y), 200)
        eps_x = (T_x_plot - Tc_x) / Tc_x
        eps_y = (T_y_plot - Tc_y) / Tc_y
        # before we plot we look at the y limits in the case that we dont plot the critical amplitude
        upper_ylim = ax.get_ylim()[1]
        # The function that we need to use is called critical amplitude
        ax.plot(T_x_plot, critical_amplitude(eps_x, xix_ampl),
                 label=rf"$\xi_x^+ = {xix_ampl:.2f}, T_c = {Tc_x:.3f}$")
        ax.plot(T_y_plot, critical_amplitude(eps_y, xiy_ampl),
                label=rf"$\xi_y^+ = {xiy_ampl:.2f}, T_c = {Tc_y:.3f}$")
        # For this we dont use logarithmic scale I think
        # we set the limits from before plotting the fit
        ax.set_ylim(0, upper_ylim)
        ax.set_xlabel("T")
        ax.set_ylabel(r"$\xi$")
        configure_ax(fig, ax)
        # Save the plot
        plt.savefig(self.simulation_path + "/plots/T-xi.png", format="png")
        plt.show()

        # I think that is it.
        return
    def get_para_nr(self):
        # this tactic inreases the parameter number everytime get_para_nr is called so that we do not submit any job twice
        self.cur_para_nr += 1
        return self.para_nr + self.cur_para_nr - 1

class z_measurement(autonomous_measurement):
    def __init__(self, J_para, J_perp, h, eta, dt, filepath, simulation_path, exec_file, test_exec_file, Tc, nr_GPUS=6, size_min=64,
                          size_max=256, nr_sizes=3, max_steps=1e9, nr_sites=5e5, Ly_Lx = 1/8, equil_error=0.005, equil_cutoff=0.5,
                 variation_error_rate=0.011, z_guess=2, min_nr_sites=1e6, min_nr_systems=100, fold=50, cum_density=1/100,
                 test_cum_density=1/2, test_min_cum_nr=2000):
        # call the constructor of the parent classe
        super().__init__(J_para, J_perp, h, eta, dt, filepath, simulation_path, exec_file,  nr_GPUS=nr_GPUS, Ly_Lx=Ly_Lx)
        # ah maybe you should first write a concept about what this measurment should do

        self.Tc = Tc
        self.size_min = size_min
        self.size_max = size_max
        self.nr_sizes = nr_sizes
        self.max_steps = max_steps
        self.nr_sites = nr_sites
        self.equil_error = equil_error      # only useful for the first run that determines how long the runs will be
        self.test_exec_file = test_exec_file

        # What else do we need?
        self.variation_error_rate = variation_error_rate      # this will be the error of the standard deviation of the differences of the U_L values
        self.z_guess = z_guess                      # this is the guess for z that we will need to estimate how long the runs for the larger sizes should run
        self.test_size = size_min / 2               # the test size to see how long the simulations should run will be half of the minimum size
        # for this run we also need the equil error
        self.test_equil_time = 0                          # this is the equilibration time of the test system which will later be used to guess the simulation times for the other sizes

        self.min_nr_sites = min_nr_sites
        self.min_nr_systems = min_nr_systems        # same system like for quench, the U_L values that we use should be averaged out of at least 100 systems and if 100 systmes contain less than min nr of sites then so many systems that we reach teh min nr of sites
        self.nr_points_diff_variance = 20           # this is the number of consecutive datapoints from which we will calculate the differences and then their standard deviation to judge how strongly the data is variation
        # We could actually also think about just averaging all U_L values, this should already give a good measure
        self.cum_density = cum_density                  # the standard density will write 1 cum value in 100 steps
        self.test_cum_density = test_cum_density        # the test cum density will be much larger since we want to see if we are equilibrated as fast as possible
        # Okay now we basically do the same stuff as the last 3 times
        self.sizes = np.array([])                   # this is the array for sizes that we have to simulate
        self.valid_sizes = []                         # this is the list to keep track of the sizes that we simulated
        self.cur_para_nr = 0                        # since we need multiple paramter files we need the current parameter variable
        self.cur_run_nr = 0                         # we will use this one in the run_jobs logic to know at which submit we are and which parameter number we have to return
        self.para_nr_run_dic = {}                   # this is the dictionary that will know which run number has which parameter number
        self.para_nr = 200                          # again different parameter number for different kind of simulation
        self.base_fold = fold                       # for plotting we fold 10 points to one for better visibility

        self.test_min_cum_nr = test_min_cum_nr      # the minimum nr of cumulant values we calculate for the testrun
        self.test_equil_cutoff = equil_cutoff       # the values we cut off when running the test run


    def setup(self):
        # This function does the setup again, I mean we dont have to many but
        self.sizes = 2 ** np.linspace(np.log2(self.size_min), np.log2(self.size_max), self.nr_sizes, dtype=np.int64)
        # what else do we need to do?
        # for what did I need total runs again? Just for submitting jobs in this run right?
        # how many jobs do we have to submit to get the minimum of nr_sites? For the testrun it might be okay to
        # only submit one graphics card with the
        self.total_runs = 1                     # the number of total runs is not only the number of sizes sadly
        # the number of jobs for the first run_jobs will just be one as it will be the test measurement to get the equilibration time
        # do we want to construct the dict that knows the parameter numbers here?
        self.para_nr_run_dic[0] = 0     # we should not need more

    def iteration(self):
        # so here again we need to implement the pickup capability and stuff
        # first of all if this is the first iteration (or should we call that in the setup already)?
        # we could just do the test measurement every iteration, or check if the test measurement is there already so we only do it once
        # we could test if the test_equil_time is zero, if yes we look if we already have a test measurment from an earlier run and if we really do not see anything we really run it
        if self.test_equil_time == 0:
            # now we check if the simulation is available
            # I dont think that we have to implement a checkfunction since it doesnt really matter if the variation of the differences is small since we only want to guess
            # the time
            test_measurement_available = check_directory_structure([self.test_size], [self.Tc], self.simulation_path)
            if not test_measurement_available:
                # If else we need to perform the simulation
                # so we need to write parameters
                self.cur_para_nr = 0
                self.write_test_paras()
                # If we did this we reset the cur para nr?
                self.cur_para_nr = 0
                # now we need to run those the jobs
                self.run_jobs(self.test_exec_file)     # this should work like this if we have the correct current para nr and the nr of total runs
                # no it doesnt work quite like this since we need to combine multiple parameter files with multiple runs for the same parameter value...
                # sooo. what are we going to do ? I guess we need to implement some more elaborate get_para_nr logic.
                # the thing is we need different behaviors depending on if we get the parameter number to write the parameter file or if we want to
                # submit the job
                # as soon as the rn jobs is done we can extract the test_equil_time
            # if this is available, we read the parameters? or some other file and extract the last time. Sadly I dont write the quilibration time to a file so we just read the .cum file?
            # We need the filepath
            self.test_equil_time = self.get_last_time()
        # Now we have an equil time
        # This means we already have a test_equil_time
        # so we check if the simulations that are supposed to run have already run
        # we also should check if we trust the result, so if the variation error is small enough
        # this means we have to write another validation function, done
        # so now we do the same thing as in z_extraction suite, we let the computer search the valid measurements
        valid_simulations = get_avail_simulations(self.sizes, [self.Tc],
                                                  self.simulation_path,
                                                  check_cum_valid,
                                                  check_function_args=(
                                                  self.variation_error_rate,
                                                  self.nr_points_diff_variance,
                                                  self.dt / self.cum_density))
        # every valid simulation does not have to be repeated
        for valid_simulation in valid_simulations:
            # brutal but effective i guess
            valid_size = valid_simulation[0]
            sizes = list(self.sizes)
            sizes.remove(valid_size)
            self.sizes = np.array(sizes)
            # if the size is valid we can add it to the sizes array
            self.valid_sizes.append(valid_size)

        if self.sizes.size != 0:
            # this means there are still measurements to do
            # write the parameters
            self.write_para_files()
            # I guess now is the time to construct the 'para_nr_run_dic' ?
            self.para_nr_run_dic = {}       # but first we reset it
            self.total_runs = 0
            self.construct_para_nr_run_dic()    # this function also has to determine the total nr of jobs
            # now we run the jobs... I guess?
            self.run_jobs()
            self.cur_run_nr = 0
            #... we dont know whether the simulations are valid atm so we dont add them to valid_sizes?
        return self.evaluate()

    def evaluate(self):
        # so this one looks at a done measurement and decides whether if was okay or not
        # tbh to check if it is valid was already done by the get valid_measurements function, so we should be able to reuse that here?
        # what we have to write is the plotting and fitting...
        valid_simulations = get_avail_simulations(list(self.sizes) + self.valid_sizes, [self.Tc],
                                                  self.simulation_path,
                                                  check_cum_valid,
                                                  check_function_args=(
                                                  self.variation_error_rate,
                                                  self.nr_points_diff_variance,
                                                  self.dt / self.cum_density))
        if len(valid_simulations) == self.nr_sizes:
            # this means the simulation was valid and we can do the rescaling
            print("Simulation successful, evaluating")
            # add to valid sizes
            # [self.valid_sizes.append(valid_sim[0]) for valid_sim in valid_simulations]          # is this proper python to shorten the for loop? Probably a bit slower than a normal for loop but this doesnt matter here
            # so we had already a method to average the cumulant in a folder, right?
            # In the CumulantOverTime script we use two maps, one that maps the size and the temperatue to the cum array
            # and another that maps it to the time
            # We shouldnt have different temperatures in this case, but what if we do because we adjusted the critical temperature?
            # We write something that takes a simulation folder and iterates through every simulation for like the 100th time,
            # shouldnt you somehow encapsulate this? F it for now, this could be complicated
            size_cum_dic, size_times_dic = self.get_results_time_resolved()
            # Okay we have the values, now we need to do the remapping
            # we need a list of z's that we want to check, just user input parameter? Nah I think we can do it like this. We will just start at 1 and go to three or something like this?
            # It wont be outside this right? Afterwards we can for every size pair get a new z interval
            z_list = np.linspace(1, 3, 201)
            # we want to run sorted through the keys
            self.valid_sizes = sorted(self.valid_sizes)
            # It is the same thing to map the small sizes on the large sizes and vice versa, right?
            fig, ax = plt.subplots(1, 1)
            for i, size in enumerate(self.valid_sizes):
                # We prbably want a 'base fold' for the smallest size
                cum = size_cum_dic[size]
                times = size_times_dic[size]

                t_fold, cum_fold = fold(times, cum, fold=self.base_fold)
                # plot the points
                ax.plot(t_fold, cum_fold, linestyle='', marker="x",
                        markersize=5,
                        color=f"C{i}", label=f"$L_x$={size},  T = {self.Tc}")
                # the new interploation works with np.interp
                # the number of points of interp should be... i dont know at least a bit larger than the number of folded values?
                # what even happens when you try to use less values haha?
                # we could just use base_fold * len(t_fold) to get the same number of values we had before but folded?
                # The thing is for the error calculation we need the same interval for both sizes, I guess here it doesnt matter to much
                t_inter_plot = np.linspace(np.min(t_fold), np.max(t_fold), self.base_fold * len(t_fold))
                cum_inter_plot = np.interp(t_inter_plot, t_fold, cum_fold)

                ax.plot(t_inter_plot, cum_inter_plot, f"C{i}")
                # okay so much for the plotting, now the refitting
                # if there is a larger size we want to map this size onto it
                if i < len(self.valid_sizes) - 1:
                    next_size = self.valid_sizes[i + 1]
                    # cumulant and time values of the next size:
                    cum_next = size_cum_dic[next_size]
                    times_next = size_times_dic[next_size]
                    # the folding shall be calculated based on the z that we use? so that t_fold and t_fold_next have
                    # approximately the same number of points in the shared interval
                    b = next_size / size            # next size is the larger size so b > 0
                    # values to keep track wich z is the best
                    best_z = 0
                    best_msd = np.infty
                    best_t_compare_arr = []
                    best_cum_compare_arr = []
                    # we need to try out every z
                    for z in z_list:
                        # first we resacle the time of the small size?
                        times_next_rescaled = times_next / (b ** z)
                        # now we can look vor the shared interval by times_next and times_rescaled
                        t_lower_limit = np.maximum(np.min(times_next_rescaled), np.min(times_next))      # the lower limit is the greater minimum of the times
                        t_upper_limit = np.minimum(np.max(times_next_rescaled), np.max(times_next))

                        # between those bounds we need an array with dense t values to compare the two cum curves
                        # the number of values should be... i dont know how large, just large enough? The arrays we will use
                        # to calculate the interpolation will use folded values, so just len(times_next) or something like this?
                        t_compare_arr = np.linspace(t_lower_limit, t_upper_limit, len(times_next))
                        # now we need the interpolations, but we want to use the interpolation of folded values, so we fold?
                        # so now the folding of the previous will stay the same, but the folding of the next will change
                        # if we for example rescale by a factor of 4, the folding should be 4 times as large?
                        t_next_rescaled_fold, cum_next_fold = fold(times_next_rescaled, cum_next, b**z * self.base_fold)              # one with just the interpolation of the next curve without rescaling (Constant and independent of z)
                        # oh maaaaaaan you did it wrong you wanted to rescale the larger size because it has more values which can be folded to get fewer fluctuations
                        cum_next_compare_arr = np.interp(t_compare_arr, t_next_rescaled_fold, cum_next_fold)
                        cum_compare_arr = np.interp(t_compare_arr, t_fold, cum_fold)
                        # Okay from them we can now calculate an error?
                        err_arr = cum_next_compare_arr - cum_compare_arr
                        msd = np.mean(err_arr ** 2)
                        if msd < best_msd:
                            best_z = z
                            best_msd = msd
                            best_t_compare_arr = t_compare_arr
                            best_cum_compare_arr = cum_next_compare_arr
                # If we did this we want to plot it
                    ax.plot(best_t_compare_arr,
                            best_cum_compare_arr,
                            linestyle="-", label=f"$L_x$={next_size},"
                                                 f"  T = {self.Tc}  rescaled z = {best_z:.3f}", color=f"C{i+1}")
            ax.set_ylabel(r"$U_L$")
            ax.set_xlabel("t")
            ax.set_xlim(0, ax.get_xlim()[1] / 4)
            configure_ax(fig, ax)
            ax.set_title("z extraction")
            plt.savefig(self.simulation_path + f"/cum-over-time-scan.png", format="png")
            plt.show()
            return
        else:
            # else means we need some logic to see which measurement has to be repeated
            print("Error of ... to large, repeating")
            # So we need to know which simulations we have to repeat
            # we just remove the valid simulations from the self.sizes and run the iteration again
            for valid_sim in valid_simulations:
                valid_size = valid_sim[0]
                sizes = list(self.sizes)
                try:
                    # this should be okay, it tries to remove all valid sizes but some of them might already have been valid
                    # so the try only removes the ones that are newly valid
                    sizes.remove(valid_size)
                except ValueError:
                    pass
                self.sizes = np.array(sizes)
            # now just rerun=?
            return self.iteration()

    def get_results_time_resolved(self):
        size_cum_dic = {}
        size_times_dic = {}
        for size in self.valid_sizes:
            size_path = os.path.join(self.simulation_path, str(size))
            if os.path.isdir(size_path):
                temp_path = os.path.join(size_path, f"{self.Tc:.6f}")  # we only want to look at the current used critical temperature
                if os.path.isdir(temp_path):
                    # okay so here we extract the cumulant for the a certain size temp combination
                    # I guess we use the same tactic as in the ccumovertime script since I think there is no simpler
                    # method to extract both the time and the values
                    # we can use get folder avage here
                    cum, times = get_folder_average(temp_path)
                    size_cum_dic[size] = cum
                    size_times_dic[size] = times
        return size_cum_dic, size_times_dic

    def construct_para_nr_run_dic(self):
        # sooo for every size we need to determine how many jobs we need to submit
        for i, size in enumerate(self.sizes):
            nr_subsystems = int(self.nr_sites / (size ** 2 * self.Ly_Lx))   # nr of subsystems per job
            min_jobs = int(
                self.min_nr_sites / self.nr_sites)  # the minimum number of jobs is the minimum_nr of total sites divided by the number of sites per job
            if nr_subsystems >= self.min_nr_systems / min_jobs:
                # if the number of systems per job * the minimumb number of jobs given by the minimum number of system sites
                # divided by the nr of sites per job is larger than the required amount of systems, wo only do the
                # minimum number of jobs
                nr_jobs = min_jobs
            else:
                # else we do as many jobs as we need to exceed the min_nr_systems
                nr_jobs = int(np.ceil(self.min_nr_systems / nr_subsystems))
            self.total_runs += nr_jobs      # we also need to know how many jobs
            # okay so nr of jobs is how many jobs we need to submit to reach the minimum requirements for the current size
            # so the keys of the para_nr_run_dic for the next nr_jobs entries shall be the parameter number corresponding to this size
            # we need to know how many keys are alread occupied
            nr_previous_jobs = len(self.para_nr_run_dic)
            for j in range(nr_jobs):
                # j ranges from 0, .., nr_jobs - 1
                self.para_nr_run_dic[nr_previous_jobs + j] = i + self.para_nr
        # this should be it again
    def get_last_time(self):
        temp_folder_path = f"{self.simulation_path}/{self.test_size:.0f}/{self.Tc:6f}"
        for file in os.listdir(temp_folder_path):
            file_path = os.path.join(temp_folder_path, file)
            if file_path.endswith(".cum"):
                # then we have the correct file
                df = pd.read_csv(file_path)
                # how do we get the last line here?
                return df.iloc[-1][0]  # or something like this?

    def run(self):
        self.setup()
        self.iteration()

    def write_test_paras(self):
        print("Writing the parameter files for testfile...")
        # We now need to construct the parameterfile with the appropriate temperature
        # Is it okay if we construct all files in the beginning and deal with the threading of the gpus later?
        # we also need to know how many subsystems we can initialize
        sys_per_job = int(np.ceil(self.nr_sites / (self.test_size ** 2 * self.Ly_Lx)))
        with open(self.filepath + "/parameters/para_set_" + str(self.get_write_para_nr()) + '.txt', 'w') as f:
                f.write(self.simulation_path)
                f.write(f"\nend_time, {1e7} \n"
                        f"dt, {self.dt} \n"
                        f"J, {self.J_para} \n"
                        f"Jy, {self.J_perp} \n"
                        f"alpha, {self.h} \n"
                        f"eta, {self.eta} \n"
                        f"nr_saves, 2 \n"           # We dont know how long the simulation will go so we could either use a density observer or weeeee just dont care
                        f"nr_repeat, 0 \n"
                        f"min_temp, {self.Tc} \n"
                        f"max_temp, {self.Tc} \n"
                        f"nr_runs, 0.0 \n"
                        f"random_init, 1.0 \n"      # for the amplitude we want to initialize randomly
                        f"curand_random, 1 \n"
                        f"subsystem_min_Lx, {self.test_size} \n"
                        f"subsystem_max_Lx, {self.test_size} \n"
                        f"nr_subsystem_sizes, 0  \n"
                        f"nr_subsystems, {sys_per_job} \n"    # The number of subsystems will be one, we use large systems that will run long to eliminate the statistical deivations
                        f"x_y_factor, {self.Ly_Lx} \n"
                        f"nr_corr_values, 0 \n"     # We need a new corr observer that just observes with density and doesnt switch after quench     
                        f"nr_ft_values, 0 \n"       # Ah we still wanted to check whether the values of the ft and fit and python or direct fit in c++ are the same, but they should be fairly similar
                        f"equil_error, {self.equil_error}\n"
                        f"min_cum_nr, {self.test_min_cum_nr}\n"
                        f"cum_write_density, {self.test_cum_density}\n"
                        f"equil_cutoff, {self.test_equil_cutoff}")
        # we need to copy the files to hemera
        rsync_command = ["rsync", "-auv", "--rsh", "ssh",
                         f"{self.filepath}/parameters/",
                         "hemera:~/Code/Master-Arbeit/CudaProject/parameters/"]
        subprocess.run(rsync_command, cwd=pathlib.Path.home())
    def write_para_files(self):
        # you ..., you know that you have to construct the parameter file at hemera?
        # and you need to do rsync after the jobs are finished!
        print("Writing the parameter files for the actual sim...")
        for i, size in enumerate(self.sizes):
            # We now need to construct the parameterfile with the appropriate temperature and size
            # Is it okay if we construct all files in the beginning and deal with the threading of the gpus later?
            # to construct the para set we need to know how many subsystems we should initialize
            nr_subsystems = int(self.nr_sites / (size ** 2 * self.Ly_Lx))
            # the endtime can now be guessed using the test equilibration time and the z guess
            endtime = (size / self.test_size) ** self.z_guess * self.test_equil_time
            # depending on the endtime we also need a density of cumulant values.
            nr_cum_values = endtime / self.dt * self.cum_density            # endtime / dt is the number of steps, times the density is the number of cum values
            with open(self.filepath + "/parameters/para_set_" + f"{self.get_write_para_nr()}" + '.txt', 'w') as f:
                f.write(self.simulation_path)
                f.write(f"\nend_time, {endtime} \n"
                        f"dt, {self.dt} \n"
                        f"J, {self.J_para} \n"
                        f"Jy, {self.J_perp} \n"
                        f"alpha, {self.h} \n"
                        f"eta, {self.eta} \n"
                        f"nr_saves, 2 \n"       # standard number of saves is 2... will this be changed anytime?
                        f"nr_repeat, 0 \n"
                        f"min_temp, {self.Tc} \n"
                        f"max_temp, {self.Tc} \n"
                        f"nr_runs, 0.0 \n"
                        f"random_init, 1.0 \n"
                        f"curand_random, 1 \n"
                        f"subsystem_min_Lx, {size} \n"
                        f"subsystem_max_Lx, {size} \n"
                        f"nr_subsystem_sizes, 0  \n"
                        f"nr_subsystems, {nr_subsystems} \n"
                        f"x_y_factor, {self.Ly_Lx} \n"
                        f"nr_cum_values, {nr_cum_values} \n"
                        f"nr_corr_values, 0 \n"
                        f"nr_ft_values, 0 \n"
                        f"equil_error, {self.equil_error}")
        # we need to copy the files to hemera
        rsync_command = ["rsync", "-auv", "--rsh", "ssh",
                         f"{self.filepath}/parameters/",
                         "hemera:~/Code/Master-Arbeit/CudaProject/parameters/"]
        subprocess.run(rsync_command, cwd=pathlib.Path.home())

    def get_write_para_nr(self):
        # this is the simple tactic that just advances the parameter number everytime it is called
        self.cur_para_nr += 1
        return self.para_nr + self.cur_para_nr - 1

    def get_para_nr(self):
        # this logic has to know how often each parameter file has to run and which parameter number to return
        # I guess the easiest way would be to construct a map or something like this that maps the cur_para_nr to the
        self.cur_run_nr += 1
        # we just increment the cur run nur and then use the dictionary that we constructed somewhere else
        return self.para_nr_run_dic[self.cur_run_nr - 1]
def main():
    # okay what is the first thing we need to do?
    # we need parameters like the number of gpus we are able to use
    nr_gpus = 15
    # we somehow need the relevant parameters
    # The model defining parameters are J_perp J_para h eta
    # the simulation defining parameters are dt
    #J_para = -150000
    J_para = -3.11
    #J_perp = -2340
    J_perp = -0.1
    #h = 1.7e6
    h = 0.5
    eta = 1.5
    #dt = 0.00001
    dt = 0.01
    max_size_Tc = 80
    min_size_Tc = 48
    nr_sizes_Tc = 3
    filepath = "/home/andi/Studium/Code/Master-Arbeit/CudaProject"
    #filepath = "/home/weitze73/Documents/Master-Arbeit/Code/Master-Arbeit/CudaProject"
    simulation_path = "../../Generated content/Silicon/Subsystems/Suite/Test8/"

    Tc_exec_file = "AutoCumulant.cu"
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

    max_rel_intersection_error = 0.02

    # Quench parameters
    max_size = 2048
    min_nr_sites = 1e6


    # Amplitude parameters
    amplitude_size = 2048
    equil_error = 0.03
    equil_cutoff = 0.1

    # z parameters
    size_min = 64
    size_max = 256
    nr_sizes = 3
    z_min_nr_sites = 1e6
    z_min_nr_systems = 500
    z_equil_error = 0.005

    # Enter which calculations are supposed to run here
    measurements = {
        "Tc": False,
        "Quench": False,
        "Amplitude": True,
        "z": False,
    }

    # I honestly have no idea on how to account h, that is really a problem
    # the Scanned interval
    if measurements["Tc"]:
        sim = crit_temp_measurement(J_para, J_perp, h, eta, dt, filepath, simulation_path + "Tc", Tc_exec_file, nr_GPUS=nr_gpus,
                                    size_min=min_size_Tc, size_max=max_size_Tc, nr_sizes=nr_sizes_Tc,
                                    intersection_error=max_rel_intersection_error)
        T_c, T_c_error = sim.routine()
    else:
        T_c = float(input("Enter critical temperature:"))
        T_c_error = 0
    if measurements["Quench"]:
        quench = quench_measurement(J_para, J_perp, h, eta, dt, filepath, simulation_path + "Quench", quench_exec_file, T_c, nr_GPUS=nr_gpus, size_max=max_size, min_nr_sites=min_nr_sites )
        quench.run()
    if measurements["Amplitude"]:
        ampl = amplitude_measurement(J_para, J_perp, h, eta, dt, filepath, simulation_path + "Amplitude",
                                     amplitude_exec_file, T_c, nr_GPUS=nr_gpus, size=amplitude_size,
                                     equil_error=equil_error, equil_cutoff=equil_cutoff)
        ampl.run()
    if measurements["z"]:
        z_measure = z_measurement(J_para, J_perp, h, eta, dt, filepath, simulation_path + "z",
                                        z_exec_file, z_test_exec_file, T_c, nr_GPUS=nr_gpus, size_min=64, size_max=256,
                                        nr_sizes=nr_sizes, min_nr_sites=z_min_nr_sites, min_nr_systems=z_min_nr_systems,
                                     equil_error=z_equil_error)
        z_measure.run()

if __name__ == '__main__':
    main()