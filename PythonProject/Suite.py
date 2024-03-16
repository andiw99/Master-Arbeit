import os

import re
from fabric import Connection
import time
from FunctionsAndClasses import *
import subprocess
import pathlib
from scipy.interpolate import CubicSpline, PchipInterpolator
from glob import glob
import statsmodels.tsa.stattools as st

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
    found_simulations = []
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

def get_avail_simulations(sizes, temperatures, directory_path, check_function, check_function_args, temp_tolerance=0.001):
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
            available_temps = [float(temp) for temp in os.listdir(size_path) if (temp[0]!="." and temp!="plots")]
            for temp in temperatures:
                for avail_temp in available_temps:
                    if abs(avail_temp - float(temp)) < temp_tolerance * float(temp):
                        use_temp = avail_temp
                        break
                    use_temp = float(temp)
                temp_path = os.path.join(size_path, f"{use_temp:.6f}")
                # Check if the temperature folder exists and is valid based on the custom function
                if (os.path.exists(temp_path) and os.path.isdir(temp_path)):
                        sim_valid = check_function(temp_path, *check_function_args)
                        if sim_valid:
                            valid_folders.append((size, temp))

    return valid_folders

def get_avail_jobs(jobs_to_do, directory_path, check_function, check_function_args):
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
    for size, temp in jobs_to_do:
        size_path = os.path.join(directory_path, str(size))

        # Check if the size folder exists
        if os.path.exists(size_path) and os.path.isdir(size_path):
            temp_path = os.path.join(size_path, f"{temp:.6f}")

            # Check if the temperature folder exists and is valid based on the custom function
            if (os.path.exists(temp_path) and os.path.isdir(temp_path)):
                    sim_valid = check_function(temp_path, *check_function_args)
                    if sim_valid:
                        valid_folders.append((size, temp))

    return valid_folders

def check_corr_valid(folderpath, equil_error, equil_cutoff, max_moving_factor=0.01):
    # This function is supposed to check whether the .corr file in the folder (There should only be one)
    # has a low enough error, the error is specified by equil error in the class
    # We also need the threshold to calculate the accurate values with the accurate errors
    # If we want to calculate the autocorrletioin time here we
    # would either have to rewrite process temp folder, which wouldnt even
    # be that dumb since this error calculation is not up to date anyway?
    # How would we actually calculate the error?
    # Weighted average other the file averages. The file averages
    # have to calculate their own autocorrlation time etc?
    # yeah probably since we cannot be sure that we can average the cumulant
    # values beforehand as we dont know if we have the same number
    # We want to change this to use a process function but we also need the moving factor
    # If we use only process file we ignore other files in the folder and
    # also we have to construct the filename from the foldername
    xix_avg, xix_error, moving_factors_x, nr_values = process_temp_folder(folderpath,
                                equil_cutoff, value="xix", file_ending="corr")
    xiy_avg, xiy_error, moving_factors_y, nr_values = process_temp_folder(folderpath,
                                equil_cutoff, value="xiy", file_ending="corr")

    # We have the error, so we can check if it is small enough
    # How do we check in the observer again? is not written actually...
    # We need the relative errors here
    moving_factor = np.mean(moving_factors_x + moving_factors_y)
    xix_rel_error = xix_error / xix_avg
    xiy_rel_error = xiy_error / xiy_avg
    avg_error = 1/ 2 * (xix_rel_error + xiy_rel_error)

    if (avg_error < equil_error) and (moving_factor < max_moving_factor):
        return True
    else:
        if (avg_error > equil_error):
            print(f"Avg error {avg_error} too large for {equil_error}")
        elif (moving_factor > max_moving_factor):
            print(f"Moving factor {moving_factor} too large")
        return False

def check_cum_valid(folderpath, equil_error, equil_cutoff, max_moving_factor, value_name="U_L", file_ending="cum",
                    process_file_func=process_file, min_cum_nr=25):
    """
    function that checks whether a cumulant measurement was valid
    :param folderpath:
    :param equil_error:
    :param equil_cutoff:
    :return:
    """
    # I think this is outdated
    # cum_avg, error = process_temp_folder(folderpath, equil_cutoff, "cum", "U_L")
    # cum, times = get_folder_average(folderpath)
    # # we have this cutoff...
    # cutoff_ind = int(equil_cutoff * len(times))
    # cum = cum[cutoff_ind:]
    # times = times[cutoff_ind:]
    # cum_avg = np.mean(cum)
    # # The folder average should in this usecase only average one file so
    # # if we somehow can extract the autocorrelation time, you know this would
    # # be great
    # txt_file = find_first_txt_file(folderpath)
    # parameters = read_parameters_txt(txt_file)
#
    # ds = times[1] - times[0]
    # try:
    #     autocorr_time = parameters["autocorrelation_time"]
    # except KeyError:
    #     autocorr_time = integrated_autocorr_time(cum, ds)
#
    # N = len(times)
    # cum_dist_var = np.var(cum)
    # variance = autocorr_time / (N * ds) * cum_dist_var
    # error = np.sqrt(variance)
#

    # Now for the moving factor

    cum_avg, rel_error, moving_factors, nr_values = process_temp_folder(folderpath,
                                equil_cutoff, value=value_name, file_ending=file_ending, process_file_func=process_file_func)

    # We have the error, so we can check if it is small enough
    # How do we check in the observer again? is not written actually...
    # We need the relative errors here
    moving_factor = np.mean(moving_factors)
    nr_vals = np.mean(nr_values)

    print(f"rel_error = {rel_error}  vs   equil_error = {equil_error}")
    if (rel_error <= equil_error) and (moving_factor <= max_moving_factor):
        return True
    else:
        if (rel_error <= equil_error) and (moving_factor > max_moving_factor):
            print(f"Moving Factor {moving_factor:.5f} is too large ! Otherwise {folderpath} would have been valid")
        elif (rel_error > equil_error) and (moving_factor <= max_moving_factor):
            print(f"Equil_error {rel_error:.5f} is too large for {equil_error}! Otherwise {folderpath} would have been valid")
        if nr_vals < min_cum_nr:
            print("Not enough values, but I hope we will just use that one as base?")
        return False

def check_cum_variation_valid(folderpath, variation_error_rate, num_per_var_val, ds):
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

def check_quench_valid(folderpath, min_nr_systems, min_nr_sites):
    """
    just checks if the simulation has more than the minimum nr of sites and systems
    :param folderpath:
    :param min_nr_systems:
    :param min_nr_sites:
    :return:
    """
    files = os.listdir(folderpath)
    total_nr_sites = 0
    total_nr_systems = 0
    for file in files:
        if file.endswith(".csv"):
            para_file = folderpath + f"/{pathlib.Path(file).stem}.txt"
            paras = read_parameters_txt(para_file)
            nr_sites = paras["total_size"]
            nr_systems = paras["nr_subsystems"]
            total_nr_sites += nr_sites
            total_nr_systems += nr_systems
    #print(folderpath)
    #print(total_nr_sites, total_nr_systems)

    if (total_nr_sites >= min_nr_sites) and (total_nr_systems >= min_nr_systems):
        return True
    else:
        return False

def check_exists(folderpath, file_ending=".cum"):
    """
    trivial check function that just returns true if the folder already exists
    :param foderpath:
    :return: true if at least one file with the required ending exists
    """
    files = [f for f in os.listdir(folderpath) if f.endswith(file_ending)]
    return bool(files)

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
    def __init__(self, J_para, J_perp, h, eta, p, dt, filepath, simulation_path, exec_file, nr_GPUS=6, Ly_Lx = 1/8,
                 host="hemera", user="weitze73", runfile="run_cuda.sh", para_nr=100, wait=60):
        # This class is supposed to encapsulate some of the functionality that the following classes share
        # For every simulation I need the basic simulation parameters
        self.J_para = J_para
        self.J_perp = J_perp
        self.h = h
        self.eta = eta
        self.p = p
        self.dt = dt
        self.filepath = filepath
        self.simulation_path = simulation_path
        self.nr_GPUS = nr_GPUS
        self.Ly_Lx = Ly_Lx

        # also some parameters for the cluster
        self.host = host                            # adress of the cluster
        self.user = user                            # user on the cluster
        self.walltime = "16:00:00"
        self.file = exec_file
        self.folder = "simulations"
        self.wait = 30

        # Besides the external simulation parameters that I have to provide there are some other attributes that
        # every autonomous suite needs
        # The bookkeeping variables for commiting the jobs
        self.running_jobs = set()
        self.completed_jobs = set()
        self.connection = None
        self.total_runs = 0
        self.para_nr = para_nr                            # every mearsurement has at least one current parameter number
        self.runfile = runfile

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
                    self.call_rsync()

                else:
                    # if it is not completed but not in jobs_on_hemera anymore,
                    # we have a problem
                    print(
                        "Oups! The Job vanished. Please check what happend")
        # remove the completed jobs from the running jobs list
        for job_id in just_completed_jobs:
            self.running_jobs.remove(job_id)

    def call_rsync(self):
        running = subprocess.call("./rsync.sh", cwd=pathlib.Path.home())
        if running:
            time.sleep(20)
            return self.call_rsync()
        else:
            return

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
                         f' logs/errors/%j.err {self.runfile} {file} {self.folder} {para_nr}'
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
    def __init__(self, J_para, J_perp, h, eta, p, dt, filepath, simulation_path, exec_file, runfile="run_cuda.sh", nr_GPUS=6, nr_Ts=5, size_min=48,
                          size_max=80, nr_sizes=3, max_steps=1e9, nr_sites=5e5, Ly_Lx = 1/8, equil_error=0.004, min_equil_error=0.001,
                 intersection_error=0.02, equil_cutoff=0.1, T_min=None, T_max=None, para_nr=100, max_moving_factor=0.005,
                 min_val_nr=2500, value_name="U_L",  file_ending="cum",
                 process_file_func=process_file, random_init=0, second=False, observed_direction=0, value_write_density=1/100):
        # call the constructor of the parent classe
        super().__init__(J_para, J_perp, h, eta, p, dt, filepath, simulation_path, exec_file,  nr_GPUS=nr_GPUS,
                         Ly_Lx=Ly_Lx, para_nr=para_nr, runfile=runfile)
        self.nr_Ts = nr_Ts
        self.size_min = size_min
        self.size_max = size_max
        self.nr_sizes = nr_sizes
        self.max_steps = max_steps
        self.nr_sites = nr_sites
        self.max_moving_factor = max_moving_factor

        # We can also specify the temperature interval ourselves If we want that
        self.T_min = T_min
        self.T_max = T_max

        self.T_arr = np.array([])
        self.max_T_step = 0.1               # fraction of the critical temperature that is the maximum stepsize for accepted measurements
        self.sizes = np.array([])
        self.random_init = random_init
        self.max_time = 0
        self.total_runs = 0

        self.all_T_arr = np.array([])       # Bookkeeping for all the temperatures we have simulated in this setting
        self.all_T_dic = {}                 # Bookkeeping dictionary for the different simulation spheres
        self.max_rel_intersec_error = intersection_error  # standard maximum error of 2%
        self.equil_error = equil_error           # standard equilibration error for the U_L runs
        self.min_equil_error = min_equil_error      # If we allow arbitrarily small errors, the runs wont end
        self.maximum_iterations = 4
        self.iteration_nr = 0
        self.repeat = False             # variable that is set to true if we have to repeat a simulation
        self.min_val_nr = min_val_nr
        self.val_write_density = value_write_density
        self.equil_cutoff = equil_cutoff
        self.jobs_to_do = None          # Bookkeeping parameter to know whether I already did some jobs / specific size / temp pairs or not

        self.cur_para_nr = 0
        self.check_function = check_cum_valid
        self.value_name = value_name
        self.second = second
        self.file_ending = file_ending
        self.process_file_func = process_file_func
        self.observed_direction = observed_direction
    def init(self):
        # this somehow needs the parameters, where do we put them? In a file? On the moon? User input?
        T_min = T_c_est(np.abs(self.J_para), np.abs(self.J_perp), self.h)[0]
        # TODO really crude approximation of the maximum T, I really should think about something better
        print(f"T_min = {T_min}\n"
              f"J_perp = {self.J_perp}\n"
              f"h={self.h}")
        # TODO this is still all a bit fishy but...
        T_max = T_min + self.h
        T_max = 2 * T_min
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
        print(f"Initializing Simulation with T_min = {self.T_min} and T_max = {self.T_max}")
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
        T_arr = np.array(self.all_T_dic[np.max(list(self.all_T_dic.keys()))])

        results = self.construct_results(self.simulation_path, self.equil_error, T_arr, self.sizes, value_name=self.value_name, file_ending=self.file_ending,
                                         process_file_func=self.process_file_func)

        # interpolate and minimize is deprecated, we use the technique we also use in iteration
        intersections, intersections_y = get_first_intersections(results, self.value_name)

        T_c = np.mean(intersections)
        U_L_intersection = np.mean(intersections_y)
        T_c_error = np.ptp(intersections)

        self.plot_value_curve(self.simulation_path, results, crit_point=(T_c, U_L_intersection), equil_error=self.equil_error)
        plt.show()

        # constructing cum dic
        cum_dic = {}
        for size in results:
            cum_dic[size] = results[size][self.value_name]

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
        # simulation_available = check_directory_structure(self.sizes, self.T_arr, self.simulation_path)
        # We check how many fitting temp-size pairs are already available
        valid_simulations = get_avail_simulations(self.sizes,
                                                  self.T_arr,
                                                  self.simulation_path,
                                                  self.check_function, check_function_args=(self.equil_error, self.equil_cutoff,
                                                                                            self.max_moving_factor,
                                                                                            self.value_name,
                                                                                            self.file_ending,
                                                                                            self.process_file_func)
                                                  )
        # Okay how do we now keep track of the sizes and temps we do not need to run anymore?
        # In the other two cases it was easy because we only had one size or one temperature always?
        # I mean we could construct a parameter with 'jobs to run' and remove
        # valid jobs from this one, but that would be one more parameter to keep
        # track of
        print("valid simulations:", valid_simulations)
        self.jobs_to_do = list(product(self.sizes, self.T_arr))

        for valid_sim in valid_simulations:
            self.jobs_to_do.remove(valid_sim)       # should work? Or should we check whether the temperatue is close enough to the requested value? I think this is already done in the get avail sim function

        if self.jobs_to_do or self.repeat:
            self.repeat = False                 # we set repeat to false again
            self.cur_para_nr = 0                # reset the parameter number
            self.write_para_files()  # setting up the parameter files for every simulation
            self.cur_para_nr = 0                # reset the parameter number
            exit()
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
        # indeed we also only want the expected sizes. The sizes luckily dont change in one simulation
        print(self.process_file_func)
        print(T_arr)
        results = self.construct_results(self.simulation_path, self.equil_error, T_arr, self.sizes, value_name=self.value_name,
                                         file_ending=self.file_ending, process_file_func=self.process_file_func)
        print(results)
        val_min_T_min_size = results[np.min(self.sizes)][self.value_name][0]
        val_min_T_max_size = results[np.max(self.sizes)][self.value_name][0]
        val_max_T_min_size = results[np.min(self.sizes)][self.value_name][-1]
        val_max_T_max_size = results[np.max(self.sizes)][self.value_name][-1]
        # we say we have an intersection if U_L_min_T_min > U_L_max_T_min
        # and U_L_min_T_max < U_L_max_T_max
        # This is we have an intersection at all, but we dont know if we have an intersection in the current T_simulation.
        intersection = ((val_min_T_min_size >= val_min_T_max_size) & (
                val_max_T_min_size <= val_max_T_max_size)) or self.intersection_anyway(val_min_T_min_size, val_min_T_max_size, val_max_T_min_size, val_max_T_max_size)       # TODO this is a bit fishy but will probably work in 99% of times, catches the case that the maximum temperature is far in the high temperaturef region and therefore inflicated with strong fluctiations
        dT = T_arr[1] - T_arr[0]
        if intersection:
            # now the usual stuff, estimate the postion of the intersection
            intersections, intersections_y = get_first_intersections(results, self.value_name)

            T_c = np.mean(intersections)

            print(f"Found an intersection at T_c / J = {T_c / self.J_para}")

            print("intersections: ", intersections)

            T_c_error = np.ptp(intersections)
            # the error will now be the minimum of this and a fourth of the stepsize
            T_c_error = max(T_c_error, dT / 2)
            print(f"T_c_error = {T_c_error}")
            rel_intersec_error = T_c_error / T_c
            print(f"rel_intersec_error = {rel_intersec_error}")
            if rel_intersec_error < self.max_rel_intersec_error:
                # In this case we are done since we covered the case of large stepsizes, small errors with the dT / 4
                print(f"Determined crit. Temp T_c = {T_c} +- {rel_intersec_error}")
                return T_c, T_c_error
            else:
                # We still want to plot the stuff to see where we at
                val_intersection = np.mean(intersections_y)
                self.plot_value_curve(self.simulation_path, results, crit_point=(T_c, val_intersection), equil_error=self.equil_error)
                plt.show()
                # If the error is too large we do a child measurement with smaller error
                self.equil_error = max(self.equil_error/ 2, self.min_equil_error)
                # We wanted to have our edgecases of within 5 or 20 percent of the interval edges...
                T_interval_low = np.max(T_arr[T_arr <= T_c])
                T_interval_up = np.min(T_arr[T_arr >= T_c]) # this should hopefully guaranteed to be always dT larger than lower interval?
                # This doesnt make sense as we saw
                #if T_c < (T_interval_low + 0.05 * dT):
                #    # in this case we want to half the new interval
                #    T_interval_up = T_interval_low + dT / 2
                #    T_interval_low -= 0.02 * dT
                #elif T_c > (T_interval_up - 0.2 * dT):
                #    T_interval_low = T_interval_up - dT / 2
                #    T_interval_up += 0.02 * dT
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
            self.plot_value_curve(self.simulation_path, results)
            plt.show()
            if val_max_T_min_size > val_max_T_max_size:
                print("The maximum temperature is too low")
                # this means we are below the critical temperature.
                # TODO can we somehow approximate how far away we are?
                # for now we just double the range i would say?
                T_range = np.ptp(self.T_arr)
                T_min = np.max(self.T_arr) + (dT)  # the next T_min is the current maximum plus the current stepsize
                T_max = np.max(self.T_arr) + T_range + dT    # TODO add dT!!
                self.T_arr = np.linspace(T_min, T_max, self.nr_Ts)
                # Okay we updated the temperature array... now we just run everything again?
            elif val_min_T_min_size < val_min_T_max_size:
                print("The minimum temperature is too high")
                T_range = np.ptp(self.all_T_arr)
                T_min = np.maximum(np.min(self.T_arr) - T_range, 0.0)  # We should not consider negative temperatures
                T_max = np.min(self.T_arr) - (self.T_arr[1] - self.T_arr[0])
                self.T_arr = np.linspace(T_min, T_max, self.nr_Ts)
            # we add the new stuff to the dictionary, but we already overwrote the T_arr, but the T_arr should still be
            # fine
            self.all_T_dic[sim_hierarchy_nr] = np.concatenate((T_arr, self.T_arr))

            return self.iteration()

    def intersection_anyway(self, val_min_T_min, val_max_T_min, val_min_T_max, val_max_T_max):
        return (val_max_T_max / val_max_T_min > 2.9)

    def get_para_nr(self):
        # this tactic inreases the parameter number everytime get_para_nr is called so that we do not submit any job twice
        self.cur_para_nr += 1
        return self.para_nr + self.cur_para_nr - 1

    def write_para_files(self):
        # you ..., you know that you have to construct the parameter file at hemera?
        # and you need to do rsync after the jobs are finished!
        print("Writing the parameter files...")
        self.total_runs = len(self.jobs_to_do)
        for i, (size, T) in enumerate(self.jobs_to_do):
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
                        f"p, {self.p} \n"
                        f"nr_saves, 4 \n"
                        f"nr_repeat, 0 \n"
                        f"min_temp, {T} \n"
                        f"max_temp, {T} \n"
                        f"nr_runs, 0.0 \n"
                        f"random_init, {self.random_init} \n"
                        f"curand_random, 1 \n"
                        f"subsystem_min_Lx, {size} \n"
                        f"subsystem_max_Lx, {size} \n"
                        f"nr_subsystem_sizes, 0  \n"
                        f"nr_subsystems, {nr_subsystems} \n"
                        f"x_y_factor, {self.Ly_Lx} \n"
                        f"nr_corr_values, 0 \n"
                        f"nr_ft_values, 0 \n"
                        f"equil_error, {self.equil_error}\n"
                        f"equil_cutoff, {self.equil_cutoff}\n"
                        f"{self.file_ending}_write_density, {self.val_write_density}\n"
                        f"min_{self.file_ending}_nr, {self.min_val_nr}\n"
                        f"moving_factor, {self.max_moving_factor}\n"
                        f"corr_second, {int(self.second)} \n"
                        f"observed_direction, {self.observed_direction}")
        # we need to copy the files to hemera
        rsync_command = ["rsync", "-auv", "--rsh", "ssh",
                         f"{self.filepath}/parameters/",
                         "hemera:~/Code/Master-Arbeit/CudaProject/parameters/"]
        subprocess.run(rsync_command, cwd=pathlib.Path.home())
    @staticmethod
    def construct_results(simulation_path, threshold, selected_temps=None, selected_sizes=None, value_name="U_L", file_ending="cum",
                          process_file_func=process_file):
        results = {}
        for size_folder in os.listdir(simulation_path):
            if size_folder[0] != 0 and size_folder != "plots":
                size_folder_path = os.path.join(simulation_path,
                                                size_folder)
                if os.path.isdir(size_folder_path):
                    if selected_sizes is not None:
                        size = int(size_folder)
                        if size not in selected_sizes:
                            continue
                    size_result = process_size_folder(size_folder_path,
                                                          threshold, selected_temperatures=selected_temps,
                                                          value=value_name, file_ending=file_ending,
                                                      process_file_func=process_file_func)
                    results[int(size_folder)] = size_result
        return results
    @staticmethod
    def plot_value_curve(simulation_path, results, crit_point=None, value_name="U_L", title="Binder Cumulant on T",
                         plotname="cum_time_avg", equil_error=None, config=None):
        fig, ax = plt.subplots(1, 1)
        y_upper_lim = 0
        y_lower_lim = np.infty
        shown_inds = np.linspace(0, len(results), len(results) + 1, endpoint=True,
                                 dtype=np.int64)
        ind = 0
        # The T array should be deduceable from the results that I am feeding
        # Every size should have the same Ts
        # results is dict in dict
        # Before we plot anything we want to know how large J is
        parapath = find_first_txt_file(simulation_path)
        parameters = read_parameters_txt(parapath)
        J_para = np.abs(parameters["J"])
        sizes = list(results.keys())
        size_min = np.min(sizes)
        size_max = np.max(sizes)
        first_dict = next(iter(results.values()))
        T_arr = first_dict["T"]
        max_T = np.max(T_arr) * 1.01 / J_para
        min_T = np.min(T_arr) * 0.99 / J_para

        for i, size in enumerate(sorted(results.keys())):
            if i in shown_inds:
                T = np.array(results[size]["T"]) / J_para
                val = np.array(results[size][value_name])
                ax.plot(T, val, marker="s", **get_point_kwargs_color(colors[5 * ind]), label=rf"$L_\parallel = {size}$")
                ax.plot(T, val, linestyle="-", color=colors[5 * ind], alpha=0.5)
                ind += 1
                if max_T:
                    y_upper_lim = np.maximum(
                        np.max(val[(min_T <= T) & (T <= max_T)]), y_upper_lim)
                    y_lower_lim = np.minimum(
                        np.min(val[(min_T <= T) & (T <= max_T)]), y_lower_lim)
        y_span = y_upper_lim - y_lower_lim
        ax.set_xlabel("$T~/~J_\parallel$")
        ax.set_ylabel(rf"${value_name}$")
        #ax.set_title(title)
        if min_T:
            ax.set_xlim(min_T, ax.get_xlim()[1])
            ax.set_ylim(y_lower_lim - 0.1 * y_span, y_upper_lim + 0.1 * y_span)
        if max_T:
            ax.set_xlim(ax.get_xlim()[0], max_T)
            ax.set_ylim(y_lower_lim - 0.1 * y_span, y_upper_lim + 0.1 * y_span)
        if crit_point:
            T_c = crit_point[0] / J_para
            val_intersection = crit_point[1]
            mark_point(ax, T_c, val_intersection,
                       label=rf"$T_c / J_\parallel = {T_c:.4f}$")
        configure_ax(fig, ax, config)
        create_directory_if_not_exists(f"{simulation_path}/plots/")
        fig.savefig(simulation_path + f"/plots/{plotname}-{size_min}-{size_max}-{equil_error}.png", format="png",
                    dpi=300, transparent=False)

        return fig, ax
class crit_temp_measurement_corr(crit_temp_measurement):
    def __init__(self, *args, **kwargs):
        # call the constructor of the parent classe
        super().__init__(*args, **kwargs)
        self.check_function = check_corr_valid
        self.file_ending = "corr"

    def plot_value_curve(self, simulation_path, results, crit_point=None):
        super().plot_value_curve(self.simulation_path, results, crit_point=crit_point, value_name=self.value_name, title="L/xi on T",
                                 plotname="corr_time_avg")

    def construct_results(self, threshold, selected_temps=None, selected_sizes=None):
        results = super().construct_results(selected_temps, selected_sizes, self.value_name, file_ending="corr")
        # And now we somehow need to find out L? The Ls are sadly different in the two directions
        # We make it easy...
        for size in results:
            if self.value_name == "xix":
                L = size
            else:
                L = size * self.Ly_Lx
            results[size][self.value_name] = L / results[size][self.value_name]
        return results

    def intersection_anyway(self, val_min_T_min, val_max_T_min, val_min_T_max, val_max_T_max):
        return False

class efficient_crit_temp_measurement(autonomous_measurement):
    def __init__(self, J_para, J_perp, h, eta, p, dt, filepath, simulation_path, exec_file, runfile="run_cuda.sh", nr_GPUS=6, size_min=48,
                 size_max=80, max_steps=1e9, nr_sites=5e5, Ly_Lx = 1/8, equil_error=0.01, min_equil_error=0.0025,
                 intersection_error=0.02, equil_cutoff=0.1, T_min=None, T_max=None, para_nr=100,
                 random_init=0, max_moving_factor=0.005, min_val_nr=5000, value_name="U_L", file_ending="cum",
                 process_file_func=process_file, val_write_density=1/100, second=False):
        # call the constructor of the parent classe
        super().__init__(J_para, J_perp, h, eta, p, dt, filepath, simulation_path, exec_file,  nr_GPUS=nr_GPUS,
                         Ly_Lx=Ly_Lx, para_nr=para_nr, runfile=runfile)

        self.value_name = value_name

        self.size_min = size_min
        self.size_max = size_max
        self.max_steps = max_steps
        self.nr_sites = nr_sites
        self.max_moving_factor = max_moving_factor

        # We can also specify the temperature interval ourselves If we want that
        self.T_min = T_min
        self.T_max = T_max
        self.random_init = random_init

        self.T_arr = np.array([])
        self.max_T_step = 0.1               # fraction of the critical temperature that is the maximum stepsize for accepted measurements
        self.sizes = np.array([])
        self.max_time = 0
        self.total_runs = 0

        self.all_T_arr = np.array([])       # Bookkeeping for all the temperatures we have simulated in this setting
        self.all_T_dic = {}                 # Bookkeeping dictionary for the different simulation spheres
        self.max_rel_intersec_error = intersection_error  # standard maximum error of 2%
        self.equil_error = equil_error           # standard equilibration error for the U_L runs
        self.min_equil_error = min_equil_error      # If we allow arbitrarily small errors, the runs wont end
        self.maximum_iterations = 4
        self.iteration_nr = 0
        self.repeat = False             # variable that is set to true if we have to repeat a simulation
        self.min_val_nr = min_val_nr
        self.val_write_density = val_write_density
        self.equil_cutoff = equil_cutoff
        # jobs to do will be a set so we cannot possible do one job accidentally twice
        self.jobs_to_do = set()          # Bookkeeping parameter to know whether I already did some jobs / specific size / temp pairs or not
        self.jobs_done = set()           # set with triplets of done jobs?

        self.cur_para_nr = 0
        self.check_function = check_cum_valid
        self.file_ending = file_ending
        self.process_file_func = process_file_func
        self.second = second
    def init(self):
        # this somehow needs the parameters, where do we put them? In a file? On the moon? User input?
        T_min = T_c_est(np.abs(self.J_para), np.abs(self.J_perp), self.h)[0]
        # TODO really crude approximation of the maximum T, I really should think about something better
        print(f"T_min = {T_min}\n"
              f"J_perp = {self.J_perp}\n"
              f"h={self.h}")
        # TODO this is still all a bit fishy but...
        T_max = T_min + self.h
        T_max = 2 * T_min
        # or just 10 %?


        if self.T_min is None:
            # If we do not specify the critical temperature, we use the critical temperature estimation
            self.T_min = T_min
        if self.T_max is None:
            self.T_max = T_max
            # If T_max is smaller, equal or slitghly larger than our specified T_min, what do we do then?

        # We use nr_Ts datapoints
        self.T_arr = np.array([self.T_min, self.T_max])
        # the T_array has to be added to the hierarchy
        # I would say in this case the keys should be the errors at which they equilibrated
        self.all_T_dic[self.equil_error] = self.T_arr  # This first array is on level 0

        self.sizes = np.array([self.size_min, self.size_max])
        print(f"Initializing Simulation with T_min = {T_min} and T_max = {T_max}")
        self.max_time = self.dt * self.max_steps
        self.total_runs = len(self.sizes) * len(self.T_arr)  # every temp size combo will be a seperate simulation
        self.jobs_to_do = set(product(self.sizes, self.T_arr))      # for the beginning this is still valid
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
        # In the new conclude only the simulation with the smallest error
        T_arr = self.all_T_dic[np.min(list(self.all_T_dic.keys()))]

        results = self.construct_results(self.equil_cutoff, T_arr, self.sizes,  value_name=self.value_name,
                                         file_ending=self.file_ending, process_file_func=self.process_file_func)

        # interpolate and minimize is deprecated, we use the technique we also use in iteration
        intersections, intersections_y = get_first_intersections(results, self.value_name)

        T_c = np.mean(intersections)
        val_intersection = np.mean(intersections_y)
        T_c_error = np.ptp(intersections)

        self.plotValueCurve(self.simulation_path, self.equil_error, (T_c, val_intersection), results, value_name=self.value_name,
                       title="Binder Cumulant on T", plotname=f"{self.file_ending}_time_avg")

        # constructing cum dic
        val_dic = {}
        for size in results:
            val_dic[size] = results[size][self.value_name]

        diff_arr, size_arr = calc_diff_at(T_c,
                                          list(results.values())[0]["T"],
                                          val_dic)

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
        ax.set_ylabel(r"$\frac{d" + f"{self.value_name}" + r"}{d \varepsilon}$")
        ax.legend()
        ax.set_title(
            r"$\frac{d" + f"{self.value_name}" + r"}{d \varepsilon}$ for different System sizes $L$")
        configure_ax(fig, ax)
        # save_plot(root, "/critical_exponent.pdf", format="pdf")
        fig.savefig(self.simulation_path + "/critical_exponent_time_avg.png",
                    format="png", dpi=250, transparent=False)
        plt.show()

    def plotValueCurve(self, simulation_path, equil_error, crit_point, results, value_name="U_L",
                       title="Binder Cumulant on T", plotname="cum_time_avg", config=None):
        fig, ax = plt.subplots(1, 1)
        y_upper_lim = 0
        y_lower_lim = np.infty
        shown_inds = np.linspace(0, len(results), len(results) + 1, endpoint=True,
                                 dtype=np.int64)
        ind = 0
        # The T array should be deduceable from the results that I am feeding
        # Every size should have the same Ts
        # results is dict in dict
        # Before we plot anything we want to know how large J is
        parapath = find_first_txt_file(simulation_path)
        parameters = read_parameters_txt(parapath)
        J_para = np.abs(parameters["J"])
        sizes = list(results.keys())
        size_min = np.min(sizes)
        size_max = np.max(sizes)
        first_dict = next(iter(results.values()))
        T_arr = first_dict["T"]
        max_T = np.max(T_arr) * 1.01 / J_para
        min_T = np.min(T_arr) * 0.99 / J_para

        for i, size in enumerate(sorted(results.keys())):
            if i in shown_inds:
                T = np.array(results[size]["T"]) / J_para
                val = np.array(results[size][value_name])
                ax.plot(T, val, marker="s", **get_point_kwargs_color(colors[5 * ind]), label=rf"$L_\parallel = {size}$")
                ax.plot(T, val, linestyle="-", color=colors[5 * ind], alpha=0.5)
                ind += 1
                if max_T:
                    y_upper_lim = np.maximum(
                        np.max(val[(min_T <= T) & (T <= max_T)]), y_upper_lim)
                    y_lower_lim = np.minimum(
                        np.min(val[(min_T <= T) & (T <= max_T)]), y_lower_lim)
        y_span = y_upper_lim - y_lower_lim
        ax.set_xlabel("$T~/~J_\parallel$")
        ax.set_ylabel(rf"${value_name}$")
        #ax.set_title(title)
        if min_T:
            ax.set_xlim(min_T, ax.get_xlim()[1])
            ax.set_ylim(y_lower_lim - 0.1 * y_span, y_upper_lim + 0.1 * y_span)
        if max_T:
            ax.set_xlim(ax.get_xlim()[0], max_T)
            ax.set_ylim(y_lower_lim - 0.1 * y_span, y_upper_lim + 0.1 * y_span)
        if crit_point:
            T_c = crit_point[0] / J_para
            val_intersection = crit_point[1]
            mark_point(ax, T_c, val_intersection,
                       label=rf"$T_c / J_\parallel = {T_c:.4f}$")
        configure_ax(fig, ax, config)
        create_directory_if_not_exists(f"{simulation_path}/plots/")
        fig.savefig(simulation_path + f"/plots/{plotname}-{size_min}-{size_max}-{equil_error}.png", format="png",
                    dpi=300, transparent=False)

        return fig, ax
    def iteration(self):
        self.iteration_nr += 1
        # Here I want to have something that checks whether there is already a measurement
        # We check how many fitting temp-size pairs are already available
        # also for corr we just check here that the required error is small enough. The stuff with the ''if it is to close together'' comes later
        # so we only have to adjust the check_function
        valid_simulations = get_avail_jobs(self.jobs_to_do,
                                                  self.simulation_path,
                                                  self.check_function, check_function_args=(self.equil_error,
                                                                                        self.equil_cutoff,
                                                                                        self.max_moving_factor,
                                                                                            self.value_name,
                                                                                            self.file_ending,
                                                                                            self.process_file_func)
                                                  )
        # Okay so the thing is even if we have a simulation that is not valid,
        # but there is at least already a simulation, we want to continue it instead of
        # rerunning so we somehow need to get this simulation into jobs done.
        # But I dont really want to change the behavior of get_avail_jobs for that
        # I think check directory structure is perfect for this? No but we can run
        # get_avail_jobs with a trival check function that returns true if it just finds a file
        avail_simulations = get_avail_jobs(self.jobs_to_do,
                                                  self.simulation_path,
                                                  check_exists, check_function_args=[self.file_ending]        # what do we do with the file ending here
                                                  )
        # Okay how do we now keep track of the sizes and temps we do not need to run anymore?
        # In the other two cases it was easy because we only had one size or one temperature always?
        # I mean we could construct a parameter with 'jobs to run' and remove
        # valid jobs from this one, but that would be one more parameter to keep
        # track of
        print("valid simulations:", valid_simulations)
         # TODO does this work with set?

        for valid_sim in valid_simulations:
            self.jobs_to_do.remove(valid_sim)       # should work? Or should we check whether the temperatue is close enough to the requested value? I think this is already done in the get avail sim function
        for avail_sim in avail_simulations:
            self.jobs_done.add(avail_sim)      # I think we dont even need the error here, if those done jobs are relevant again they by definition have had a larger error, otherwise they would have been picked up by valid simulations and would not be in jobs to do

        if self.jobs_to_do or self.repeat:
            self.repeat = False                 # we set repeat to false again
            self.cur_para_nr = 0                # reset the parameter number
            self.write_para_files()  # setting up the parameter files for every simulation
            self.cur_para_nr = 0                # reset the parameter number
            self.run_jobs()
            self.jobs_done.update(self.jobs_to_do) # add the jobs we just did to the done jobs
        else:
            print("Found valid simulation, evaluating")
        # after running the jobs, we need to
        # calculate the binder cumulant for every run
        return self.evaluate_simulation()

    def evaluate_simulation(self):
        # The simulation we are at right now is the one with the highest key
        smallest_error = np.min(list(self.all_T_dic.keys()))
        T_arr = np.array(sorted(list(self.all_T_dic[smallest_error])))
        # for this one we want to check if it has an intersection
        # indeed we also only want the expected sizes. The sizes luckily dont change in one simulation
        results = self.construct_results(self.equil_cutoff, T_arr, self.sizes,  value_name=self.value_name,
                                         file_ending=self.file_ending, process_file_func=self.process_file_func)     # If i just implement this for xi and for U_L which makes almost no difference, am I not almost done? it will just use L / xi instead of U
        # Ah another problem is with the error and the validation function, but besides that I think we are basically fine.
        val_min_T_min = np.min(results[np.min(self.sizes)][self.value_name])
        val_max_T_min = np.min(results[np.max(self.sizes)][self.value_name])
        val_min_T_max = np.max(results[np.min(self.sizes)][self.value_name])
        val_max_T_max = np.max(results[np.max(self.sizes)][self.value_name])
        # we say we have an intersection if U_L_min_T_min > U_L_max_T_min
        # and U_L_min_T_max < U_L_max_T_max
        # This is we have an intersection at all, but we dont know if we have an intersection in the current T_simulation.
        intersection = ((val_min_T_min > val_max_T_min) & (
                val_min_T_max < val_max_T_max)) or self.intersection_anyway(val_min_T_min, val_max_T_min, val_min_T_max, val_max_T_max)      # TODO this is a bit fishy but will probably work in 99% of times, catches the case that the maximum temperature is far in the high temperaturef region and therefore inflicated with strong fluctiations
        if intersection:
            # now the usual stuff, estimate the postion of the intersection
            # now we actually want all intersections, not only the one with the lowest index
            # the efficient thing is that we only use 2 sizes? x)
            T_range = results[self.size_min]["T"]
            val_small_size = results[self.size_min][self.value_name]
            val_large_size = results[self.size_max][self.value_name]

            # before we do anything we will check now if the values that we calculated lie in each others errorbars
            # and redo thos calculations but only if they are far enough away from the low and high temperature phase
            # how about we just advance all in the current error hierarchy if it is the case for one point...
            # but that would mean that we always advance all points...
            # how about we only advance the points that are part of the intersection...

            intersections, intersections_y = find_intersections(T_range, val_small_size, val_large_size)

            # If the length of the intersections is larger than 1 so if we have more than one intersection, we instantly have a problem
            if len(intersections) > 1:
                # so in what case are we here? We can only have more than one intersection if we are already
                # in the second iteration, otherwise we only have two points
                # This means we should probably redo all points that take part in the intersections
                # Remember that we could also have 3 or more intersections, although i find that quite unprobable with small
                # enough errors, we will for now deal with the two intersection case?
                # we can instantly reduce the error if we see more than two intersections
                self.equil_error = max(self.equil_error / 2, self.min_equil_error)
                # now we redo the points in question, the points in question would be the points at the edges
                # of the intersections and also all points between, although I dont think there would be any with this
                # kind of algorithm
                self.all_T_dic[self.equil_error] = set()
                for intersection in intersections:
                    lower_T_bound = np.max(T_arr[T_arr < intersection])
                    upper_T_bound = np.min(T_arr[T_arr > intersection])
                    for size in self.sizes:
                        self.jobs_to_do.add((size, lower_T_bound))
                        self.jobs_to_do.add((size, upper_T_bound))
                    self.all_T_dic[self.equil_error].add(lower_T_bound)
                    self.all_T_dic[self.equil_error].add(upper_T_bound)
                # then we can instantly retake the measurement?

                return self.iteration()
            else:
                # else means we have exactly 1 intersection, (or zero if the 2.95 rule catches)
                if len(intersections) == 1:
                    intersection = intersections[0]
                    lower_T_bound = np.max(T_arr[T_arr < intersection])
                    upper_T_bound = np.min(T_arr[T_arr > intersection])
                    # since T is ordered it is just always zero in this case
                    # lower_T_bound_ind = np.argmax(T_arr[T_arr < intersection])
                    # upper_T_bound_ind = np.argmin(T_arr[T_arr > intersection])

                    val_lower_bound_small_size = results[self.size_min][self.value_name][T_arr < intersection][-1]
                    val_lower_bound_large_size = results[self.size_max][self.value_name][T_arr < intersection][-1]
                    val_upper_bound_small_size = results[self.size_min][self.value_name][T_arr > intersection][0]
                    val_upper_bound_large_size = results[self.size_max][self.value_name][T_arr > intersection][0]

                    errors_overlap = self.check_errors_overlap(lower_T_bound, upper_T_bound, val_lower_bound_large_size,
                                                               val_lower_bound_small_size, val_upper_bound_large_size,
                                                               val_upper_bound_small_size)
                    if errors_overlap:
                        return self.iteration()
                    T_c = np.mean(intersections)
                    print(f"Found an intersection at T_c = {T_c}")
                    print("intersections: ", intersections)
                    T_c_error = np.ptp(intersections)
                elif not intersections:
                    # means intersection is empty so we measured far away U_L points that do not intersection
                    T_c = np.mean(T_arr)
                    T_c_error = np.infty

                T_interval_low = np.max(T_arr[T_arr <= T_c])
                T_interval_up = np.min(T_arr[T_arr >= T_c]) # this should hopefully guaranteed to be always dT larger than lower interval?
                dT = T_interval_up - T_interval_low
                T_c_error = max(T_c_error, dT / 2)
                print(f"T_c_error = {T_c_error}")
                rel_intersec_error = T_c_error / T_c
                print(f"rel_intersec_error = {rel_intersec_error}")
                if rel_intersec_error < self.max_rel_intersec_error:
                    # In this case we are done since we covered the case of large stepsizes, small errors with the dT / 4
                    print(f"Determined crit. Temp T_c = {T_c} +- {rel_intersec_error}")
                    return T_c, T_c_error
                else:
                    val_intersection = np.mean(intersections_y)
                    self.plotValueCurve(self.simulation_path, self.equil_error, (T_c, val_intersection), results, value_name=self.value_name, title="Binder Cumulant on T", plotname="cum_mag_avg")

                    T_low = T_interval_low + dT/3
                    T_up = T_interval_low + 2 * dT / 3
                    new_Ts = np.linspace(T_low, T_up, 2)
                    self.all_T_dic[self.equil_error] = set(np.sort(np.concatenate((T_arr, new_Ts))))
                    self.T_arr = new_Ts
                    # and update the jobs to do
                    self.jobs_to_do = set(product(self.sizes, self.T_arr))
                    # If we are here this directly means that we started a new child a a new level in the hierarchy

                    print(f"Error was too large: Temp T_c = {T_c} +- {T_c_error} \n"
                          f"Starting new run with T_min = {T_low}, T_max = {T_up}")
                    return self.iteration()
        else:
            # This means we missed the intersection so we want to restart a simulation with the same error
            print("We do not see a intersection")
            self.plotValueCurve(self.simulation_path, self.equil_error, None, results, self.value_name,
                                title="No Intersection", plotname=f"{self.file_ending}_time_avg")
            if val_min_T_max > val_max_T_max:
                print("The maximum temperature is too low")
                # this means we are below the critical temperature.
                # TODO can we somehow approximate how far away we are?
                # for now we just double the range i would say?
                T_range = np.ptp(self.T_arr)
                dT = T_arr[1] - T_arr[0]
                T_min = np.max(self.T_arr) + (dT)  # the next T_min is the current maximum plus the current stepsize
                T_max = np.max(self.T_arr) + T_range + dT    # TODO add dT!!
                self.T_arr = np.linspace(T_min, T_max, 2)
                # Okay we updated the temperature array... now we just run everything again?
            elif val_min_T_min < val_max_T_min:
                print("The minimum temperature is too high")
                T_range = np.ptp(self.T_arr)
                T_min = np.maximum(np.min(self.T_arr) - T_range, 0.0)  # We should not consider negative temperatures
                T_max = np.min(self.T_arr) - (self.T_arr[1] - self.T_arr[0])
                self.T_arr = np.linspace(T_min, T_max, 2)
            # we add the new stuff to the dictionary, but we already overwrote the T_arr, but the T_arr should still be
            # fine
            self.all_T_dic[self.equil_error] = set(np.sort(np.concatenate((T_arr, self.T_arr))))
            for size in self.sizes:
                self.jobs_to_do.add((size, T_min))
                self.jobs_to_do.add((size, T_max))
            return self.iteration()

    def check_errors_overlap(self, lower_T_bound, upper_T_bound, val_lower_bound_large_size, val_lower_bound_small_size,
                             val_upper_bound_large_size, val_upper_bound_small_size):
        # Those numbers are kind of willkürlihc
        errors_overlap = False
        if (self.check_if_values_to_close(val_upper_bound_small_size, val_lower_bound_large_size)
                and self.equil_error > self.min_equil_error):
            # otherwise we say that the U_L span is so large that we aller wahrscheinlichkeit nach have an intersection between
            # for the correlation length we here have the problem that the error on L/xi is not just the error on xi
            # and we (kinda) lost the information on xi
            if (val_lower_bound_small_size * (1 - self.equil_error) <
                val_lower_bound_large_size * (1 + self.equil_error)) or (
                    val_upper_bound_small_size * (1 + self.equil_error) >
                    val_upper_bound_large_size * (1 - self.equil_error)):
                # this means the error bars overlap and we redo the corresponding jobs with half the error
                self.equil_error = max(self.equil_error / 2, self.min_equil_error)
                for size in self.sizes:
                    self.jobs_to_do.add((size, lower_T_bound))
                    self.jobs_to_do.add((size, upper_T_bound))
                self.all_T_dic[self.equil_error] = {lower_T_bound, upper_T_bound}
                errors_overlap = True
        return errors_overlap

    def intersection_anyway(self, val_min_T_min, val_max_T_min, val_min_T_max, val_max_T_max):
        return (val_max_T_max / val_max_T_min > 2.6)

    def check_if_values_to_close(self, val_upper_bound_small_size, val_lower_bound_large_size):
        return (val_upper_bound_small_size < 2.95) and (val_lower_bound_large_size > 1.01)

    def construct_results(self, threshold, selected_temps=None, selected_sizes=None,
                          value_name="U_L", file_ending="cum", process_file_func=process_file):
        results = {}
        for size_folder in os.listdir(self.simulation_path):
            if size_folder != "plots":
                size_folder_path = os.path.join(self.simulation_path,
                                                size_folder)
                if os.path.isdir(size_folder_path) and size_folder[0] != ".":
                    if selected_sizes is not None:
                        size = int((size_folder))
                        if size not in selected_sizes:
                            continue
                    if (size_folder[0] != ".") & (size_folder != "plots"):
                        size_result = process_size_folder(size_folder_path,
                                                          threshold, selected_temperatures=selected_temps,
                                                          value=value_name, file_ending=file_ending, process_file_func=process_file_func)
                        results[int((size_folder))] = size_result
        return results

    def get_para_nr(self):
        # this tactic inreases the parameter number everytime get_para_nr is called so that we do not submit any job twice
        self.cur_para_nr += 1
        return self.para_nr + self.cur_para_nr - 1

    def write_advance_file(self, size, T, ind, mode=-1, done_path=None):
        # we need to find the name of the file...
        if not done_path:
            folder_path = self.simulation_path + f"/{size}/{T:.6f}"
            csv_path = find_first_csv_file(folder_path)
            name = pathlib.Path(csv_path).stem
            filepath = f"{folder_path}/{name}"
        else:
            filepath = done_path
        nr_subsystems = int(self.nr_sites / (size ** 2 * self.Ly_Lx))
        with open(self.filepath + "/parameters/para_set_" + str(self.para_nr + ind) + '.txt', 'w') as f:
            f.write(filepath)
            f.write(f"\nend_time, {self.max_time} \n"
                    f"dt, {self.dt} \n"
                    f"J, {self.J_para} \n"
                    f"Jy, {self.J_perp} \n"
                    f"alpha, {self.h} \n"
                    f"eta, {self.eta} \n"
                    f"p, {self.p} \n"
                    f"nr_saves, 4 \n"                    
                    f"nr_repeat, 0 \n"
                    f"min_temp, {T} \n"
                    f"max_temp, {T} \n"
                    f"nr_runs, 0.0 \n"
                    f"random_init, {mode} \n"     # important, this has to be -1
                    f"curand_random, 1 \n"
                    f"subsystem_min_Lx, {size} \n"
                    f"subsystem_max_Lx, {size} \n"
                    f"nr_subsystem_sizes, 0  \n"
                    f"nr_subsystems, {nr_subsystems} \n"
                    f"x_y_factor, {self.Ly_Lx} \n"
                    f"nr_corr_values, 0 \n"
                    f"nr_ft_values, 0 \n"
                    f"equil_error, {self.equil_error}\n"
                    f"equil_cutoff, {self.equil_cutoff}\n"
                    f"{self.file_ending}_write_density, {self.val_write_density}\n"
                    f"min_{self.file_ending}_nr, {self.min_val_nr}\n"
                    f"moving_factor, {self.max_moving_factor} \n"
                    f"corr_second, {int(self.second)}")

    def write_new_file(self, size, T, ind):
        # We now need to construct the parameterfile with the appropriate temperature and size
        # Is it okay if we construct all files in the beginning and deal with the threading of the gpus later?
        # to construct the para set we need to know how many subsystems we should initialize
        nr_subsystems = int(self.nr_sites / (size ** 2 * self.Ly_Lx))
        with open(self.filepath + "/parameters/para_set_" + str(self.para_nr + ind) + '.txt', 'w') as f:
            f.write(self.simulation_path)
            f.write(f"\nend_time, {self.max_time} \n"
                    f"dt, {self.dt} \n"
                    f"J, {self.J_para} \n"
                    f"Jy, {self.J_perp} \n"
                    f"alpha, {self.h} \n"
                    f"eta, {self.eta} \n"
                    f"p, {self.p} \n"
                    f"nr_saves, 4 \n"
                    f"nr_repeat, 0 \n"
                    f"min_temp, {T} \n"
                    f"max_temp, {T} \n"
                    f"nr_runs, 0.0 \n"
                    f"random_init, {self.random_init} \n"
                    f"curand_random, 1 \n"
                    f"subsystem_min_Lx, {size} \n"
                    f"subsystem_max_Lx, {size} \n"
                    f"nr_subsystem_sizes, 0  \n"
                    f"nr_subsystems, {nr_subsystems} \n"
                    f"x_y_factor, {self.Ly_Lx} \n"
                    f"nr_corr_values, 0 \n"
                    f"nr_ft_values, 0 \n"
                    f"equil_error, {self.equil_error}\n"
                    f"equil_cutoff, {self.equil_cutoff}\n"
                    f"{self.file_ending}_write_density, {self.val_write_density}\n"
                    f"min_{self.file_ending}_nr, {self.min_val_nr}\n"
                    f"moving_factor, {self.max_moving_factor}\n"
                    f"corr_second, {int(self.second)}")
    def write_para_files(self):
        # you ..., you know that you have to construct the parameter file at hemera?
        # and you need to do rsync after the jobs are finished!
        print("Writing the parameter files...")
        self.total_runs = len(self.jobs_to_do)
        for i, (size, T) in enumerate(self.jobs_to_do):
            if (size, T) in self.jobs_done:
                # we already have this job with a too large error
                self.write_advance_file(size, T, i, -1)
            elif self.jobs_done:
                # If jobs done is not empty this means we alread have a measurment
                # with the same parameters but at a different temperature
                best_done_T = 0
                min_T_diff = np.infty
                for ind, (done_size, done_temp) in enumerate(self.jobs_done):
                    # the size has to coincide
                    if done_size == size:
                        temp_diff = np.abs(T - done_temp)
                        if temp_diff < min_T_diff:
                            best_done_T = done_temp
                            min_T_diff = temp_diff
                if best_done_T:
                    done_folder_path = self.simulation_path + f"/{size}/{best_done_T:.6f}"
                    done_csv_path = find_first_csv_file(done_folder_path)
                    name = pathlib.Path(done_csv_path).stem
                    done_path = f"{done_folder_path}/{name}"

                    self.write_advance_file(size, T, i, -2, done_path)
                else:
                    self.write_new_file(size, T, i)
            else:
                self.write_new_file(size, T, i)
        # we need to copy the files to hemera
        rsync_command = ["rsync", "-auv", "--rsh", "ssh",
                         f"{self.filepath}/parameters/",
                         "hemera:~/Code/Master-Arbeit/CudaProject/parameters/"]
        subprocess.run(rsync_command, cwd=pathlib.Path.home())

class efficient_crit_temp_measurement_corr(efficient_crit_temp_measurement):
    def __init__(self, *args, **kwargs):
        # call the constructor of the parent classe
        super().__init__(*args, **kwargs)
        self.check_function = check_corr_valid
        self.file_ending = "corr"


    def plotValueCurve(self, T_c, val_intersection, results):
        super().plotValueCurve(self.simulation_path, self.equil_error, (T_c, val_intersection), results, value_name=self.value_name, title="L/xi on T",
                               plotname="corr_time_avg")

    def construct_results(self, threshold, selected_temps=None, selected_sizes=None):
        results = super().construct_results(threshold, selected_temps, selected_sizes, self.value_name,
                                            file_ending="corr")
        # And now we somehow need to find out L? The Ls are sadly different in the two directions
        # We make it easy...
        for size in results:
            if self.value_name == "xix":
                L = size
            else:
                L = size * self.Ly_Lx
            results[size][self.value_name] = L / results[size][self.value_name]
        return results

    def intersection_anyway(self, val_min_T_min, val_max_T_min, val_min_T_max, val_max_T_max):
        return False

    def check_errors_overlap(self, lower_T_bound, upper_T_bound, val_lower_bound_large_size, val_lower_bound_small_size,
                             val_upper_bound_large_size, val_upper_bound_small_size):
        # Those numbers are kind of willkürlihc
        errors_overlap = False
        if self.equil_error > self.min_equil_error:
            # otherwise we say that the U_L span is so large that we aller wahrscheinlichkeit nach have an intersection between
            # for the correlation length we here have the problem that the error on L/xi is not just the error on xi
            # and we (kinda) lost the information on xi

            low_overlap = (val_lower_bound_small_size * (1 - self.get_error(val_lower_bound_small_size, self.size_min))
                           < val_lower_bound_large_size * (1 + self.get_error(val_lower_bound_large_size, self.size_max)))
            up_overlap = (val_upper_bound_small_size * (1 + self.get_error(val_upper_bound_small_size, self.size_min)) >
                          val_upper_bound_large_size * (1 - self.get_error(val_upper_bound_large_size, self.size_max)))
            if low_overlap or up_overlap:
                # this means the error bars overlap and we redo the corresponding jobs with half the error
                self.equil_error = max(self.equil_error / 2, self.min_equil_error)
                for size in self.sizes:
                    self.jobs_to_do.add((size, lower_T_bound))
                    self.jobs_to_do.add((size, upper_T_bound))
                self.all_T_dic[self.equil_error] = {lower_T_bound, upper_T_bound}
                errors_overlap = True
        return errors_overlap

    def get_error(self, L_xi, L):
        # okay so  1 / L_xi * L is xi
        xi = 1 / L_xi * L
        abs_error = xi * self.equil_error
        abs_errro_L_xi = L_xi / xi * abs_error
        rel_error_L_xi = abs_errro_L_xi / L_xi
        return rel_error_L_xi

    def get_L(self, L):
        if self.value_name == "xix":
            return L
        else:
            return int(L * self.Ly_Lx)
class quench_measurement(autonomous_measurement):
    def __init__(self, J_para, J_perp, h, eta, p, dt, filepath, simulation_path, exec_file, runfile, Tc,
                 nr_GPUS=6, size_min=64, size_max=4096, nr_sites=5e5, Ly_Lx=1/8,
                 min_quench_steps=100, min_nr_sites=1e6, min_nr_systems=10,
                 host="hemera", user="weitze73", wait=30, max_nr_steps=1e7, para_nr=100, tau_max=np.infty):
        super().__init__(J_para, J_perp, h, eta, p, dt, filepath, simulation_path, exec_file,  nr_GPUS=nr_GPUS,
                         Ly_Lx=Ly_Lx, wait=wait, para_nr=para_nr, runfile=runfile)
        self.size_min = size_min        # The starting size at which we do go higher
        self.size_max = size_max        # maximum size, if xi = size_max / 10 we stop the simulation
        self.nr_sites = nr_sites        # nr of sites for one run
        self.Tc = Tc                    # The critical temperature that was calculated in the measurement before
        self.min_quench_steps = min_quench_steps    # the minimum number of steps done during a quench, influences tau_min
        self.min_nr_sites = min_nr_sites            # the minimum nr of sites i want to simulate for every tau
        self.min_nr_systems = min_nr_systems        # the minimum nr of systems to go through the qunech, guaranteeing an approximately precice xi value
        self.max_nr_steps = max_nr_steps            # for system whose correlation lengths rise very slowly only

        self.T_start = Tc
        self.T_end = Tc
        self.tau_min = 1
        self.tau_max = tau_max
        self.tau_factor = 1                    # current tau
        self.tau_list = []                 # empty list to keep track of the taus that we used
        self.size = size_min            # current size, starts at size_min
        self.equil_error = 0.05        # the equilibration error, so the error of U_L for which we assume that we are approximately equilibrated, doesnt need to be as small as for the T_c measurement
        self.para_nr = 100
        self.min_nr_corr_values = 500

        self.cut_zero_impuls = True     # will probably always be true since we are quenching
        self.fitfunc = lorentz_offset      # We usually use this atm?
        self.nr_measured_values = 300   # standard number of measured values during the quench
        self.min_tau_scaling_fit = 10   # this is whack, it should automatically find the best regression

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
                    f"min_cum_nr, 50\n"
                    f"min_corr_nr, {self.min_nr_corr_values}\n"
                    f"corr_second, {int(self.second)}")
        # we need to copy the files to hemera
        rsync_command = ["rsync", "-auv", "--rsh", "ssh",
                         f"{self.filepath}/parameters/",
                         "hemera:~/Code/Master-Arbeit/CudaProject/parameters/"]
        subprocess.run(rsync_command, cwd=pathlib.Path.home())

    def iteration(self):
        # Check directory structure again
        sims_available = get_avail_simulations([self.size], [self.tau()], self.simulation_path,
                                              check_quench_valid,
                                              check_function_args=(self.min_nr_systems, self.min_nr_sites))      # the check directory structure takes lists in
        # one of the first things we do in the iteration is writing the parameter file
        if not sims_available:
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
            # Calculate how many steps this would mean
            nr_steps = (self.T_start - self.T_end) * self.tau() / self.dt
            if nr_steps > self.max_nr_steps:
                print(f"tau = {self.tau()} would require more than {self.max_nr_steps}, aborting here")
                return
            elif self.tau() > self.tau_max:
                print(f"Reached maximum tau {self.tau()} > {self.tau_max}, aborting here")
                return
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
        # Now we want to do some fitting
        # I have the feeling that it will be faster if i just rewrite it
        # Okay so we will use from the size_tau_xi_dics for every size
        # all but the last value. Then we should have a continous tau-xi curve
        # About the minimum tau, I don't know how to deal with that atm,
        # we will just set it to 10 or let it be a parameter
        # We dont need a max tau in this case as the largest tau value is excluded and all other are basically guaranteed to be accurate
        # xix_scaling = np.array(xix_scaling)[
        #     np.array(tau_scaling) > self.min_tau_scaling_fit]
        # xiy_scaling = np.array(xiy_scaling)[
        #     np.array(tau_scaling) > self.min_tau_scaling_fit]
        # tau_scaling = np.array(tau_scaling)[
        #     np.array(tau_scaling) > self.min_tau_scaling_fit]
        #
        # # Do the fitting
        # popt_x, _ = curve_fit(linear_fit, np.log(tau_scaling),
        #                       np.log(xix_scaling))
        # popt_y, _ = curve_fit(linear_fit, np.log(tau_scaling),
        #                       np.log(xiy_scaling))
        # quench_exp_x = popt_x[0]
        # quench_ampl_x = popt_x[1]
        # quench_exp_y = popt_y[0]
        # quench_ampl_y = popt_y[1]
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
        size_tau_xix_dic, size_tau_xiy_dic = self.get_size_quench_results(self.simulation_path,
                                                                          self.cut_zero_impuls, self.fitfunc)

        tau_scaling, xix_scaling, reg_x, max_tau_ind_x, min_tau_ind_x = quench_measurement.fit_kzm(
            size_tau_xix_dic)

        tau_scaling, xiy_scaling, reg_y, max_tau_ind_y, min_tau_ind_y = quench_measurement.fit_kzm(
            size_tau_xiy_dic)


        quench_measurement.plot_kzm_scaling(tau_scaling, size_tau_xix_dic, reg_x, max_tau_ind_x, direction="parallel")
        plt.savefig(self.simulation_path + f"/plots/tau-xi-parallel.png", format="png")
        plt.show()

        self.plot_kzm_scaling(tau_scaling, size_tau_xiy_dic, reg_y, max_tau_ind_y, min_tau_ind_y)
        plt.savefig(self.simulation_path + f"/plots/tau-xi-perp.png", format="png")
        plt.show()

        self.plot_ratio_after_quench(tau_scaling, xix_scaling, xiy_scaling, min_tau_ind_x)
        plt.savefig(self.simulation_path + "/plots/xix_xiy.png", format="png")
        plt.show()


    @staticmethod
    def plot_quench_process(simulation_path, taus, xi_ampl, Tc, direction="parallel", cut_around_peak=True,
                            cut_zero_impuls=True, peak_cut_threshold=0.2, set_fts_to_zero=False, min_points_fraction=0.2,
                             fitfunc=lorentz_offset):
        # what is this going to do, just plotting the quench process and
        avail_taus = get_avail_tau_dic(simulation_path)

        tau_sizes = []
        tau_paths = []
        for tau in taus:
            try:
                largest_avail_size = np.max(avail_taus[tau])
                tau_sizes.append(largest_avail_size)
                tau_paths.append(f"{simulation_path}/{largest_avail_size}/{tau:.6f}")
            except KeyError:
                print(f"tau = {tau} not available")

        # now the stuff with the twofold plot that you made sometime
        fig, axes = plt.subplots(1, 2, figsize=(10, 7), gridspec_kw={'width_ratios': [1, 3]})
        ax_equil = axes[0]
        ax_quench = axes[1]

        for i, (tau, taupath) in enumerate(zip(taus, tau_paths)):
            if i == 0:
                parafile = find_first_txt_file(taupath)
                parameters = read_parameters_txt(parafile)

            t_equil, t_process, xix_equil, xiy_equil, xix_process, xiy_process = quench_measurement.get_time_resolved_xi(taupath, cut_around_peak,
                            cut_zero_impuls, peak_cut_threshold, set_fts_to_zero, min_points_fraction, fitfunc)
            if direction == "parallel":
                xi_equil = xix_equil
                xi_process = xix_process
            else:
                xi_equil = xix_equil
                xi_process = xix_process

            ax_equil.plot(t_equil, xi_equil, linestyle="", markersize=5, marker=markers[2 * i],
                               color=colors[2 * i], markerfacecolor="none")
            t_process /= tau
            ax_quench.plot(t_process, xi_process, linestyle="", markersize=5, marker=markers[2 * i],
                               color=colors[2 * i], markerfacecolor="none")

       # plot the divergence
        # get the limits beforehand
        y_limits = ax_quench.get_ylim()
        T_start = parameters["starting_temp"]
        T = T_start - t_process
        eps = (T - Tc) / Tc

        pos_eps = eps[eps >= 0]
        neg_eps = eps[eps < 0]
        U_xi = 1.95     # I think this is the Ising amplitude ratio? negative should be larger

        nu = 1
        xi_equilibrium_pos = xi_div(pos_eps, xi_ampl, nu)
        xi_equilibrium_neg = xi_div(neg_eps, xi_ampl * U_xi, nu)

        t_pos_temp = t_process[eps >= 0]
        t_neg_temp = t_process[eps < 0]

        ax_quench.plot(t_pos_temp, xi_equilibrium_pos, color=colors[0])
        ax_quench.plot(t_neg_temp, xi_equilibrium_neg, color=colors[0])

        ax_quench.set_ylim(y_limits)
        ax_equil.set_ylim(y_limits)
        ax_equil.set_xlim(ax_equil.get_xlim()[0], 0)
        ax_quench.set_xlim(0, ax_quench.get_xlim()[1])
        ax_equil.set_ylabel(rf"$\xi_\{direction} / a_\{direction}$")
        ax_equil.set_xlabel("t/ns")
        ax_quench.set_xlabel(
            r"$t / \tau_Q$")

        quench_config = {
            "titlesize": 0,
            "ytickfontsize": 0,
            "ylabelsize": 0,
            "y_tickwidth": 0,
            "y_ticklength": 0
        }
        equil_config = {
            "nr_x_major_ticks": 2
        }
        configure_ax(fig, ax_quench, quench_config)
        configure_ax(fig, ax_equil, equil_config)
        fig.subplots_adjust(wspace=0.01)


        # plt.savefig(simulation_path + f"quench-process-{direciton}.png")
        # plt.show()
        return fig, axes

    @staticmethod
    def get_time_resolved_xi(folderpath, cut_around_peak=True, cut_zero_impuls=True,
                                             peak_cut_threshold=0.2, set_fts_to_zero=False, min_points_fraction=0.2,
                             fitfunc=lorentz_offset):
        # are you going to get clapped because you again didnt divide?
        ft_k, ft_l = average_ft_unequal_times(folderpath)

        xix = []
        ts = []
        xiy = []

        parapath = find_first_txt_file(folderpath)
        parameters = read_parameters_txt(parapath)
        Lx = parameters["subsystem_Lx"]
        Ly = parameters["subsystem_Ly"]
        tau = parameters["tau"]
        T_start = parameters["starting_temp"]
        T_end = parameters["end_temp"]
        quench_time = (T_start - T_end) * tau


        for t in ft_k:
            ft_k_fit = ft_k[t]
            ft_k_fit, p_k = prepare_fit_data(cut_around_peak, cut_zero_impuls, ft_k_fit,
                                             peak_cut_threshold, set_fts_to_zero, min_points_fraction)

            popt_x, perr_x = fit_lorentz(p_k, ft_k_fit,
                                         fitfunc=fitfunc)
            # print("offset = ", popt_x[1])
            xi = np.minimum(np.abs(popt_x[0]), Lx)
            ts.append(t)
            xix.append(xi)
        for t in ft_l:
            ft_l_fit = ft_l[t]
            ft_l_fit, p_l = prepare_fit_data(cut_around_peak, cut_zero_impuls, ft_l_fit,
                                             peak_cut_threshold, set_fts_to_zero, min_points_fraction)

            popt_y, perr_y = fit_lorentz(p_l, ft_l_fit,
                                         fitfunc=fitfunc)
            xi = np.minimum(np.abs(popt_y[0]), Ly)
            xiy.append(xi)

        ts = np.array(ts)
        xix = np.array(xix)
        xiy = np.array(xiy)
        t_end = np.max(ts)

        equil_time = t_end - quench_time

        xix_equil = xix[ts <= equil_time]
        xiy_equil = xiy[ts <= equil_time]

        xix_process = xix[ts > equil_time]
        xiy_process = xiy[ts > equil_time]

        t_equil = ts[ts <= equil_time] - equil_time
        t_process = ts[ts > equil_time] - equil_time

        return t_equil, t_process, xix_equil, xiy_equil, xix_process, xiy_process

    @staticmethod
    def plot_equilibration_phase(ax, t, xi):
        ax.plot(t, xi)
    @staticmethod
    def plot_process_process(ax, t, xi):
        # what is this supposed to do, take an ax and plot the quench, just the quench and it shoud just take in
        # the correlation length and the times?
        ax.plot(t, xi)
    @staticmethod
    def fit_kzm(size_tau_xix_dic):
        tau_scaling = []
        xix_scaling = []
        for size in sorted(list(size_tau_xix_dic.keys())):
            tau_xix_dic = size_tau_xix_dic[size]
            # from those dictionaries I want to extract the largest tau. Are the taus floats here?
            for tau in np.sort(list(tau_xix_dic.keys()))[:-1]:
                # By the -1 we exclude hopefully the largest tau
                tau_scaling.append(tau)
                xix_scaling.append(tau_xix_dic[tau])
        tau_scaling = np.array(tau_scaling)
        log_tau = np.log(tau_scaling)
        xix_scaling_log = np.log(xix_scaling)
        reg_x, min_tau_ind_x, max_tau_ind_x = best_lin_reg(log_tau, xix_scaling_log,
                                                           min_r_squared=0.95, more_points=True,
                                                           min_points=3,
                                                           require_end=True)  # I mean it would be sad if we spend the most time on the last datapoint and werent to use it?
        return tau_scaling, xix_scaling, reg_x, max_tau_ind_x, min_tau_ind_x

    @staticmethod
    def plot_kzm_scaling(tau_scaling, size_tau_xi_dic, reg, max_tau_ind,
                         min_tau_ind, direction="parallel"):
        quench_exp = reg.slope
        quench_ampl = np.exp(reg.intercept)
        # xiy scaling
        figy, ax = plt.subplots(1, 1)
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_xlabel(r"$\tau$")
        ax.set_ylabel(rf"$\xi_\{direction}$")


        for i, size in enumerate(size_tau_xi_dic):
            # Plotting, every size should get its own color and/or symbol?
            # construct the lists to plot
            tau = list(size_tau_xi_dic[size].keys())
            xi = list(size_tau_xi_dic[size].values())

            ax.plot(tau, xi, marker=markers[i], linestyle="None", label=rf"$L_\{direction}$ = ${size}", color="C1")
        prev_y_low = ax.get_ylim()[0]
        prev_y_up = ax.get_ylim()[1]
        ax.plot(tau_scaling[min_tau_ind:max_tau_ind],
                 poly(tau_scaling[min_tau_ind:max_tau_ind], quench_exp, (quench_ampl)),
                 color="black", alpha=0.5, linestyle="dashed",
                 label=r"$\frac{\nu}{1 + \nu z} =$" + f"{quench_exp:.2f}")
        ax.set_ylim(prev_y_low, prev_y_up)
        configure_ax(figy, ax)


    @staticmethod
    def plot_ratio_after_quench(tau_scaling, xix_scaling, xiy_scaling, min_tau_ind_x):
        # ratio
        fig, ax = plt.subplots(1, 1)
        ax.set_yscale("log")
        ax.set_xscale("log")
        ax.set_xlabel(r"$\tau$")
        ax.set_ylabel(r"$\xi_y$")
        xix_xiy_ratio = np.array(xix_scaling) / np.array(xiy_scaling)
        ax.plot(tau_scaling, xix_xiy_ratio, color="C0", linestyle="", marker="s", markerfacecolor=None,
                markeredgecolor="C0", label=r"$\hat{\xi_x} / \hat{\xi_y}$")
        xix_xiy = np.mean(xix_xiy_ratio[min_tau_ind_x:])
        xix_xiy_fast = np.mean(xix_xiy_ratio[:min_tau_ind_x])
        ax.plot([], [], linestyle="", label=rf"$\langle \xi_x / \xi_y \rangle = {xix_xiy}$")
        ax.plot([], [], linestyle="",
                label=r"$\langle \xi_x / \xi_y (\tau < \tau_{min}) \rangle = $" + f"{xix_xiy_fast}")
        configure_ax(fig, ax)

    @staticmethod
    def get_size_quench_results(simulation_path, cut_zero_impuls, fitfunc):
        size_tau_xix_dic = {}
        size_tau_xiy_dic = {}
        for size in os.listdir(simulation_path):
            tau_xix = {}
            tau_xiy = {}
            if (size != "plots") & (size[0] != "."):
                sizepath = os.path.join(simulation_path, size)
                if os.path.isdir(sizepath):
                    for tau in os.listdir(sizepath):
                        if (tau != "plots") & (tau[0] != "."):
                            taupath = os.path.join(sizepath, tau)
                            if os.path.isdir(taupath):

                                ft_k, ft_l = average_lastline_ft(taupath)
                                p_k = get_frequencies_fftw_order(len(ft_k))
                                p_l = get_frequencies_fftw_order(len(ft_l))

                                if cut_zero_impuls:
                                    p_k, ft_k = cut_zero_imp(p_k, ft_k)
                                    p_l, ft_l = cut_zero_imp(p_l, ft_l)

                                popt_x, perr_x = fit_lorentz(p_k, ft_k, fitfunc=fitfunc)
                                popt_y, perr_y = fit_lorentz(p_l, ft_l, fitfunc=fitfunc)
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
    def __init__(self, J_para, J_perp, h, eta, p, dt, filepath, simulation_path, exec_file, runfile, Tc, nr_GPUS=6, nr_Ts=6, size=1024,
                 max_steps=1e9, Ly_Lx = 1/8, equil_error=0.01, equil_cutoff=0.1, T_range_fraction=0.05, T_min_fraction=0.01,
                 max_moving_factor=0.005, min_nr_sites=5e5, para_nr=150, second=False, observed_direction=0, min_corr_nr=5000):
        super().__init__(J_para, J_perp, h, eta, p, dt, filepath, simulation_path, exec_file,
                         nr_GPUS=nr_GPUS, Ly_Lx=Ly_Lx, para_nr=para_nr, runfile=runfile)

        self.nr_Ts = nr_Ts                      # nr of temperatures used to fit
        self.T_range_fraction = T_range_fraction    # This is the fraction of Tc that is used to determine the interval of Ts in [Tc, (1+T_range_raction) * Tc]
        self.T_min_fraction = T_min_fraction
        self.Tc = Tc                            # We need to know the critical temperature that we determined in previous measurements
        # I think we only use one size, being the largest possible
        # Is there some way of guessing how long a simulation will take?
        # 1000 is a save bet but I guess we would like to go up to even larger sizes?
        self.size = size
        self.min_nr_sites = min_nr_sites
        self.max_steps = max_steps
        self.max_moving_factor = max_moving_factor

        self.T_arr = np.array([])
        self.total_runs = 0                 # just a parameter to keep track of how many jobs we want to run in this iteration


        self.all_T_arr = np.array([])       # Bookkeeping for all the temperatures we have simulated in this setting
        self.equil_error = equil_error           # standard equilibration error for the xi runs
        self.maximum_iterations = 2         # we first look at the interval [Tc, 1.05Tc] and If this doesnt work we inrease to [Tc, 1.1Tc] and If this doesnt work we abort
        self.iteration_nr = 0
        self.min_corr_nr = min_corr_nr
        self.corr_write_density = 1 / 100          # We expect a long simulation with very long correlation times and the fit of xi takes a whole lot of time
        self.equil_cutoff = equil_cutoff             # This is the values that we cut off because we think we are still equilibrating. Since we definitely want the values in equilibration we use a relatively large cutoff here
        self.max_time = 0
        self.Tc_fit_tolerance = 0.15        # 5% tolerance for the Tc obtained from the linear regression around the critical point. If its further away, we do not accept the fit
        self.min_r_sqaured = 0.98           # The r²-value of the linear regression should be fairly high so that we can be sure that the interval that we fit is really linear

        self.cur_para_nr = 0                # same method as in Tc?
        self.second = second
        self.observed_direction = observed_direction
    def setup(self):
        # This function will determine the initial T range
        # we want nr_Ts temperatures between
        Tmax = self.Tc * (1 + self.T_range_fraction)
        Tmin = self.Tc * (1 + self.T_min_fraction)
        # I think we can drop Tc itself since it is to close to the critical point
        # and everything that i am trying to calculate here is not accurate and
        # takes a very long time
        self.T_arr = np.linspace(Tmin, Tmax, num=self.nr_Ts)
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
                                                  check_function_args=(self.equil_error, self.equil_cutoff),        # TODO add moving factor
                                                  temp_tolerance=0.002)
        # For every valid simulation we do not have to do this simulation in the following
        for valid_simulation in valid_simulations:
            valid_temp = valid_simulation[1]        # valid simulations is tuple of (size, temp), we only need the temp here
            Ts = list(self.T_arr)                   # change to list as we then can easier use remove
            Ts.remove(valid_temp)                   # remove the valid T from the Ts we still have to do
            self.T_arr = np.array(Ts)               # overwrite the original array with the new Ts
            self.total_runs = len(self.T_arr)

        if self.T_arr.size != 0:
            # if the array is not empty this means that there are still simulations to do
            # write the parameter files
            # I also need to reset this current parameter number here...
            self.cur_para_nr = 0
            self.write_para_files()
            # submit the jobs, should be handlede by the super method
            self.cur_para_nr = 0
            # If the array is empty this means that we have all measurements done already and we can call evaluate
            self.run_jobs()
        # We want to evaluate anyways if we have a valid simulation or if we just ran one
        return self.evaluate()
    def write_para_files(self):
        print("Writing the parameter files...")
        nr_subsystems = np.ceil(self.min_nr_sites / (self.size ** 2 * self.Ly_Lx))
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
                        f"p, {self.p} \n"                        
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
                        f"nr_subsystems, {nr_subsystems} \n"    # The number of subsystems will be one, we use large systems that will run long to eliminate the statistical deivations
                        f"x_y_factor, {self.Ly_Lx} \n"
                        f"nr_corr_values, 0 \n"     # We need a new corr observer that just observes with density and doesnt switch after quench     
                        f"nr_ft_values, 0 \n"       # Ah we still wanted to check whether the values of the ft and fit and python or direct fit in c++ are the same, but they should be fairly similar
                        f"equil_error, {self.equil_error}\n"
                        f"equil_cutoff, {self.equil_cutoff}\n"
                        f"min_corr_nr, {self.min_corr_nr}\n"
                        f"corr_write_density, {self.corr_write_density}\n"
                        f"moving_factor, {self.max_moving_factor}\n"
                        f"corr_second, {int(self.second)}\n"
                        f"observed_direction, {int(self.observed_direction)}")
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
        T_xix, xix_arr = amplitude_measurement.prep_folder_data(self.equil_cutoff, size_path, "xix")
        T_xiy, xiy_arr = amplitude_measurement.prep_folder_data(self.equil_cutoff, size_path, "xiy")

        xix_inv = 1 / xix_arr
        xiy_inv = 1 / xiy_arr
        # Now we have the inverse arrays for the two directions.
        # For the linear regression we need the critical temperature which we can access by self.Tc
        # so we can call the best_fit_inv function that we wrote for the evaluation outside of the suite
        reg_x, T_include_start_x, T_include_end_x = best_fit_inv(T_xix, xix_inv, self.Tc, self.Tc_fit_tolerance, self.min_r_sqaured)
        reg_y, T_include_start_y, T_include_end_y = best_fit_inv(T_xiy, xiy_inv, self.Tc, self.Tc_fit_tolerance,
                                                                 self.min_r_sqaured)
        # The best fit_inv function returns None if we didnt find a fitting segment
        if (reg_x is None) or (reg_y is None):
            return
            print("No fitting fit.")
            fig, ax = plt.subplots(1, 1)
            # First plot the data points
            ax.plot(T_xix,xix_arr, label=rf"$\xi_x$", linestyle="",
                    marker="x")
            ax.plot(T_xiy, xiy_arr, label=rf"$\xi_y$", linestyle="",
                    marker="x")
            ax.set_xlabel("T")
            ax.set_ylabel(r"$\xi$")
            configure_ax(fig, ax)
            plt.show()

            fig, ax = plt.subplots(1, 1)
            # First plot the data points
            ax.plot(T_xix, 1 / xix_arr, label=rf"$1 / \xi_x$", linestyle="",
                    marker="x")
            ax.plot(T_xiy, 1 / xiy_arr, label=rf"$1 / \xi_y$", linestyle="",
                    marker="x")
            ax.set_xlabel("T")
            ax.set_ylabel(r"$\xi$")
            configure_ax(fig, ax)
            plt.show()

            # if we are doing to many iterations we have to return here
            if self.iteration_nr > self.maximum_iterations:
                print("Maximum iterations reached. Aborting")
                return
            else:
                print("We have to repeat the simulation")
            # so we extend the temperatures that we investigate
            stepsize = self.T_arr[1] - self.T_arr[0]
            self.T_arr = np.linspace(np.max(self.T_arr) + stepsize, (1 + 2 * self.T_range_fraction + self.T_min_fraction) * self.Tc, self.nr_Ts)
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
        # Seems this was wrong? But I dont believe it tbh
        # xix_ampl = 1 / reg_x.slope
        # xiy_ampl = 1 / reg_y.slope
        # Tc_x = - reg_x.intercept * xix_ampl
        # Tc_y = - reg_y.intercept * xiy_ampl

        xix_ampl = - 1 / reg_x.intercept
        xiy_ampl = - 1 / reg_y.intercept
        Tc_x = - reg_x.intercept / reg_x.slope
        Tc_y = - reg_y.intercept / reg_y.slope

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
        ax.plot([], [], label=rf"$\xi_x / \xi_y  = {(xix_ampl / xiy_ampl):.2f}$", linestyle="", marker="")
        # like I said if we want them to be in one plot we need to scale the y axix to be logarthimic, you cant do that, the fits will cross zero and so wont look good logarithmic
        ax.set_xlabel("T")
        ax.set_ylabel(r"$1 / \xi$")
        # We set the lower limit of the y axis to be zero since negative xi are not sensible but the fit can become negative
        ax.set_ylim(0, ax.get_ylim()[1])
        configure_ax(fig, ax)
        # saving the plot
        create_directory_if_not_exists(self.simulation_path + "/plots/")
        plt.savefig(self.simulation_path + f"/plots/xi-inv-{self.size}.png", format="png")
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
        plt.savefig(self.simulation_path + f"/plots/T-xi-{self.size}.png", format="png")
        plt.show()

        # I think that is it.
        return


    @staticmethod
    def plot_divergence(simpath, equil_cutoff=0.1, Tc=None, direction="both"):
        """
        member function that top level mäßig does everythin gon its own with a default behavior
        If different behavior is required we have to build out of the constituents funcitons
        It is a staticmethod and is only here form some encapsulation
        :param filepath:
        :return:
        """
        fig, axx = plt.subplots(1, 1, figsize=(6.4 * 1.5, 4.8 * 1.5))

        if direction == "xix":
            results_x = amplitude_measurement.prep_sim_data(equil_cutoff, simpath, "xix")
            results = [results_x]
            axes = [axx]
        elif direction == "xiy":
            results_y = amplitude_measurement.prep_sim_data(equil_cutoff, simpath, "xiy")
            results = [results_y]
            axes = [axx]

        else:
            results_x = amplitude_measurement.prep_sim_data(equil_cutoff, simpath, "xix")
            results_y = amplitude_measurement.prep_sim_data(equil_cutoff, simpath, "xiy")
            results = [results_x, results_y]
            axy = axx.twinx()
            axes = [axx, axy]

        # first of all we want to plot the inverse stuff
        fit_size = 0
        max_size_len = 0
        for j,result in enumerate(results):
            for i, size in enumerate(result):
                # We definitely plot the points
                T_xix, xix_arr = result[size]
                amplitude_measurement.plot_inverse_points(fig, axes[j], T_xix, xix_arr, marker=markers[i], color=5 * j)
                # Find out which size has the largest number of points and for that we do the fit
                if len(xix_arr) > max_size_len:
                    max_size_len = len(xix_arr)
                    fit_size = size

            # The thing is one could actually think about performing the fit on
            # combined measurement points, but that is for another time
            # one would have to replace smart the small size values for the large size
            T_xix, reg_x = amplitude_measurement.perform_fit_on_sizes(
                result, fit_size, Tc)

            amplitude_measurement.plot_xi_inv_fit(axes[j], T_xix, reg_x)

            # Now we still have to configure the inverse plot
        amplitude_measurement.configure_xi_inv_plot(axes, fig)
        create_directory_if_not_exists(simpath + "/plots")
        plt.savefig(simpath + f"/plots/xiy-inv-toplevel.png",
                    format="png")
        plt.show()

        fig, axx = plt.subplots(1, 1, figsize=(6.4 * 1.5, 4.8 * 1.5))
        if direction == "xix":
            axes = [axx]
        elif direction == "xiy":
            axes = [axx]
        else:
            axy = axx.twinx()
            axes = [axx, axy]
        # Now the same stuff for the not inverse part
        for j,result in enumerate(results):
            for i, size in enumerate(result):
                # We definitely plot the points
                T_xix, xix_arr = result[size]

                amplitude_measurement.plot_xi_points(axes[j], T_xix,
                                                     xix_arr, marker=markers[i], size=size, color=5 * j)

                if len(xix_arr) > max_size_len:
                    max_size_len = len(xix_arr)
                    fit_size = size

            # The thing is one could actually think about performing the fit on
            # combined measurement points, but that is for another time
            # one would have to replace smart the small size values for the large size
            T_xix, reg_x = amplitude_measurement.perform_fit_on_sizes(
                result, fit_size, Tc)
            T_xix, xix_arr = result[fit_size]

            amplitude_measurement.plot_xi_fit(axes[j], reg_x, T_xix)

        amplitude_measurement.configure_xi_plot(axes, fig)

        plt.savefig(simpath + f"/plots/T-xi-toplevel.png",
                    format="png")
        plt.show()

    @staticmethod
    def plot_inverse_points(fig, ax, T, xi,
                            marker="s", color=0):
        """
        plots inverse correlation length to an already existing plot
        :param fig:
        :param axx:
        :param axy:
        :param T_xix:
        :param T_xiy:
        :param xix:
        :param xiy:
        :return:
        """
        xix_inv = 1 / xi
        ax.plot(T, xix_inv, label=rf"$1 / \xi_\parallel$",
                 linestyle="", marker=marker, markerfacecolor="none",
                 markeredgecolor=colors[color])


    @staticmethod
    def configure_xi_plot(axes, fig, direction=None):
        # For the limits we want to know the largest plotted point, but how do we
        # get that? I just want to get it out of the axes
        if len(axes) == 2:
            axx, axy = axes[0], axes[1]
            upper_ylim_parallel = get_largest_value_with_linestyle(axx,
                                                                   linestyle='None') * 1.05
            upper_ylim_perp = get_largest_value_with_linestyle(axy,
                                                               linestyle='None') * 1.15
            axx.set_ylim(0, upper_ylim_parallel)
            axy.set_ylim(0,
                         upper_ylim_perp * 1.5)
            axx.set_xlabel("T")
            axx.set_ylabel(r"$\xi_\parallel / a_\parallel$")
            axx.set_xlabel(r"$T~/~J_\parallel$")
            axy.set_ylabel(r"$\xi_\perp / a_\perp$")
            config_x = {
                "labelrotation": 90,
                "labelhorizontalalignement": "right",
                "grid": True,
                "tight_layout": False,
                "legend": False,
                "increasefontsize": 0.75,
            }
            config_y = {
                "labelrotation": 90,
                "labelhorizontalalignement": "right",
                "grid": False,
                "legend": False,
                "increasefontsize": 0.75,
            }
            configure_ax(fig, axx, config_x)
            configure_ax(fig, axy, config_y)
            lines, labels = axx.get_legend_handles_labels()
            lines2, labels2 = axy.get_legend_handles_labels()
            axy.legend(lines + lines2, labels + labels2, loc=0,
                       fontsize=int(PLOT_DEFAULT_CONFIG["legendfontsize"] * (
                               1 + config_x["increasefontsize"])))
        else:
            axx = axes[0]
            upper_ylim_parallel = get_largest_value_with_linestyle(axx,
                                                                   linestyle='None') * 1.05
            axx.set_ylim(0, upper_ylim_parallel)
            axx.set_xlabel("T")
            if direction == "xix":
                axx.set_ylabel(r"$\xi_\parallel / a_\parallel$")
            else:
                axx.set_ylabel(r"$\xi_\perp / a_\perp$")
            axx.set_xlabel(r"$T /$meV")
            config_x = {
                "labelrotation": 90,
                "labelhorizontalalignement": "right",
                "grid": True,
                "tight_layout": False,
                "legend": False,
                "increasefontsize": 0.5,
            }
            configure_ax(fig, axx, config_x)
            lines, labels = axx.get_legend_handles_labels()


    @staticmethod
    def plot_xi_fit(axx, reg_x, T_xix, color=0, direction="parallel", Tc_label=True):
        T_x_plot = np.linspace(np.min(T_xix), np.max(T_xix), 200)
        xix_ampl = - 1 / reg_x.intercept
        Tc_x = - reg_x.intercept / reg_x.slope
        eps_x = (T_x_plot - Tc_x) / Tc_x
        # before we plot we look at the y limits in the case that we dont plot the critical amplitude
        upper_ylim_parallel = axx.get_ylim()[1]
        # The function that we need to use is called critical amplitude
        if Tc_label:
            axx.plot(T_x_plot, critical_amplitude(eps_x, xix_ampl),
                     label=rf"$\xi_\{direction}^+ = {xix_ampl:.2f}, T_c = {Tc_x:.3f}$",
                     color=colors[color])
        else:
            axx.plot(T_x_plot, critical_amplitude(eps_x, xix_ampl),
                     label=rf"$\xi_\{direction}^+ = {xix_ampl:.2f}$",
                     color=colors[color])


    @staticmethod
    def plot_xi_points(axx, T_xix, xix_arr, marker="s", size="", color=0, direction="parallel"):
        axx.plot(T_xix, xix_arr, label=rf"$\xi_\{direction}" + "^{" + str(size) + "}$", linestyle="",
                 marker=marker, markerfacecolor="none",
                 markeredgecolor=colors[color])


    @staticmethod
    def configure_xi_inv_plot(axes, fig, direction=None):
        if len(axes) == 2:
            axx, axy = axes
            axx.set_ylabel(r"$(\xi_\parallel / a_\parallel)^{-1}$")
            axx.set_xlabel(r"$T /$meV")
            axy.set_ylabel(r"$(\xi_\perp / a_\perp)^{-1}$")
            axy.set_xlabel(r"$T /$ meV")
            # Okay the lower y limit should be zero for both
            axx.set_ylim(0, axx.get_ylim()[
                1] * 4 / 3)  # the x values  are smaller so they should look smaller, increasing the upper bound by one third?
            axy.set_ylim(0, axy.get_ylim()[1])
            axx.set_title(r"Inverse Correlation Length $1 / \xi$")
            axy.set_title(r"Inverse Correlation Length $1 / \xi$")
            config_x = {
                "labelrotation": 90,
                "labelhorizontalalignement": "right",
                "grid": True,
                "tight_layout": False,
                "legend": False,
                "increasefontsize": 0.5,
            }
            config_y = {
                "labelrotation": 90,
                "labelhorizontalalignement": "right",
                "grid": False,
                "legend": False,
                "increasefontsize": 0.5,
            }
            configure_ax(fig, axx, config_x)
            configure_ax(fig, axy, config_y)
            # We need to deal with the leglend ourself I think
            lines, labels = axx.get_legend_handles_labels()
            lines2, labels2 = axy.get_legend_handles_labels()
            axy.legend(lines + lines2, labels + labels2, loc=0, fontsize=int(
                PLOT_DEFAULT_CONFIG["legendfontsize"] * (
                        1 + config_x["increasefontsize"])))
        else:
            axx = axes[0]
            if direction == "xix":
                axx.set_ylabel(r"$(\xi_\parallel / a_\parallel)^{-1}$")
            else:
                axx.set_ylabel(r"$(\xi_\perp / a_\perp)^{-1}$")
            axx.set_xlabel(r"$T /$meV")
            # Okay the lower y limit should be zero for both
            axx.set_ylim(0, axx.get_ylim()[
                1] * 4 / 3)  # the x values  are smaller so they should look smaller, increasing the upper bound by one third?
            axx.set_title(r"Inverse Correlation Length $1 / \xi$")
            config_x = {
                "labelrotation": 90,
                "labelhorizontalalignement": "right",
                "grid": True,
                "tight_layout": False,
                "legend": False,
                "increasefontsize": 0.5,
            }

            configure_ax(fig, axx, config_x)
            # We need to deal with the leglend ourself I think
            lines, labels = axx.get_legend_handles_labels()

    @staticmethod
    def perform_fit_on_sizes(results_x, fit_sizes, Tc=None, min_points=0):
        if isinstance(fit_sizes, list):
            T_xix = np.array([])
            xix_arr = np.array([])
            for size in fit_sizes:
                cur_T, cur_xix = results_x[size]
                T_xix = np.concatenate((T_xix, cur_T))
                xix_arr = np.concatenate((xix_arr, cur_xix))
            xix_inv = 1 / xix_arr
        else:
            T_xix, xix_arr = results_x[fit_sizes]
            xix_inv = 1 / xix_arr
        if not Tc:
            Tc = np.min(T_xix)  # then we guess it to be this
        reg_x, T_include_start_x, T_include_end_x = best_fit_inv(T_xix, xix_inv,
                                                                 Tc, 0.5,
                                                                 min_r_squared=0.9, min_points=min_points)
        return T_xix, reg_x

    @staticmethod
    def plot_xi_inv_fit(axx, T_xix, reg_x):
        """
        does what the name says
        :param axx:
        :param axy:
        :param T_xiy:
        :param T_xix:
        :param reg_x:
        :param reg_y:
        :return:
        """

        xix_ampl = - 1 / reg_x.intercept
        Tc_x = - reg_x.intercept / reg_x.slope
        # Plot the fit
        axx.plot(T_xix, reg_x.intercept + reg_x.slope * T_xix,
                 label=rf"$\xi_\parallel^+ = {xix_ampl:.2f}, T_c = {Tc_x:.3f}$",
                 color=colors[0])


    @staticmethod
    def prep_sim_data(equil_cutoff, simpath, value):
        results_x = {}
        size_folders = find_size_folders(simpath)
        for size_folder in size_folders:
            if (size_folder != "plots") & (size_folder[0] != "."):
                size_folder_path = os.path.join(simpath, size_folder)
                if os.path.isdir(size_folder_path):
                    T_xix, xix_arr = amplitude_measurement.prep_folder_data(equil_cutoff, size_folder_path, value)
                    results_x[int(size_folder)] = (T_xix, xix_arr)
        return results_x

    @staticmethod
    def prep_folder_data(equil_cutoff, folderpath, value):
        xix_dic = process_size_folder(folderpath, threshold=equil_cutoff, key="T",
                                      value=value,
                                      file_ending="corr")  # the selected temperatures are just all temperatures

        T_xix = xix_dic["T"]
        xix_arr = xix_dic[value]
        # They should be sorted and we are actually interested in the inverse correlation lengths
        xix_arr = xix_arr[np.argsort(T_xix)]
        T_xix = T_xix[np.argsort(T_xix)]
        return T_xix, xix_arr

    def get_para_nr(self):
        # this tactic inreases the parameter number everytime get_para_nr is called so that we do not submit any job twice
        self.cur_para_nr += 1
        return self.para_nr + self.cur_para_nr - 1

class z_measurement(autonomous_measurement):
    def __init__(self, J_para, J_perp, h, eta, p, dt, filepath, simulation_path, exec_file, test_exec_file, runfile, Tc, nr_GPUS=6, size_min=64,
                          size_max=256, nr_sizes=3, max_steps=1e9, nr_sites=5e5, Ly_Lx = 1/8, equil_error=0.005, equil_cutoff=0.5,
                 variation_error_rate=0.011, z_guess=2, min_nr_sites=1e6, min_nr_systems=100, fold=50, cum_density=1/100,
                 test_cum_density=1/2, test_min_cum_nr=2000):
        # call the constructor of the parent classe
        super().__init__(J_para, J_perp, h, eta, p, dt, filepath, simulation_path, exec_file,
                         nr_GPUS=nr_GPUS, Ly_Lx=Ly_Lx, runfile=runfile)
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
                                                  check_cum_variation_valid,
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
                                                  check_cum_variation_valid,
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
                        f"p, {self.p} \n"
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
                        f"p, {self.p} \n"                        
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
    nr_gpus = 20
    # we somehow need the relevant parameters
    # The model defining parameters are J_perp J_para h eta
    # the simulation defining parameters are dt
    J_para = -110000
    #J_para = -9
    J_perp = -3500
    #J_perp = -0.1
    h = 10000
    #h = 0.5
    eta = 1
    p = 2.5
    dt = 0.00001
    # dt = 0.01
    max_size_Tc = 128
    min_size_Tc = 64
    nr_sizes_Tc = 2
    nr_Ts = 10
    cum_error = 0.0015
    equil_cutoff_Tc = 0.5
    value_name = "U_L"
    file_ending = "mag"
    process_file_func = recalculate_mag_file_to_U_L
    value_write_density = 0.001


    random_init = 0.0
    #filepath = "/home/weitze73/Documents/Master-Arbeit/Code/Master-Arbeit/CudaProject"
    filepath = "/home/andi/Studium/Code/Master-Arbeit/CudaProject"
    simulation_path = "../../Generated content/Final/CriticalTemperature/"

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
    runfile = "run_cuda_gpu_a100_low.sh"

    # T- parameters?
    max_rel_intersection_error = 0.02
    min_cum_nr = 5000
    moving_factor = 0.001
    T_min = 29071.961123
    T_max = 31494.624550

    # Quench parameters
    max_size = 1024
    min_nr_sites = 1e6


    # Amplitude parameters
    amplitude_size = 1024
    equil_error = 0.05
    equil_cutoff = 0.01

    # z parameters
    size_min = 64
    size_max = 256
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
                                    value_write_density=value_write_density, runfile=runfile,
                                    T_min=T_min, T_max=T_max, para_nr=para_nr)
        T_c, T_c_error = sim.routine()
    elif measurements["efficient Tc"]:
        sim = efficient_crit_temp_measurement(J_para, J_perp, h, eta, p, dt, filepath, simulation_path + "Tc", Tc_exec_file, nr_GPUS=nr_gpus,
                                              size_min=min_size_Tc, size_max=max_size_Tc,
                                              intersection_error=max_rel_intersection_error, random_init=random_init,
                                              max_moving_factor=moving_factor, min_val_nr=min_cum_nr)
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
                                     equil_error=equil_error, equil_cutoff=equil_cutoff)
        ampl.run()
    if measurements["z"]:
        z_measure = z_measurement(J_para, J_perp, h, eta, p, dt, filepath, simulation_path + "z",
                                        z_exec_file, z_test_exec_file, T_c, nr_GPUS=nr_gpus, size_min=64, size_max=256,
                                        nr_sizes=nr_sizes, min_nr_sites=z_min_nr_sites, min_nr_systems=z_min_nr_systems,
                                     equil_error=z_equil_error)
        z_measure.run()

if __name__ == '__main__':
    main()