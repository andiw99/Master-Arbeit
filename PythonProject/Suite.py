import numpy as np
from scipy.optimize import fsolve
from itertools import product
import re
from fabric import Connection
import time
from FunctionsAndClasses import *

def T_c_XY(T, J_parallel, J_perp):
    return 2 * T / J_parallel * np.log(2 * T / J_perp) - 1

def T_c_est(J_para, J_perp, h):
    # we need to solve the transcendet equation of the XY model
    # scipy f_solve needs a guess of the critical temperature which is difficult but lets try with the mean of the Js
    T_c_est = (J_para + J_perp) / 2
    T_c_est = fsolve(T_c_XY, T_c_est, args=(J_para, J_perp))
    return T_c_est


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

def start_cum_measurement(J_para, J_perp, h, eta, dt, filepath, simulation_path, nr_GPUS=6, nr_Ts=5, size_min=48,
                          size_max=80, nr_sizes=3, max_steps=1e9, nr_sites=5e6, Ly_Lx = 1/8):
    # this somehow needs the parameters, where do we put them? In a file? On the moon? User input?
    T_min = T_c_est(J_para, J_perp, h)
    # TODO really crude approximation of the maximum T, I really should think about something better
    T_max = T_min + (h/J_perp) * T_min
    # We use nr_Ts datapoints
    T_arr = np.linspace(T_min, T_max, nr_Ts)
    sizes = np.linspace(size_min, size_max, nr_sizes, endpoint=True, dtype=np.int32)
    max_time = dt * max_steps
    total_runs = nr_sizes * nr_Ts   # every temp size combo will be a seperate simulation
    # okay weve written the parameter file, what do we do now?
    # now I have to think about how to commit the jobs and how to keep track of runnign jobs
    for i, (T, size) in enumerate(product(T_arr, sizes)):
        # We now need to construct the parameterfile with the appropriate temperature and size
        # Is it okay if we construct all files in the beginning and deal with the threading of the gpus later?
        # to construct the para set we need to know how many subsystems we should initialize
        nr_subsystems = int(nr_sites / (size ** 2 * Ly_Lx))
        with open(filepath + "/para_set_100" + str(i) + '.txt', 'w') as f:
            f.write(simulation_path)
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
                    f"nr_subsystem_sizes, 0  \n"
                    f"nr_subsystems, {nr_subsystems} \n"
                    f"x_y_factor, {Ly_Lx} \n"
                    f"nr_corr_values, 0 \n"
                    f"nr_ft_values, 0 \n")
    running_jobs = []   # this is a list that will have all the job ids
    # actually I think this should probably all be a class



class crit_temp_measurement():
    def __init__(self, J_para, J_perp, h, eta, dt, filepath, simulation_path, nr_GPUS=6, nr_Ts=5, size_min=48,
                          size_max=80, nr_sizes=3, max_steps=1e9, nr_sites=5e6, Ly_Lx = 1/8):
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
    def init(self):
        # this somehow needs the parameters, where do we put them? In a file? On the moon? User input?
        T_min = T_c_est(self.J_para, self.J_perp, self.h)
        # TODO really crude approximation of the maximum T, I really should think about something better
        T_max = T_min + (self.h / self.J_perp) * T_min
        # We use nr_Ts datapoints
        self.T_arr = np.linspace(T_min, T_max, self.nr_Ts)
        self.sizes = np.linspace(self.size_min, self.size_max, self.nr_sizes, endpoint=True,
                            dtype=np.int32)
        self.max_time = self.dt * self.max_steps
        self.total_runs = self.nr_sizes * self.nr_Ts  # every temp size combo will be a seperate simulation
    def routine(self, walltime='24:00:00', file='SubsystemRelaxation.cu', folder='simulations', user='weitze73', wait=20):
        """
        the outer routine for the T_c calculation
        :param walltime: walltime per job
        :param file: file to be compiled, we may want to have a seperate file for the autonomous exection
        :param folder: folder for the executables
        :param wait: sleep for wait seconds before checking again if a job is finished        
        :return: -
        """
        self.init()                 # initializing determines the T range
        self.iteration(file, folder, user, wait, walltime)

    def iteration(self, file, folder, user, wait, walltime):
        para_nr = self.write_para_files()  # setting up the parameter files for every simulation
        self.run_jobs(file, folder, para_nr, user, wait, walltime)
        # after running the jobs, we need to
        # calculate the binder cumulant for every run
        threshold = 0.1  # in the simulation we calculated the mean from the last 90% of the values and achieved a small error
        results = {}
        for size_folder in os.listdir(self.simulation_path):
            if (size_folder[0] != ".") & (size_folder != "plots"):
                size_folder_path = os.path.join(self.simulation_path,
                                                size_folder)
                if os.path.isdir(size_folder_path):
                    size_result = process_size_folder(size_folder_path,
                                                      threshold)
                    results[int(size_folder)] = size_result
        # how does results look again precisely?
        # I have for every size a dictionary with {'T' : [T_min, ..., T_max], 'U_L': [U_L(T_min), ..., U_L(T_max)]}
        # once we have the U_L values the first thing we should do is check whether we have an intersection
        # something should be recursively if we need a second or third iteration...
        U_L_min_T_min = np.min(results[np.min(self.sizes)]['U_L'])
        U_L_max_T_min = np.min(results[np.max(self.sizes)]['U_L'])
        U_L_min_T_max = np.max(results[np.min(self.sizes)]['U_L'])
        U_L_max_T_max = np.max(results[np.max(self.sizes)]['U_L'])
        # we say we have an intersection if U_L_min_T_min > U_L_max_T_min
        # and U_L_min_T_max < U_L_max_T_max
        intersection = (U_L_min_T_min > U_L_max_T_min) & (
                    U_L_min_T_max < U_L_max_T_max)
        if intersection:
            # good sign
            # determine Tc AND its error
            T_range, U_L_intersection, T_intersection, U_L_interpolated = interpolate_and_minimize(
                results)
            # I think T_range is the common Ts for the sizes, I think here every size should
            # definetely have the same temperatures
            # U_L_intersection is the U_L_value at the intersection
            # T_intersection is the temperature at the intersection
            # U_L_interpolated is a list with the numerical interpolated values at the stützstellen with resolution of
            # 1000, the resolution can be 
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
                T_min = np.maximum(np.min(self.T_arr) - T_range, 0.0)   # We should not consider negative temperatures
                T_max = np.max(self.T_arr) - (self.T_arr[1] - self.T_arr[0])
                self.T_arr = np.linspace(T_min, T_max, self.nr_Ts)
            self.iteration(file, folder, user, wait, walltime)


    def run_jobs(self, file, folder, para_nr, user, wait, walltime):
        # how do we make this routine? First we can make one cycle and submit nr gpus jobs?
        # or just the routine that will be waiting for the jobs to finish instantly?
        # establish the connection to hemera
        self.connection = Connection('hemera')
        next_job = 0
        while (len(self.completed_jobs)) < self.total_runs:
            # after this the set of running jobs is guaranteed to be empty
            # now we should check wheter some jobs are completed
            # just by getting which jobs are pending etc?
            queue_command = f"squeue -u {user}"
            queue_feedback = self.connection.run(queue_command)
            jobs_on_hemera = extract_numbers_before_newline(
                queue_feedback.stdout)

            # now we just check if all ids in running jobs are still in jobs_on_hemera
            # if not, we check if the id is completed
            # if thats true, we move the job id to the completed jobs

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
                        self.running_jobs.remove(job_id)
                        self.completed_jobs.add(job_id)
                    else:
                        # if it is not completed but not in jobs_on_hemera anymore,
                        # we have a problem
                        print(
                            "Oups! The Job vanished. Please check what happend")

            # now we need to determine how many jobs are currently running
            # and we need to know how to submit jobs
            # while the number running jobs is smaller than the number of GPUs
            # to use, we submit new jobs
            # indeed it would be optimal if the long running jobs would start first
            # but that is hard to do? For now we just submit the jobs in arbitrary
            # order
            # we know which job is next through the next job variable,
            # but we still need to know the number at which to start
            while len(self.running_jobs) < self.nr_GPUS:
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
                # I think we can be sure that the job is running if we just commited it
                # better double check?
                self.running_jobs.add(job_id)
            # now we just wait some time before we do the next check?
            time.sleep(wait)
        # if we are here that means that all runs are done
        # we add the currently simulated temperatures to the bookkeeping variable
        self.all_T_arr = np.concatenate((self.all_T_arr, self.T_arr))

    def write_para_files(self, para_nr=100):
        for i, (T, size) in enumerate(product(self.T_arr, self.sizes)):
            # We now need to construct the parameterfile with the appropriate temperature and size
            # Is it okay if we construct all files in the beginning and deal with the threading of the gpus later?
            # to construct the para set we need to know how many subsystems we should initialize
            nr_subsystems = int(self.nr_sites / (size ** 2 * self.Ly_Lx))
            with open(self.filepath + "/para_set_" + str(para_nr) + str(i) + '.txt', 'w') as f:
                f.write(self.simulation_path)
                f.write(f"end_time, {self.max_time} \n"
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
                        f"nr_ft_values, 0 \n")

        return para_nr
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