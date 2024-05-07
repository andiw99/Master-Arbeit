//
// Created by andi on 25.02.24.
//

#ifndef CUDAPROJECT_MAIN_CUDA_CUH
#define CUDAPROJECT_MAIN_CUDA_CUH
#include "../LearningProject/Header/Helpfunctions and Classes.h"
#include <functional>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <boost/typeof/typeof.hpp>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/random.h>
#include <thrust/reduce.h>
#include <thrust/transform.h>
#include <thrust/transform_reduce.h>
#include <thrust/sequence.h>
#include <thrust/execution_policy.h>
#include <map>
#include <boost/asio/ip/host_name.hpp>
#include <cufft.h>
#include <numeric>
#include <thrust/generate.h>
#include <thrust/for_each.h>
#include <thrust/execution_policy.h>
#include <thrust/random.h>
using namespace std;

enum Parameter {
    dim_size_x, dim_size_y, total_size, J, Jy,dt, end_time, start_time, starting_temp, end_temp, T, alpha, beta, tau,
    eta, nr_saves, nr_repeat, min_temp, max_temp, nr_runs, K, tol, random_init, x0, p0, min_tau_factor, max_tau_factor,
    equil_time, logspaced, step_nr, run_nr, min_lat_factor, max_lat_factor, nr_ner_values, m, p, curand_random,
    subsystem_Lx, subsystem_Ly, subsystem_min_Lx, subsystem_max_Lx, nr_subsystem_sizes, nr_subsystems, x_y_factor,
    nr_cum_values, nr_mag_values, nr_corr_values, nr_ft_values, equil_error, min_cum_nr, min_corr_nr, equil_cutoff,
    cum_write_density, corr_write_density, ft_write_density, moving_factor, min_mag_nr, mag_write_density, corr_second,
    observed_direction, equil_time_end, gamma,
};

map<Parameter, string> parameter_names {
        {dim_size_x, "dim_size_x"}, {dim_size_y, "dim_size_y"}, {total_size, "total_size"}, {J, "J"}, {Jy,"Jy"},
        {dt,"dt"}, {end_time,"end_time"}, {start_time,"start_time"}, {starting_temp,"starting_temp"},
        {end_temp,"end_temp"}, {T,"T"}, {alpha,"alpha"}, {Parameter::beta, "beta"}, {tau,"tau"}, {eta,"eta"},
        {nr_saves,"nr_saves"}, {nr_repeat,"nr_repeat"}, {min_temp,"min_temp"}, {max_temp,"max_temp"},
        {nr_runs,"nr_runs"}, {K,"K"}, {tol,"tol"}, {random_init,"random_init"}, {x0,"x0"}, {p0,"p0"},
        {min_tau_factor,"min_tau_factor"}, {max_tau_factor,"max_tau_factor"}, {equil_time,"equil_time"}, {
        logspaced,"logspaced"}, {step_nr,"step_nr"}, {run_nr,"run_nr"}, {min_lat_factor,"min_lat_factor"},
        {max_lat_factor,"max_lat_factor"}, {nr_ner_values, "nr_ner_values"}, {m, "m"}, {p, "p"},
        {curand_random, "curand_random"}, {subsystem_Lx, "subsystem_Lx"}, {subsystem_Ly, "subsystem_Ly"},
        {subsystem_min_Lx,"subsystem_min_Lx"}, {subsystem_max_Lx,"subsystem_max_Lx"},
        {nr_subsystem_sizes,"nr_subsystem_sizes"}, {nr_subsystems,"nr_subsystems"}, {x_y_factor,"x_y_factor"},
        {nr_cum_values, "nr_cum_values"}, {nr_mag_values, "nr_mag_values"}, {nr_corr_values, "nr_corr_values"},
        {nr_ft_values, "nr_ft_values"}, {equil_error, "equil_error"}, {min_cum_nr, "min_m_nr"},
        {equil_cutoff, "equil_cutoff"}, {cum_write_density, "cum_write_density"}, {corr_write_density, "corr_write_density"},
        {min_corr_nr, "min_corr_nr"}, {moving_factor, "moving_factor"}, {mag_write_density, "mag_write_density"},
        {min_mag_nr, "min_mag_nr"}, {ft_write_density, "ft_write_density"}, {corr_second, "corr_second"},
        {observed_direction, "observed_direction"}, {equil_time_end, "equil_time_end"}, {gamma, "gamma"},


};

map<string, Parameter> string_to_parameter {
        {"dim_size_x", dim_size_x}, {"dim_size_y", dim_size_y}, {"total_size", total_size}, {"J", J}, {"Jy", Jy},
        {"dt", dt}, {"end_time", end_time}, {"start_time", start_time}, {"starting_temp", starting_temp},
        {"end_temp", end_temp}, {"T", T}, {"alpha", alpha}, {"beta", Parameter::beta}, {"tau", tau}, {"eta", eta},
        {"nr_saves", nr_saves}, {"nr_repeat", nr_repeat}, {"min_temp", min_temp}, {"max_temp", max_temp},
        {"nr_runs", nr_runs}, {"K", K}, {"tol", tol}, {"random_init", random_init}, {"x0", x0}, {"p0", p0},
        {"min_tau_factor", min_tau_factor}, {"max_tau_factor", max_tau_factor}, {"equil_time", equil_time},
        {"logspaced", logspaced}, {"step_nr", step_nr}, {"run_nr", run_nr}, {"min_lat_factor", min_lat_factor},
        {"max_lat_factor", max_lat_factor}, {"m", m}, {"p", p}, {"nr_ner_values", nr_ner_values},
        {"curand_random", curand_random}, {"subsystem_Lx", subsystem_Lx}, {"subsystem_Ly", subsystem_Ly},
        {"subsystem_min_Lx", subsystem_min_Lx}, {"subsystem_max_Lx", subsystem_max_Lx},
        {"nr_subsystem_sizes", nr_subsystem_sizes}, {"nr_subsystems", nr_subsystems}, {"x_y_factor", x_y_factor},
        {"nr_cum_values", nr_cum_values}, {"nr_mag_values", nr_mag_values}, {"nr_corr_values", nr_corr_values},
        {"nr_ft_values", nr_ft_values}, {"equil_error", equil_error}, {"min_m_nr", min_cum_nr},
        {"equil_cutoff", equil_cutoff}, {"cum_write_density", cum_write_density}, {"corr_write_density", corr_write_density},
        {"min_corr_nr", min_corr_nr}, {"moving_factor", moving_factor}, {"mag_write_density", mag_write_density},
        {"min_mag_nr", min_mag_nr}, {"ft_write_density", ft_write_density}, {"corr_second", corr_second},
        {"observed_direction", observed_direction}, {"equil_time_end", equil_time_end}, {"gamma", gamma},


};

std::map<Parameter, double> readTxtFileToParameterMap(const path& filename, int startLine = 1) {
    std::map<Parameter, double> paramMap;
    std::ifstream file(filename);

    if (!file.is_open()) {
        std::cerr << "Error: Could not open the file." << std::endl;
        return paramMap;
    }
    int currentLine = 0; // Track the current line number
    std::string line;
    while (std::getline(file, line)) {
        // Skip lines until the desired startLine is reached
        currentLine++;
        if (currentLine < startLine) {
            continue;
        }


        std::istringstream iss(line);
        std::string parameter, value;

        if (std::getline(iss, parameter, ',') && std::getline(iss, value)) {
            paramMap[string_to_parameter[parameter]] = stod(value);
        }

    }

    file.close();
    return paramMap;
}

string readTxtFileToString(const std::string& filepath, int lineToRead = 1) {
    std::ifstream file(filepath);
    std::string line;
    int currentLine = 0; // Track the current line number

    if (file.is_open()) {
        while (currentLine < lineToRead && std::getline(file, line)) {
            currentLine++;
        }
        file.close();
    }

    return line;
}

/*
 * Now we need a "Thrust Algebra", will this be doable? I still need to figure out what my statetype is and how I
 * iterate over the lattice. I hope in the coupled oscillator example this will get clear
 */
struct thrust_algebra {
    template<class S1, class S2, class S3, class S4, class Op>
    static void for_each(S1 &s1, S2 &s2, S3 &s3, S4 &s4, Op op) {
        thrust::for_each(
                thrust::make_zip_iterator( thrust::make_tuple(
                        s1.begin(), s2.begin(), s3.begin(), s4.begin()) ),
                thrust::make_zip_iterator( thrust::make_tuple(
                        s1.end(), s2.end(), s3.end(), s4.begin())),
                op);
    }
    template<class S1, class S2, class S3, class Op>
    static void for_each(S1 &s1, S2 &s2, S3 &s3, Op op) {
        thrust::for_each(
                thrust::make_zip_iterator( thrust::make_tuple(
                        s1.begin(), s2.begin(), s3.begin()) ),
                thrust::make_zip_iterator( thrust::make_tuple(
                        s1.end(), s2.end(), s3.end())),
                op);
    }

    template<class S1, class S2, class S3, class Op>
    static void for_each(S1 s1, S2 s2, S3 s3, size_t n, Op op) {
        thrust::for_each(
                thrust::make_zip_iterator( thrust::make_tuple(
                        s1, s2, s3) ),
                thrust::make_zip_iterator( thrust::make_tuple(
                        s1 + n, s2 + n, s3 + n)),
                op);
    }

    template<class S1, class S2, class Op>
    static void for_each(S1 &s1, S2 &s2, Op op) {
        thrust::for_each(
                thrust::make_zip_iterator( thrust::make_tuple(
                        s1.begin(), s2.begin()) ),
                thrust::make_zip_iterator( thrust::make_tuple(
                        s1.end(), s2.end())),
                op);
    }
};



/*
 * Now thrust operations that have to apply em with only one factor but we have to make a sqrt transform somewhere
 */
struct thrust_operations {
    template<class time_type = double>
    struct apply_drift {
        const time_type dt;
        apply_drift(time_type dt)
                : dt(dt) { }
        template< class Tuple >
        __host__ __device__ void operator()(Tuple tup) const {
            thrust::get<0>(tup) = thrust::get<1>(tup) +
                                  dt * thrust::get<2>(tup);

        }
    };

    template<class time_type = double>
    struct apply_bbk_v1 {
        const time_type pref, dt;
        apply_bbk_v1(time_type dt, time_type eta): dt(dt), pref(1.0 - 0.5 * eta * dt){}
        template< class Tuple >
        __host__ __device__ void operator()(Tuple tup) const {
            thrust::get<0>(tup) = pref * thrust::get<0>(tup) +             // v^n
                                  0.5 *  (dt * thrust::get<1>(tup) +        // F
                                          sqrt(dt) * thrust::get<2>(tup));      // theta

        }
    };

    template<class time_type = double>
    struct apply_bbk_v2 {
        const time_type pref, dt;
        apply_bbk_v2(time_type dt, time_type eta): dt(dt), pref(1.0  / (1.0 + 0.5 * eta * dt)){}
        template< class Tuple >
        __host__ __device__ void operator()(Tuple tup) const {
            thrust::get<0>(tup) = pref * (thrust::get<0>(tup) +             // v^n
                                          0.5 *  (dt * thrust::get<1>(tup) +            // F
                                                  sqrt(dt) * thrust::get<2>(tup)));      // theta

        }
    };

};


string construct_filename(int run_nr) {
    string job_id_str;
    if (const char* job_id = std::getenv("SLURM_JOB_ID")){
        job_id_str = string(job_id);
    } else {
        job_id_str = "local";
    }
    return to_string(run_nr) + "-" + get_current_time().substr(0, 4) + "-" + job_id_str + "-" + boost::asio::ip::host_name();
}

struct checkpoint_timer{
    map<string, long> times;
    map<string, chrono::time_point<chrono::high_resolution_clock>> startpoints;
    map<string, chrono::time_point<chrono::high_resolution_clock>> endpoints;
public:
    // we now make a function, that sets the startpoint of checkpoint i
    void set_startpoint(const string& checkpoint) {
        startpoints[checkpoint] = chrono::high_resolution_clock::now();
    }

    void set_endpoint(const string& checkpoint) {
        endpoints[checkpoint] = chrono::high_resolution_clock::now();
        // add the duration
        times[checkpoint] += std::chrono::duration_cast<std::chrono::microseconds>(
                endpoints[checkpoint] - startpoints[checkpoint]).count();
    }

    map<string, long> get_times() {
        return times;
    }

    ~checkpoint_timer() {
        // print the durations
        for(const auto& timepair : times) {
            cout << timepair.first << " took " << (double)timepair.second * 0.001 << " ms" << endl;
        }
    }
};

class Singleton_timer {
    // this is supposed to be a singleton timer, but I dont know how this will work with members
    // I want to have a checkpoint timer as member.
    // does this member have to be static? We try first not to because I dont know which implications this has
public:
    Singleton_timer(const Singleton_timer&) = delete;   // this prevents the Singleton from being copied

    static Singleton_timer& Get() {
        static Singleton_timer instance;        // I don't really get what happens here, but something with the timer being instantiated the first time Get is called?
        return instance;
    }
    // now the static functionality that is internally handled by the checkpoint timer
    static void set_startpoint(const string& checkpoint) {
        Get().Iset_startpoint(checkpoint);      // until now it does not complain?
    }
    static void set_endpoint(const string& checkpoint) {
        Get().Iset_endpoint(checkpoint);
    }

    static void write(fs::path simulation_path, int file_nr) {
        // again this should work as a static file since i am getting the instance by Get()
        // and this instance then knows the checkpoint timer
        Get().Iwrite(simulation_path, file_nr);
    }
private:
    // private constructor
    Singleton_timer() {}

    // a checkpoint timer as private member that handles the calculation of the timings
    checkpoint_timer Timer = checkpoint_timer();        // I just construct it here?

    void Iset_startpoint(const string& checkpoint) {
        Timer.set_startpoint(checkpoint);
    }

    void Iset_endpoint(const string& checkpoint) {
        Timer.set_endpoint(checkpoint);
    }

    void Iwrite(fs::path simulation_path, int file_nr) {
        // I think I need to get the file nr from outside
        fs::path filepath;
        string filename = construct_filename(file_nr);
        if(fs::is_directory(simulation_path)) {
            // If it is a directory we have th normal case so not pickup
            filepath = simulation_path / (filename + ".perf" );
        } else {
            // pickup so filepath is the simulation parent path (folder to the simulation and the filename)
            filepath = simulation_path.parent_path() / (filename + ".perf");
        }

        ofstream ofile;
        ofile.open(filepath);
        if(ofile.is_open()) {
            cout << "Singleton timer successfully opened file";
        } else {
            cout << "Singleton timer failed to open file";
        }
        // sure the folder is already created but it cannot harm
        // now we can write, i think we stay simple for now?
        // write the header of the file
        ofile << "Part,ms,s" << endl;
        // We implement the write logic here i would say, so lets get the times from the timer
        map<string, long> times = Timer.get_times();
        // so for now we just write the name and the time it took in... ms?
        // or both in ms and in seconds so it is more readable
        // the longs in times are in microseconds
        for(const auto& name_time : times) {
            int ms = (int)((double)name_time.second * 1e-3);        // TODO does this work with 1e-3?
            // since it so hard to f*ing format strings that you write to a file we do some rounding
            // we drop the last digit of the ms to get 10ms as the last digit and then move the dot by two places
            double s = (double)((int)(ms * 1e-1) * 1e-2);
            // I would like to format this but somehow this is impossible in c++
            ofile << name_time.first << ":" << ms << ", " << s << endl;
        }
        // once we have written we want to reset the times, so we just replace the timer?
        Timer = checkpoint_timer();     // I mean this calls the constructor so this should be a new timer

    }

    // Now we should be able to use this Singleton instance everywhere where main is imported?

};

struct rand_normal_values
{
    double mu, sigma, ampl;

    __host__ __device__
    rand_normal_values(double ampl, double mu = 0.0, double sigma = 1.0) : ampl(ampl), mu(mu), sigma(sigma) {};

    __host__ __device__
    float operator()(const unsigned int ind) const
    {
        thrust::default_random_engine rng;
        thrust::normal_distribution<double> dist(mu, sigma);
        rng.discard(ind);

        // dist * ampl zu returnen ist wie... aus dist mit std ampl zu ziehen: b * N(m, o) = N(m, b*o)
        double rnr = (double)dist(rng);
        // if the drawn number is greater than 3 * sigma we cut it
/*        printf(
        "%f  ", rnr
        );
        printf("\n");*/
        if(abs(rnr) > 2.0 * sigma) {
            rnr = (rnr < 0) ? (-1.0) * 2.0 * sigma : 2.0 * sigma;
/*            printf(
                    "rnr changed to %f", rnr
                    );*/
        }
        return (double)ampl * rnr;
    }
};

struct rand_uni_values
{
    float xmin, xmax;

    __host__ __device__
    rand_uni_values(float xmin, float xmax) : xmin(xmin), xmax(xmax) {};

    __host__ __device__
    float operator()(const unsigned int ind) const
    {
        thrust::default_random_engine rng;
        thrust::uniform_real_distribution<float> dist(xmin, xmax);
        rng.discard(ind);
        auto rnr = dist(rng);
        return rnr;
    }
};

struct thrustDivideBy
{
    double dividor;
    thrustDivideBy(double dividor): dividor(dividor) {}

    __host__ __device__
    double operator()(double x) const
    {
        return x / dividor;
    }
};

template <class value_type>
__host__ __device__
inline value_type positive_modulo(value_type i, value_type n) {
    return fmod((fmod(i, n) + n), n);
}

std::pair<std::vector<double>, double*> cut_data_around_peak(
        const std::vector<double>& x_values,
        double* y_values,
        double threshold_fraction = 0.1,
        double min_points_fraction = 0.2
) {
    size_t L = x_values.size();
    // Find the minimum y_value (offset)
    double y_offset = *std::min_element(y_values, y_values + L);

    // Remove the offset from y_values
    for (size_t i = 0; i < L; ++i) {
        y_values[i] -= y_offset;
    }

    // Find the index of the peak
    size_t peak_index = std::max_element(y_values, y_values + L) - y_values;
    double peak_value = y_values[peak_index];

    // Calculate the threshold based on the fraction of the peak value
    double threshold = threshold_fraction * peak_value;

    // Find the range of indices where the y-values are above the threshold
    // basically np where
    std::vector<size_t> above_threshold_indices;
    for (size_t i = 0; i < L; ++i) {
        if (!(y_values[i] < threshold)) {
            above_threshold_indices.push_back(i);
        }
    }
    // If the length is too small, reduce the threshold fraction
    if (above_threshold_indices.size() < min_points_fraction * L) {
        double new_threshold = 0.9 * threshold_fraction;
        return cut_data_around_peak(x_values, y_values, new_threshold, min_points_fraction);
    }

    // Extract the subset of data around the peak
    std::vector<double> x_cut;
    // I will initialize y_cut as double* instantly I think. What happens if I return the pointer to the first element in
    // the array? will I have memory issues because the array only lives in this function? No right? Should be fine if I return it?
    // for sure this works would be really weird if not
    double* y_cut = new double[above_threshold_indices.size()];     // we allocate the memory

    for (size_t i = 0; i < above_threshold_indices.size(); i++) {
        size_t ind = above_threshold_indices[i];        // index in the y_values and x values
        x_cut.push_back(x_values[ind]);
        // now we just set the i-th value
        y_cut[i] = y_values[ind];       // what is it doing here exactly? isnt y_cut[i] the same as y_cut + i? Seems not so, since you checked out pointers last time the syntax seems to have changed like always
    }

    // Add the offset back to y_cut
    for (size_t i = 0; i < x_cut.size(); ++i) {
        y_cut[i] += y_offset;
    }

    return std::make_pair(x_cut, y_cut);
}

double get_autocorrtime(double* f, int f_size, double ds) {
    // okay we could get the standard dev here already becasue we already computed it, but to be honest this wont be
    // computationally intensive in comparison to the rest that the function does
    // indeed for the integration we need the timestep (between the xi measurements)
    // Well the fing timestep is deviating in the cum equilibration observer...
    // maybe with the new error we dont need adaptive stepsizes anymore? because the error will be large for
    // slowly relaxing systems?
    // We write the function now for a constant dt
    // first we could calculate C(0) or the standarddev of xi (why are you talking about xi, it works for every observable
    // what do we do about the min ind? this should be part of this function, but copying a vector will also be expensive
    // tbh for this double* arrays would be insane?
    // we could construct a doulbe* vec_ptr = &vector[min_pos], that way we do not have to allocate more memory
    // cout << "f_size = " << f_size << endl;
    double avg_f = accumulate(f, f + f_size, 0.0) / (double)(f_size);
    // I think the average cannot be just the average of all used values, it might have to be the one in the interval until T - dist
    // cout << "avg_f = " << avg_f << endl;
    std::vector<double> diff(f_size);       // I mean we can use vectors here again if we want, then we dont have to deal with memory management
    std::transform(f, f + f_size, diff.begin(), [avg_f](double x) { return x - avg_f; });
    double variance_f = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0) / (double)f_size;
    variance_f = max(variance_f, 1e-5);
    // cout << "variance_f = " << variance_f << endl;
    // okay this is the variance
    // okay we have the simplest one of infinity values so we are basically done
    // the total duration of the simulation is, or of the sample that we are considering
    double T = f_size * ds;     // and we integrate from -T/2 to T/2
    vector<double> norm_autocorrelation_function{};     // I think we will use the normalized autocorrelation function so that we dont have to modify this thing if we want to calculate the autocorrelation time
    // we have discrete possible time distances
    for(int dist = 0; dist < f_size; dist++) {
        // dist is the time distance measured in timesteps ds, the actual time distance is
        double t = (double)dist * ds;
        // for this t we somehow need to calculate the autocorrelation function
        // I think you have to stop here because otherwise you wont be goin in an hour...
        // the easiest way to calculate the autocorrleation function will just be through its definition?
        // definition is C(t) = < f(s) f(s + t) > - <f>^2 ... <f> we already know
        // for < f(s) f(s+t) > we need to collect every f value pair that is seperated by the time t
        // I think for now we can think just of the discrete number of values in between
        // so we use differences from -fsize/ 2 to fsize/2.... Using the edge values is probably not useful since
        // there is only one value pair and this will be inflicted with strong statistical flucutations? This is probably where
        // the cutoff comes to play that this other book was talking about
        // so where should we cut? maybe for now just at fsize/4, later we can employ this self consistency thingy

        // for every dist we have the integration, but what are you talking about with f_size / 2?  You have C(t) and this
        // t should run from -fsize/2 to fsize/2
        // what we are doing now is integrate over the simulation time, which is this finding of the value pairs that
        // you were talking about
        // They should be easy to find because we just start at the first value that we have, its partner is 0 + dist,
        // then we go to 1 with the partner 1 + dist and at the end we have f_size - dist and f_size
        // we call the partners f_left and f_right or f_early and f_late
        // f_early is just the normal f pointer and f_late starts at f + dist
        // okay the thing is that the distance could strongly speaking also be negative, but does that matter?
        // and isnt just C(-t) = C(t) ?
        // ah damn I think I would have to write that down, for now we just use the absolute value of t
        double* f_early = f;
        double* f_late = f + dist;

        // the number of pairs we have is just
        size_t nr_integration_values = f_size - dist;
        // lets try it with the average until t - dist * ds
        double avg_early = accumulate(f_early, f_early + nr_integration_values, 0.0) / (double)nr_integration_values;
        double avg_late =  accumulate(f_late, f_late + nr_integration_values, 0.0) / (double)nr_integration_values;
        // Now we fill up a vector with the multiplied values
        vector<double> prod(nr_integration_values);
/*        cout << "f_early:" << endl;
        print_array(f_early, nr_integration_values);
        cout << "f_late:" << endl;
        print_array(f_late, nr_integration_values);*/
        // cout << "dist = " << dist << endl;
        // cout << "avg_early = " << avg_early << endl;
        std::transform(f_early, f_early + nr_integration_values, f_late, prod.begin(), std::multiplies<double>());
        // cout << "prod" << endl;
        // print_vector(prod);
        // I hope this is correct syntax, the aim is to have the multiplied values f(t) * f(t +dist) inside the prod vector
        // the prod vector can then be reduced and multiplied with the time distance between the measured values (ds)
        // Afterwards we have to average the integral with 1 / T, but T is just ns * ds so we can in this case just average with ns which is f_size
        double autocorr_value = accumulate(prod.begin(), prod.end(), 0.0) / (double)nr_integration_values -  avg_f * avg_early - avg_f * avg_late + avg_f * avg_f;
        if(autocorr_value < 0) {
            break;      // I am not sure if this is legal but who cares...
        }
        //cout << "autocorr_value " << autocorr_value << endl;
        // I think we should save the autocorrelation values somewhere
        norm_autocorrelation_function.push_back(autocorr_value / variance_f);
        //cout << "norm_autocorrelation value = " << autocorr_value / variance_f << endl;
    }
    // print_vector(norm_autocorrelation_function);
    // now we have all C(t) values, then we integrate them to obtain the autocorrelation time
    // this time the actual time difference is relevant
    double autocorrelation_time = accumulate(norm_autocorrelation_function.begin(), norm_autocorrelation_function.end(), 0.0) * ds;
    // cout << "autocorr time = " << autocorrelation_time << endl;
    // cout << "ds = " << ds << endl;
    // exit(0);
    // and this should be it. I think through the fact that we always know directly where the partner is and dont have many loops the
    // computation time will not be small but also not huge, probably okay if it runs every 1000 steps or so
    return autocorrelation_time;
}

struct SquaredDifferenceFromMean
{
    float mean;

    SquaredDifferenceFromMean(float _mean) : mean(_mean) {}

    __host__ __device__
    float operator()(float x) const
    {
        float diff = x - mean;
        return diff * diff;
    }
};

struct ConjugateAndMultiply
{
    ConjugateAndMultiply() {}

    template <class valuetype>
    __host__ __device__
    valuetype operator()(valuetype x) const
    {
        return valuetype(x.x * x.x + x.y * x.y, 0.0);      // ?????
    }
};

struct SubtractMeanDivideSTD
{
    double mean, var;
    SubtractMeanDivideSTD(double mean, double var): var(var), mean(mean) {}

    __host__ __device__
    cufftDoubleComplex operator()(double x) const
    {
        return make_cuDoubleComplex((x - mean) / sqrt(var), 0.0);
    }
};


struct normalize_fft_autocorr : thrust::unary_function<thrust::tuple<cufftDoubleComplex, size_t>, double> {
int f_size;
normalize_fft_autocorr(int f_size) : f_size(f_size), thrust::unary_function<thrust::tuple<cufftDoubleComplex, size_t>, double>(){
}
template<class Tup>
__host__ __device__ double operator()(Tup tup) const {
    // first we need to know in which system we are
    double autocorr_value = thrust::get<0>(tup).x;
    int index = thrust::get<1>(tup);
    return  autocorr_value / (2.0 * f_size * (f_size - index));
}
};


double get_autocorrtime_gpu(double* f, int f_size, double ds) {
    // create a thrust device vector of f
    thrust::device_vector<double> device_f(f_size);
    thrust::copy(f, f + f_size, device_f.begin());
/*    for(int i = 0; i < 20; i++) {
        cout << f[i] << endl;
        cout << device_f[i] << endl;
    }*/
    double avg_f = thrust::reduce(device_f.begin(), device_f.end(), (double)0.0) / (double)(f_size);
    cout << "avg_f = " << avg_f << endl;
    // what is supposed to be the problem here? it should be fine but why is the IDE complaining at transform reduce it didnt before...
    double variance_f = thrust::transform_reduce(device_f.begin(), device_f.end(), SquaredDifferenceFromMean(avg_f), 0.0, thrust::plus<double>()) / (double)f_size;
    variance_f = max(variance_f, 1e-5);
    cout << "variance_f = " << variance_f << endl;
    double T = f_size * ds;     // and we integrate from -T/2 to T/2
    thrust::device_vector<double> norm_autocorrelation_function(f_size, 0.0);     // I think we will use the normalized autocorrelation function so that we dont have to modify this thing if we want to calculate the autocorrelation time
    for(int dist = 0; dist < f_size; dist++) {
        double t = (double)dist * ds;

        size_t nr_integration_values = f_size - dist;

        double avg_early = thrust::reduce(device_f.begin(), device_f.begin() + nr_integration_values, 0.0) / (double)nr_integration_values;
        double avg_late =  thrust::reduce(device_f.begin() + dist, device_f.end(), 0.0) / (double)nr_integration_values;

        //cout << "avg_early = " << avg_early << endl;
        //cout << "avg_late = " << avg_late << endl;

        thrust::device_vector<double> prod(nr_integration_values);
        thrust::transform(device_f.begin(), device_f.begin() + nr_integration_values, device_f.begin() + dist, prod.begin(), thrust::multiplies<double>());
        double autocorr_value = thrust::reduce(prod.begin(), prod.end(), (double)0.0) / (double)nr_integration_values -  avg_f * avg_early - avg_f * avg_late + avg_f * avg_f;

        //cout << "autocorr_value = " << autocorr_value << endl;
        //exit(0);
        if(autocorr_value < 0) {
            break;      // I am not sure if this is legal but who cares...
        }
        norm_autocorrelation_function[dist] = (autocorr_value / variance_f);    // is this asign fine? yes right?
    }
    double autocorrelation_time = thrust::reduce(norm_autocorrelation_function.begin(), norm_autocorrelation_function.end(), 0.0) * ds;
    return autocorrelation_time;
}

double get_autocorrtime_fft(double* f, int f_size, double ds) {
    // create a thrust device vector of f
    thrust::device_vector<double> device_f(f_size);
    thrust::copy(f, f + f_size, device_f.begin());
    // cout << "f_size = " << f_size << endl;
    double avg_f = thrust::reduce(device_f.begin(), device_f.end(), (double)0.0) / (double)(f_size);
    // cout << "avg_f = " << avg_f << endl;
    // what is supposed to be the problem here? it should be fine but why is the IDE complaining at transform reduce it didnt before...
    double variance_f = thrust::transform_reduce(device_f.begin(), device_f.end(), SquaredDifferenceFromMean(avg_f), 0.0, thrust::plus<double>()) / (double)f_size;
    variance_f = max(1e-5, variance_f);
    // cout << "variance_f = " << variance_f << endl;

    // it is just a large 1D FFT
    int fft_size = 2 * f_size + 1;

    int nr_batches = 1;
    // Do it directly with cufft? We only need one but I guess it will be also fast if we only do one large fft?
    cufftHandle plan;
    cufftCreate(&plan);
    cufftPlan1d(&plan, fft_size, CUFFT_Z2Z, nr_batches);

    // In and output? For even more performance we could indeed use an R2C transformation
    thrust::device_vector<cufftDoubleComplex> input_vector(fft_size);     // ... is this going to work? Does something like this work without cuda?
    thrust::device_vector<cufftDoubleComplex> output_vector(fft_size);

    // seems we need to use transform to fill the device vector
    thrust::transform(device_f.begin(), device_f.end(), input_vector.begin(), SubtractMeanDivideSTD(avg_f, variance_f));


    cufftExecZ2Z(plan, (cufftDoubleComplex*)thrust::raw_pointer_cast(input_vector.data()),
                 (cufftDoubleComplex*)thrust::raw_pointer_cast(output_vector.data()), CUFFT_FORWARD);
    cudaDeviceSynchronize();


    // now we somehow need to do complex multiplication
    thrust::transform(output_vector.begin(), output_vector.end(), input_vector.begin(), ConjugateAndMultiply());

    // And the inverse FFT?
    cufftHandle inv_plan;
    cufftCreate(&inv_plan);
    cufftPlan1d(&inv_plan, fft_size, CUFFT_Z2Z, nr_batches);

    // ah we probably need a vector that can use the real outpout
    cufftExecZ2Z(plan, (cufftDoubleComplex*)thrust::raw_pointer_cast(input_vector.data()),
                 (cufftDoubleComplex*)thrust::raw_pointer_cast(output_vector.data()), CUFFT_INVERSE);
    cudaDeviceSynchronize();

    thrust::device_vector<double> normalized_autocorr_function(f_size);
    auto output_index = thrust::make_zip_iterator(thrust::make_tuple(output_vector.begin(), thrust::counting_iterator<size_t>(0)));
    thrust::transform(output_index, output_index + f_size, normalized_autocorr_function.begin(), normalize_fft_autocorr(f_size));

    cufftDestroy(plan);
    cufftDestroy(inv_plan);

    // we want to add the automatic windowing algorithm. I think an exclusive scan that you convieniently found earlier could work a charm
    // no, inclusive scan will be even better. The normalized autocorr function is already double which is nice, so we just need another

    thrust::device_vector<double> autocorr_time_scan(f_size);
    // TODO I might be lacking a factor of two here
    thrust::inclusive_scan(normalized_autocorr_function.begin(), normalized_autocorr_function.end(), autocorr_time_scan.begin());

    // now we just have to select the smallest index/time that satisfies that M >= 5 * autocorr_time
    // both M and tau are not multiplied with ds yet so we can directly compare the index of the function with its value
    // I dont see much benefit from doing this on gpu so we try a simple implementation
    thrust::host_vector<double> autocorr_time_scan_host(autocorr_time_scan);
    int M_select = 0;
    double max_autocorrelation_time = 0;
    for(int M = 1; M < f_size; M++) {
        if( M >= 5 * autocorr_time_scan_host[M]) {
            M_select = M;
            break;
        } else if(autocorr_time_scan_host[M] > max_autocorrelation_time) {
            max_autocorrelation_time = autocorr_time_scan_host[M];
            M_select = M;
        }
        // If it finds none.. what then? just use f_size? Use the M that yields the largest autocorr_time_scan?
    }
    double autocorrelation_time = autocorr_time_scan_host[M_select] * ds;
    autocorrelation_time = max(autocorrelation_time, ds);   // It should not be smaller than the stepsize and especially not negative
    // double autocorrelation_time = thrust::transform_reduce(normalized_autocorr_function.begin(), normalized_autocorr_function.end(), OnlyPositive(), 0.0, thrust::plus<double>()) * ds;
    return autocorrelation_time;
}


void append_parameter(fs::path filepath, string value_name, double value) {
    path txt_path = filepath.replace_extension(".txt");
    ofstream parameter_file;
    // can there only be one output stream?
    parameter_file.open(txt_path, ios::app);
    parameter_file << value_name << "," << value << endl;
    parameter_file.close();
}

#endif //CUDAPROJECT_MAIN_CUDA_CUH
