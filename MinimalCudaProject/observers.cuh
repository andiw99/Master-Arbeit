//
// Created by andi on 15.08.23.
//

#ifndef CUDAPROJECT_OBSERVERS_CUH
#define CUDAPROJECT_OBSERVERS_CUH
#include <boost/asio/ip/host_name.hpp>


void write_parameters(ofstream& file, map<Parameter, double> paras) {
    cout << "Using the new write parameters" << endl;
    for(const auto& pair : paras) {
        file << parameter_names[pair.first] << "," << pair.second << endl;
    }
}
// Okay because of our fancy vector<observer*> stuff in the simulation, we now need new templateable observer classes
template <class system, class State>
class obsver {
public:
    ofstream ofile;

    virtual void operator()(system& sys, const State &x , double t ) {
        cout << "Base Observer withoud const gets called" << endl;
    }
    obsver() {
    }
    void write(system &sys, const State &x , double t ) {

    }

    void open_stream(fs::path filepath) {
        cout << "open stream is called with filepaht: " << filepath <<  endl;
        create_dir(filepath.parent_path());
        ofile.open(filepath);
    }

    void open_app_stream(fs::path filepath) {
        cout << "open append stream is called with filepaht: " << filepath <<  endl;
        // create_dir(filepath.parent_path()); If it is append we dont need that line anymore since by definition the folder has to exist
        ofile.open(filepath, ios::app);
    }

    void close_stream() {
        cout << "closing stream of " << this->get_name() << endl;
        ofile.close();
    }
    virtual void init(fs::path folderpath, map<Parameter, double>& paras, const system &sys) {
        cout << "dummy init is called" << endl;
    }

    virtual string get_name() {
        return "base observer";
    }

    ~obsver() {
        cout << "deleting " << this->get_name() << endl;
        close_stream();
    }

    string construct_filename(int run_nr) {
        string job_id_str;
        if (const char* job_id = std::getenv("SLURM_JOB_ID")){
            job_id_str = string(job_id);
        } else {
            job_id_str = "local";
        }
        return to_string(run_nr) + "-" + get_current_time().substr(0, 4) + "-" + job_id_str + "-" + boost::asio::ip::host_name();
    }


};

template <class system, class State>
class standard_observer : public obsver<system, State> {
public:
    double timepoint = 0;
    // Okay the observing pattern could be totally wild, i probably somehow have to initialize the observer
    // outside of the simulation class We definetely need its own constructor here
    int nr_values;
    double write_interval = 1;
    typedef obsver<system, State> obsver;
    using obsver::ofile;
    using obsver::open_stream;
    using obsver::open_app_stream;
    using obsver::close_stream;
    standard_observer(int nr_values): nr_values(nr_values) {
        // this one has just the logic that we write nr_values equidistant points, applies to quench aswell as to
        // Relaxation. The init is a bit different since we need to find out the write interval
        // It also writes the run parameters to a file
    }
    virtual void init(fs::path path, map<Parameter, double>& paras, const system &sys)  {
        // I think this will have less performance impact than an if statement catching the first observer operation
        // Make sure to use this observer only with systems that have a get_end_T method

        // open the file to write the info to, in this case it will be just run_nr.csv
        timepoint = 0.0;
        // I think we will add the run number to the paras of the run
        if(boost::asio::ip::host_name() == "andi-B550") {
            cout << "CHANGING RANDOM TO NORMAL RANDOM" << endl;
            paras[Parameter::curand_random] = 0;
        }

        int run_nr = (int)paras[Parameter::run_nr];
        bool pick_up = (paras[Parameter::random_init] == -1.0);     // if pickup is true, we have to do things a bit differently
        fs::path folderpath;
        if(pick_up) {
            // first of all the path of the folder is the parentpath of the current folderpath
            folderpath = path.parent_path();
        } else {
            folderpath = path;
        }
        // we just write the parameters first
        close_stream();
        string filename = obsver::construct_filename(run_nr);
        open_stream(folderpath / (filename + ".txt"));
        write_parameters(ofile, paras);
        ofile << "system" << "," << sys.get_name() << endl;
        // dont forget to close!;
        close_stream();
        if (pick_up) {
            path += ".csv";
            open_app_stream(path);  // with plus since we have the filename stem here?
            write_interval = (paras[end_time] - paras[start_time]) / (double)nr_values;
            cout << "standard observer write_interval " << write_interval;
            timepoint = paras[Parameter::start_time]; // + write_interval;      // We dont want to write the starting position again
            cout << "pickup, setting the next write to " << timepoint << endl;
        } else {
            open_stream(folderpath / (filename + ".csv"));
        }
    }
    void operator()(system &sys, const State &x , double t ) override {
        if(t > timepoint) {
            // write
            write_state(sys, x, t);
            // advance timeoint
            timepoint += write_interval;
        }
    }


    void write_state(const system &sys, const State &x, double t) {
        cout << "t = " << t << endl;
        // Doing some refactoring, not putting t in every row and Temperature directly after T
        double T = sys.get_cur_T();
        size_t dim_size_x = sys.get_dim_size_x();
        size_t dim_size_y = sys.get_dim_size_y();

        ofile << t << "," << T << ",";
        for(int i = 0; i < dim_size_x * dim_size_y; i++) {
            ofile << x[i] << ",";
            // cout << x[i] << endl;

        }
        // Zeilenumbruch
        ofile << "\n";
        cout << "wrote file" << endl;
    }
    string get_name() override {
        return "standard observer";
    }


};

template <class system, class State>
class relax_observer : public standard_observer<system, State> {
    // Okay the observing pattern could be totally wild, i probably somehow have to initialize the observer
    // outside of the simulation class We definetely need its own constructor here
public:
    typedef standard_observer<system, State> standard_observer;
    using standard_observer::write_interval;
    using standard_observer::nr_values;
    using standard_observer::operator();
    relax_observer(int nr_values): standard_observer(nr_values) {
        // should this be the observer that only writes a few values for the equilibrium process?
        // I am not sure at the moment, I think i am just going to write the simplest one for now, so
        // just equidistant saving values.
        // We try to initialize the observer outside the class since it probably needs
        // multiple unique-to-observer parameters
        // problem is it still might depend on system parameters. and i don't want another if for every step
        // even though i am not really sure if it impacts the performance AT ALL or if the compiler knows better
    }
    void init(fs::path folderpath, map<Parameter, double>& paras, const system &sys) override {
        // I think this will have less performance impact than an if statement catching the first observer operation
        // Make sure to use this observer only with systems that have a get_end_T method
        double end_t = paras[end_time];
        write_interval = end_t / (double)nr_values;
        cout << "relax obs init is called, write_interval is  " << write_interval << endl;
        standard_observer::init(folderpath, paras, sys);
    }

    string get_name() override {
        return "relax observer";
    }
};

template <class system, class State>
class quench_observer : public standard_observer<system, State> {
    // Okay the observing pattern could be totally wild, i probably somehow have to initialize the observer
    // outside of the simulation class We definetely need its own constructor here
public:
    typedef standard_observer<system, State> standard_observer;
    using standard_observer::write_interval;
    using standard_observer::nr_values;
    quench_observer(int nr_values): standard_observer(nr_values) {
        // should this be the observer that only writes a few values for the equilibrium process?
        // I am not sure at the moment, I think i am just going to write the simplest one for now, so
        // just equidistant saving values.
        // We try to initialize the observer outside the class since it probably needs
        // multiple unique-to-observer parameters
        // problem is it still might depend on system parameters. and i don't want another if for every step
        // even though i am not really sure if it impacts the performance AT ALL or if the compiler knows better
    }
    void init(fs::path folderpath, map<Parameter, double>& paras, const system &sys) override {
        // I think this will have less performance impact than an if statement catching the first observer operation
        // Make sure to use this observer only with systems that have a get_end_T method
        cout << "Quench obs init is called" << endl;
        double end_T = sys.get_end_t();
        cout << "sys.get_end_t():  " << sys.get_end_t() << endl;
        write_interval = end_T / (double)nr_values;
        cout << "Write interval is " << write_interval << endl;
        standard_observer::init(folderpath, paras, sys);
    }

    string get_name() override {
        return "quench observer";
    }
};

template <class system, class State>
class equilibration_observer : public standard_observer<system, State> {
    // Okay the observing pattern could be totally wild, i probably somehow have to initialize the observer
    // outside of the simulation class We definetely need its own constructor here
public:
    bool startpoint = true;         // If the startpoint is true we will write down on call
    typedef standard_observer<system, State> standard_observer;
    equilibration_observer(): standard_observer(0) {        // we just use 0 to initialize the number of points. Do we even need a constructor?
    }
    void init(fs::path folderpath, map<Parameter, double>& paras, const system &sys) override {
        // This observer is supposed to just write down once the system is equilibrated and once the system is started
        // The question is how we tell the observer that we are equilibrated. We should probably not think to much about it
        // and just check every step if the system is equilibrated. This means we do not even need a write interval
        // But we need to tell the observer somehow if it is the first iteration
        startpoint = true;      // always needs to be set true if we run a second system
        // cout << "equilibration obs init is called, write_interval is  " << endl;
        standard_observer::init(folderpath, paras, sys);
    }

    void operator()(system &sys, const State &x , double t) {
        if(startpoint) {
            // I just hope that this is optimized by the compiler...
            this->write_state(sys, x, t);
            // set startpoint to be false
            startpoint = false;
        } else if(sys.is_equilibrated()) {
            this->write_state(sys, x, t);
        }
    }

    string get_name() override {
        return "equilibration observer";
    }
};

template <class system, class State>
class m_observer : public obsver<system, State>{
    // Okay the observing pattern could be totally wild, i probably somehow have to initialize the observer
    // outside of the simulation class We definetely need its own constructor here
    typedef obsver<system, State> obsver;
    using obsver::ofile;
    using obsver::open_stream;
    using obsver::close_stream;
    int nr_values;
    double write_interval = 1;
    double timepoint = 0;
    double dt;
public:
    m_observer(int nr_values) : nr_values(nr_values) {
    }
    void init(fs::path folderpath, map<Parameter, double>& paras, const system &sys) override {
        int run_nr = (int)paras[Parameter::run_nr];
        timepoint = 0.0;
        dt = paras[Parameter::dt];

        close_stream();
        open_stream(folderpath / (obsver::construct_filename(run_nr) + ".mag"));
        cout << this->get_name() << " init called" << endl;
        ofile << "t,m_vec" << endl;

        // I think this will have less performance impact than an if statement catching the first observer operation
        // Make sure to use this observer only with systems that have a get_end_T method
        double end_t = paras[end_time];
        write_interval = end_t / (double)nr_values;
    }

    string get_name() override {
        return "m observer";
    }

    void operator()(system &sys, const State &x , double t ) override {
        if(t > timepoint - dt / 2.0) {
            thrust::host_vector<double> m_val_vec = sys.calc_m_vec(x);
            // with this observer we could think about writing at the end of the simulation, or at least rewriting
            // the file if we changed the stepsize
            // How do we do this, we cannot save every m value... I think we will try for now?
            ofile << t << ";";
            for (double m_val : m_val_vec) {
                ofile << m_val << ",";
            }
            ofile << endl;
            timepoint += write_interval;
        }
    }

};

template <class system, class State>
class m_equilibration_observer_adaptive: public obsver<system, State>{
    // does the same stuff as the normal equilibration observer but adapts the write density
    typedef obsver<system, State> obsver;
    using obsver::ofile;
    using obsver::open_stream;
    using obsver::open_app_stream;
    using obsver::close_stream;
    double write_interval = 1;
    double timepoint = 0;
    vector<double> m_vec{};
    vector<double> times{};
    int min_m_nr = 500;
    int m_nr_gpu = 2000;                 // nr of cum values at which the calculation switches to gpu
    double write_density = 1.0 / 100.0;      // we increase the write density to get a more useful value for the autocorrelation function?
    double min_write_density = 1.0 / 10000.0;   // I think 500 was already a bit slow for my taste
    double dt = 0.01;
    double dt_half;
    double equil_cutoff = 0.1;              // since the equilibration might influce the mean of m_vec a lot we cut a certain portion of m_vec values
    double max_error= 0.001;
    int m_nr = 0;                         // current number in the averageing process
    bool equilibrated = false;                      // for the usecase of the quench with dynamic equilibration
    double eval_factor = 10000.0;
    fs::path filepath;
    double max_moving_factor = 0.005;

public:
    m_equilibration_observer_adaptive(int min_m_nr) : min_m_nr(min_m_nr) {}

    m_equilibration_observer_adaptive(int min_m_nr, double write_density) : min_m_nr(min_m_nr), write_density(write_density) {}

    m_equilibration_observer_adaptive(int min_cum_nr, double write_density, double equil_cutoff):
            min_m_nr(min_cum_nr), write_density(write_density), equil_cutoff(equil_cutoff) {}

    void init(fs::path path, map<Parameter, double>& paras, const system &sys) override {
        cout << "initializing m equilibration observer" << endl;
        int run_nr = (int)paras[Parameter::run_nr];
        max_error = paras[Parameter::equil_error];
        if(paras[Parameter::moving_factor] != 0.0) {
            // check that we do not set the moving factor to be zero
            max_moving_factor = paras[Parameter::moving_factor];
        }
        if (paras[Parameter::min_cum_nr]){
            min_m_nr = (int)paras[Parameter::min_cum_nr];
        }
        timepoint = paras[Parameter::start_time];
        dt = paras[Parameter::dt];
        equilibrated = false;
        // we also need to reset m_vec and times, dont we?
        m_vec = vector<double>{};
        times = vector<double>{};
        close_stream();

        bool pick_up = (paras[Parameter::random_init]  == -1.0);
        if(!pick_up) {
            filepath = path / (obsver::construct_filename(run_nr) + ".mag");
            open_stream(filepath);
            ofile << "t;m" << endl;
        } else {
            path += ".mag";
            filepath = path;
            readMagFromFile(path, m_vec, times);
            write_interval = times[1] - times[0];
            write_density = dt / write_interval;    // We should adapt the write interval of the file that we read
            cout << "picked up write density of " << write_density << "  corresponding to a write interval of " << write_interval;
            // the other observers that dont read the values in actually have a problem with defining the timepoint?
            // problem is that you set the timepoint here to the back of m_vec, but it might be that the .csv
            // file has not written that many states because it somehow didnt finish...
            timepoint += write_interval;
            cout << "Pickup, setting timepoint of " << this->get_name()  << " to " << timepoint << endl;
            open_app_stream(path);
        }
        cout << this->get_name() << " init called" << endl;


        // I think the starting write interval should be every one hundred steps
        dt_half = dt / 2.0;
        write_interval = dt / write_density;
        cout << "I am sure the write interval is " << write_interval << endl;
    }

    string get_name() override {
        return "m equilibration observer";
    }

    void operator()(system &sys, const State &x , double t) override {
        if(t > timepoint - dt_half) {
            // advancing the cum nr
            m_nr++;
            // we calculate the cumulant and write it down
            thrust::host_vector<double> m_val_vec = sys.calc_m_vec(x);
            // with this observer we could think about writing at the end of the simulation, or at least rewriting
            // the file if we changed the stepsize
            // How do we do this, we cannot save every m value... I think we will try for now?
            ofile << t << ";";
            for (double m_val : m_val_vec) {
                 ofile << m_val << ",";
            }
            ofile << endl;
            // add the cumulant and the times to the vectors to keep track
            // m_vec.push_back(m_val);
            m_vec.insert(m_vec.end(), m_val_vec.begin(), m_val_vec.end());
            times.push_back(t);
            // now we want to see if we have to adjust the write interval
            // we have to make sure that we got 5 fresh cum values
            if(!equilibrated) {
                int nr_m_measurements = m_vec.size() / m_val_vec.size();      // check how often we measured
                int nr_m_values = m_vec.size();
                // cout << "nr_m_values: " << nr_m_values << endl;
                // now use the last avg_nr of cum values to calculate a mean m_vec
                // use transform reduce?
                // does it work like this? m_vec.end() - avg_nr is the n-th last value in the vector?
                if(nr_m_measurements >= min_m_nr){
                    int eval_every_n_steps = (int)(eval_factor * write_density) + 1;   // I think this + 1 makes everything nicer? Or does it?
                    if(m_vec.size() % eval_every_n_steps == 0) {
                        // We dont need to calculate the error and stuff everytime we write down a m_vec
                        // IF the density is small, we should write down more often
                        // the equilibration phase might influence the mean a lot, should we cut of the first x% of the values?
                        int min_ind = (int)(equil_cutoff * nr_m_values);
                        int nr_values_to_use = nr_m_values - min_ind;
                        // we calculate the Binder cumulant according to the ergodic hyptheses from all the ms we extracted
                        // TODO if you have time you can do this on GPU but it shouldnt be to useful
                        double* m = &m_vec[min_ind];
                        // cout << "min ind " << min_ind << endl;
                        // cout << "nr_values_to_use " << nr_values_to_use << endl;
                        double m_L2 = std::transform_reduce(m, m + nr_values_to_use,
                                                            0.0, // initial value for the reduction (sum)
                                                            std::plus<double>(), 
                                                            [](double m_val) -> double { return m_val * m_val; });
                        m_L2 /= (double) nr_values_to_use;
                        // TODO strictly speaking we would have to calculate autocorrelation times for m² and m⁴ but I
                        // want to get that running first
                        double m_L2_err = std::transform_reduce(m, m + nr_values_to_use,
                                                                0.0, // initial value for the reduction (sum)
                                                                std::plus<double>(), // transformation (square)
                                                                [&m_L2](double m_val) { return pow((m_val * m_val - m_L2), 2) ; });
                        m_L2_err /= (double)pow(nr_values_to_use, 2);
                        m_L2_err = sqrt(m_L2_err);
                        double m_L4 = std::transform_reduce(m, m + nr_values_to_use,
                                                            0.0, // initial value for the reduction (sum)
                                                            std::plus<>(), // transformation (square)
                                                            [](double m_val) { return (pow(m_val, 4)); });
                        m_L4 /= nr_values_to_use;
                        double m_L4_err = std::transform_reduce(m, m + nr_values_to_use,
                                                                0.0, // initial value for the reduction (sum)
                                                                std::plus<>(), // transformation (square)
                                                                [&m_L4](double m_val) { return pow(pow(m_val, 4) - m_L4, 2); });
                        m_L4_err /= (double)pow(nr_values_to_use, 2);
                        m_L4_err = sqrt(m_L4_err);
                        double U_L = m_L4 / (m_L2 * m_L2);
                        double U_L_error = sqrt(pow(1 / m_L2 / m_L2 * m_L4_err, 2) + pow(2 * m_L4 / pow(m_L2, 3) * m_L2_err, 2));

                        double U_L_variance = U_L_error * U_L_error;
                        double rel_stddev_total = U_L_error / U_L;

/*
                        std::vector<double> diff_total(nr_m_values - min_ind);

                        std::transform(m_vec.begin() + min_ind, m_vec.end(), diff_total.begin(), [U_L](double x) { return x - U_L; });

                        double dist_var = std::inner_product(diff_total.begin(), diff_total.end(), diff_total.begin(), 0.0) / (nr_m_values - min_ind);
                        // this is the variance of the distribution (in wrong because the values are correlated), so we
                        // need the autocorrelation time to adjust for this. The function works with arrays
                        double* U_L_arr = &m_vec[min_ind];
                        double autocorr_time;

                        if(nr_m_values - min_ind > m_nr_gpu) {
                            autocorr_time = get_autocorrtime_fft(U_L_arr, nr_m_values - min_ind, write_interval);
                        } else {
                            autocorr_time = get_autocorrtime(U_L_arr, nr_m_values - min_ind, write_interval);  // actually ds is just the write interval? which should be 1 or something like this
                        }
                        double U_L_variance = 2 * autocorr_time / ((nr_m_values - min_ind) * write_interval) * dist_var;
                        double rel_stddev_total = sqrt(U_L_variance) / U_L;

                        cout << "Distriubtion variance: " << dist_var << endl;
                        cout << "autocorrelation time: " << autocorr_time << endl;
                        */
                        double avg_m = std::reduce(m, m + nr_values_to_use) / (double)nr_values_to_use;
                        cout << "avg_m = " << avg_m << endl;
                        cout << "m_vec variance: " << U_L_variance << endl;
                        cout << "U_L: " << U_L << endl;
                        cout << "rel_stddev_total = " << rel_stddev_total << endl;


                        // Okay we found out that this is not the real error, so we need a function that calculates the error
                        // or first a function that calculates the integrated autocorrleation time
                        // We indded have two auto correlation times for the two directions. Could they be individually sized?
                        // the question is now if we want to extract U_L and its error if we are equilibrated
                        // we definitely should? But how and where? somehow into the parameter file?
                        if(rel_stddev_total < max_error) {
                            // I think in this case I want to see wether there is still a significant upwards or
                            // downards trend to ensure the system is equilibrated
                            // the question is what recent means, the last 20% of values? It shouldnt be to much
                            // that it wouldnt be recent anymore and it should be so few that statistical fluctuations
                            // would average out
                            // to be honest, the trend stuff is probably the prestuff we should do before even considering
                            // to calculate the error? Error calculation after all is pretty expensive
                            double moving_factor = getMovingFactor(nr_m_values, avg_m, min_ind);

                            cout << "U_L = " << U_L << " +- " << rel_stddev_total * U_L << endl;
                            cout << "MOVING FACTOR = " << moving_factor << endl;

                            if(moving_factor < max_moving_factor) {
                                // we set the system to be equilibrated
                                cout << "The system equilibrated, the equilibration lastet to t = " << t << endl;
                                sys.set_equilibration(t);
                                // write cumulant average
                                append_parameter(filepath, "U_L", U_L);
                                append_parameter(filepath, "m", avg_m);
                                append_parameter(filepath, "m2", m_L2);
                                append_parameter(filepath, "m4", m_L4);

                                // write error
                                append_parameter(filepath, "U_L_error", rel_stddev_total);
                                // write autocorrelation time
                                // append_parameter(filepath, "autocorrelation_time_U_L", autocorr_time);
                                // Also write the moving factor
                                append_parameter(filepath, "moving_factor_U_L", moving_factor);
                                // once we did this, we dont want to do that a second time?
                                equilibrated = true;
                            }
                        }
                        if(!equilibrated) { // TODO
                            // adapt_write_interval(nr_m_values, autocorr_time);
                        }
                    }
                }
            }
            timepoint += write_interval;
        }
    }

    double getMovingFactorOld(int nr_cum_values, long double avg_U_L) {
        int recent_ind = (int) (nr_cum_values * 0.9);
        double* recent_U_L = &m_vec[recent_ind];
        // the idea is to check if the differences between the values show a trend to be positive
        // or negative. But I think that comes done to calculating the deviation of the recent
        // average from the total average
        // small differences will
        // Idea would be to compare this deviation to the deviation per step
        // If we are in a high temperature state, the deviation per step will be relatively larger
        // If should not be dependent on the stepsize
        // first of all compute the avergage absolute difference per step
        int recent_size = nr_cum_values - recent_ind;
        double recent_avg_U_L = reduce(recent_U_L, recent_U_L + recent_size) / (double)recent_size;
        double mean_abs_delta_U = meanAbsDifference(recent_U_L, recent_size);
        // and we need the recent mean
        double moving_factor = fabs(avg_U_L - recent_avg_U_L) / ((double)recent_size * mean_abs_delta_U);

        cout << "recent_avg_U_L = " << recent_avg_U_L << endl;
        cout << "|U_L_recent - m_vec| = " << avg_U_L - recent_avg_U_L << endl;
        cout << "mean abs delta = " << mean_abs_delta_U << endl;
        return moving_factor;
    }

    double getMovingFactor(int nr_values, long double avg, int min_ind) {
        int recent_ind = (int) (nr_values * 0.8);
        // the idea is to check if the differences between the values show a trend to be positive
        // or negative. But I think that comes done to calculating the deviation of the recent
        // average from the total average
        // small differences will
        // Idea would be to compare this deviation to the deviation per step
        // If we are in a high temperature state, the deviation per step will be relatively larger
        // If should not be dependent on the stepsize
        // first of all compute the avergage absolute difference per step
        int recent_size = nr_values - recent_ind;
        double var_start = reduce(m_vec.begin() + min_ind, m_vec.begin() + recent_ind) / (double)(recent_ind - min_ind);
        double var_end = avg;

        double* recent_U_L = &m_vec[recent_ind];      // we need it still for the mean delta?
        double mean_abs_delta_U = meanAbsDifference(recent_U_L, recent_size);
        // and we need the recent mean
        double moving_factor = fabs(var_end - var_start) / ((double)recent_size * mean_abs_delta_U);

        cout << "var_start = " << var_start << endl;
        cout << "|var_end - var_start| = " << avg - var_start << endl;
        cout << "mean abs delta = " << mean_abs_delta_U << endl;
        return moving_factor;
    }

    void adapt_write_interval(int nr_cum_values, double autocorr_time, int values_per_autocorr_time=20) {
        // we want to adapt the write density on the autocorrelation time. But everytime we do this
        // we have to delete the values that have smaller spacing between them (the intermediate values)
        // It should be some kind of switch statement or something like this?
        // We want to write approximately 10 cumulant values during one autocorrelation time? Or maybe 100 ?
        // I think 10 should be fine
        double autocorr_time_floored = pow(10.0, floor(log10(autocorr_time)));      // we floor the autocorrelation time to the next power of 10, so 3000 would get rounded to 1000 and 101 to 100
        // we do this as we do not want to continuously change our stepping but in discrete steps
        double new_write_interval = autocorr_time_floored / (double)values_per_autocorr_time;    // if autocorr is 1000, the write interval is 100, using 100 writes per autocorr time now because ten is not satisfying somehow
        if(new_write_interval > write_interval) {
            // in this case we have to drop values accordingly to new write interval
            int keep_every_nth_val = (int)(new_write_interval / write_interval); // they should yield a glatt integer since we only work with powers of ten write interval might have been 1 and now it is 10 so we keep every 10th value
            int new_nr_cum_values = nr_cum_values / keep_every_nth_val;         // if we have 101 values with a spacing of 1, we want to keep 101, 91, ..., 1, making up for 11 values
            if (nr_cum_values % keep_every_nth_val != 0) {
                // if this is the case we can extract one more value
                new_nr_cum_values++;
            }
/*            cout << "the old number of cum values is " << nr_cum_values << endl;
            cout << "the new number of cum values is " << new_nr_cum_values << endl;*/
            // if we have 100 values with a spacing of 1, those would be the indices 0 to 99. We want to keep 99 and 89 (distance of 10 ds), 79, ..., 9 making up 10 values
            // if we have 101 values with a spacing of 1, those would have the indices 0 to 100. We would want to keep 100, 90, ..., 10, 0 so eleven values? Right now we would be missing out on 0, which would not be two bad...
            vector<double> new_U_L(new_nr_cum_values);
            vector<double> new_times(new_nr_cum_values);
            for(int i = 0; i < new_nr_cum_values; i++) {
                // TODO is the plus one correct? If you get memory errors look here?
                // okay so we work from behind because instead of keeping 0 - 100 we want to keep 1 - 101. But the order has to stay the same for the correlation function
                // cout << "i = " << i << "   nr_cum_values - keep_every_nth_val * i - 1 = " << nr_cum_values - keep_every_nth_val * i - 1 << "    m_vec = " << m_vec[nr_cum_values - keep_every_nth_val * i - 1] << endl;
                new_U_L[new_nr_cum_values - i - 1] = m_vec[nr_cum_values - keep_every_nth_val * i - 1];       // - one because in a vector of length 5 the last index is 4
                new_times[new_nr_cum_values - i - 1] = times[nr_cum_values - keep_every_nth_val * i - 1]; // if i = new_nr_cum_values - 1 = nr_cum_values / keep_ever_nth_val + 1 - 1 = nr_cum_values / keep .. so we are accessing times[-1] but why I dont get it
            }
            m_vec = new_U_L;
            times = new_times;
            cout << "ADAPTED NEW WRITE INTERVAL: write_interval_old = " << write_interval;
            write_density = max(dt / new_write_interval, min_write_density);
            write_interval =  dt / write_density;
            cout  << "  write_interval_new = " << new_write_interval << endl;
            // If we make the write_interval that much slower we also adapt the minimum cum nr
            min_m_nr = new_nr_cum_values;      // I think this is better, it should be some kind of initial minimum number of cumulant values
            cout << "new min cum nr is  min_m_nr = " << new_nr_cum_values << endl;
            // we decided to rewrite the file, it will save space and make the after simulation validation easier
            rewrite_file();
            // as soon as the autocorrelation time takes over the autocorrelation time should decide when to stop
            // this is all right?
        }
    }

    void rewrite_file() {
        int new_nr_cum_values = m_vec.size();
        close_stream();
        open_stream(filepath);
        ofile << "t,m_vec" << endl;
        for(int i = 0; i < new_nr_cum_values; i++) {
            ofile << times[i] << "," << m_vec[i] << endl;
        }
    }

};

template <class system, class State>
class m_vectorial_equilibration_observer_adaptive: public obsver<system, State> {
    // does the same stuff as the normal equilibration observer but adapts the write density
    typedef obsver<system, State> obsver;
    using obsver::ofile;
    using obsver::open_stream;
    using obsver::open_app_stream;
    using obsver::close_stream;
    double write_interval = 1;
    double timepoint = 0;
    vector<pair<double, double>> m_vec{};
    vector<double> times{};
    int min_m_nr = 500;
    int m_nr_gpu = 2000;                 // nr of cum values at which the calculation switches to gpu
    double write_density = 1.0 /
                           100.0;      // we increase the write density to get a more useful value for the autocorrelation function?
    double min_write_density = 1.0 / 10000.0;   // I think 500 was already a bit slow for my taste
    double dt = 0.01;
    double dt_half;
    double equil_cutoff = 0.1;              // since the equilibration might influce the mean of m_vec a lot we cut a certain portion of m_vec values
    double max_error = 0.001;
    int m_nr = 0;                         // current number in the averageing process
    bool equilibrated = false;                      // for the usecase of the quench with dynamic equilibration
    double eval_factor = 10000.0;
    fs::path filepath;

public:
    m_vectorial_equilibration_observer_adaptive(int min_m_nr) : min_m_nr(min_m_nr) {}

    m_vectorial_equilibration_observer_adaptive(int min_m_nr, double write_density) : min_m_nr(min_m_nr),
                                                                                      write_density(write_density) {}

    m_vectorial_equilibration_observer_adaptive(int min_cum_nr, double write_density, double equil_cutoff) :
            min_m_nr(min_cum_nr), write_density(write_density), equil_cutoff(equil_cutoff) {}

    void init(fs::path path, map<Parameter, double> &paras, const system &sys) override {
        cout << "initializing m equilibration observer" << endl;
        int run_nr = (int) paras[Parameter::run_nr];
        max_error = paras[Parameter::equil_error];
        if (paras[Parameter::min_cum_nr]) {
            min_m_nr = (int) paras[Parameter::min_cum_nr];
        }
        timepoint = paras[Parameter::start_time];
        dt = paras[Parameter::dt];
        equilibrated = false;
        // we also need to reset m_vec and times, dont we?
        m_vec = vector<pair<double, double>>{};
        times = vector<double>{};
        close_stream();

        bool pick_up = (paras[Parameter::random_init] == -1.0);
        if (!pick_up) {
            filepath = path / (obsver::construct_filename(run_nr) + ".mag");
            open_stream(filepath);
            ofile << "t;m" << endl;
        } else {
            path += ".mag";
            filepath = path;
            readMagFromFile(path, m_vec, times);
            write_interval = times[1] - times[0];
            write_density = dt / write_interval;    // We should adapt the write interval of the file that we read
            cout << "picked up write density of " << write_density << "  corresponding to a write interval of "
                 << write_interval;
            // the other observers that dont read the values in actually have a problem with defining the timepoint?
            // problem is that you set the timepoint here to the back of m_vec, but it might be that the .csv
            // file has not written that many states because it somehow didnt finish...
            timepoint += write_interval;
            cout << "Pickup, setting timepoint of " << this->get_name() << " to " << timepoint << endl;
            open_app_stream(path);
        }
        cout << this->get_name() << " init called" << endl;


        // I think the starting write interval should be every one hundred steps
        dt_half = dt / 2.0;
        write_interval = dt / write_density;
        cout << "I am sure the write interval is " << write_interval << endl;
    }

    string get_name() override {
        return "m vectorial equilibration observer";
    }

    void operator()(system &sys, const State &x, double t) override {
        if (t > timepoint - dt_half) {
            // advancing the cum nr
            m_nr++;
            // we calculate the cumulant and write it down
            thrust::host_vector<pair<double, double>> m_val_vec = sys.calc_m_vec_vectorial(x);
            // with this observer we could think about writing at the end of the simulation, or at least rewriting
            // the file if we changed the stepsize
            // How do we do this, we cannot save every m value... I think we will try for now?
            ofile << t << ";";
            for (auto m_val: m_val_vec) {
                ofile << m_val.first << "," << m_val.second << ",";
            }
            ofile << endl;
            // add the cumulant and the times to the vectors to keep track
            // m_vec.push_back(m_val);
            m_vec.insert(m_vec.end(), m_val_vec.begin(), m_val_vec.end());
            times.push_back(t);
            // now we want to see if we have to adjust the write interval
            // we have to make sure that we got 5 fresh cum values
            if (!equilibrated) {
                int nr_m_measurements = m_vec.size() / m_val_vec.size();      // check how often we measured
                int nr_m_values = m_vec.size();
                // cout << "nr_m_values: " << nr_m_values << endl;
                // now use the last avg_nr of cum values to calculate a mean m_vec
                // use transform reduce?
                // does it work like this? m_vec.end() - avg_nr is the n-th last value in the vector?
                if (nr_m_measurements >= min_m_nr) {
                    int eval_every_n_steps = (int) (eval_factor * write_density) +
                                             1;   // I think this + 1 makes everything nicer? Or does it?
                    if (m_vec.size() % eval_every_n_steps == 0) {
                        // We dont need to calculate the error and stuff everytime we write down a m_vec
                        // IF the density is small, we should write down more often
                        // the equilibration phase might influence the mean a lot, should we cut of the first x% of the values?
                        int min_ind = (int) (equil_cutoff * nr_m_values);
                        int nr_values_to_use = nr_m_values - min_ind;
                        // we calculate the Binder cumulant according to the ergodic hyptheses from all the ms we extracted
                        // TODO if you have time you can do this on GPU but it shouldnt be to useful
                        pair<double, double> *m = &m_vec[min_ind];
                        // cout << "min ind " << min_ind << endl;
                        // cout << "nr_values_to_use " << nr_values_to_use << endl;
                        double m_L2 = std::transform_reduce(m, m + nr_values_to_use,
                                                            0.0, // initial value for the reduction (sum)
                                                            std::plus<double>(),
                                                            [](pair<double, double> m_val) -> double {
                                                                return m_val.first * m_val.first +
                                                                       m_val.second * m_val.second;
                                                            });
                        m_L2 /= (double) nr_values_to_use;
                        // TODO strictly speaking we would have to calculate autocorrelation times for m² and m⁴ but I
                        // want to get that running first
                        double m_L2_err = std::transform_reduce(m, m + nr_values_to_use,
                                                                0.0, // initial value for the reduction (sum)
                                                                std::plus<double>(), // transformation (square)
                                                                [&m_L2](pair<double, double> m_val) {
                                                                    return pow((m_val.first * m_val.first +
                                                                                m_val.second * m_val.second - m_L2), 2);
                                                                });
                        m_L2_err /= (double) pow(nr_values_to_use, 2);
                        m_L2_err = sqrt(m_L2_err);
                        double m_L4 = std::transform_reduce(m, m + nr_values_to_use,
                                                            0.0, // initial value for the reduction (sum)
                                                            std::plus<>(), // transformation (square)
                                                            [](pair<double, double> m_val) {
                                                                return (pow(m_val.first * m_val.first +
                                                                            m_val.second * m_val.second, 2));
                                                            });
                        m_L4 /= nr_values_to_use;
                        double m_L4_err = std::transform_reduce(m, m + nr_values_to_use,
                                                                0.0, // initial value for the reduction (sum)
                                                                std::plus<>(), // transformation (square)
                                                                [&m_L4](pair<double, double> m_val) {
                                                                    return pow(pow(m_val.first * m_val.first +
                                                                                   m_val.second * m_val.second, 2) -
                                                                               m_L4, 2);
                                                                });
                        m_L4_err /= (double) pow(nr_values_to_use, 2);
                        m_L4_err = sqrt(m_L4_err);
                        double U_L = m_L4 / (m_L2 * m_L2);
                        double U_L_error = sqrt(
                                pow(1 / m_L2 / m_L2 * m_L4_err, 2) + pow(2 * m_L4 / pow(m_L2, 3) * m_L2_err, 2));

                        double U_L_variance = U_L_error * U_L_error;
                        double rel_stddev_total = U_L_error / U_L;


                        pair<double, double> avg_m_vec = std::reduce(m, m + nr_values_to_use, pair<double, double>(),
                                                                     add_pair_functor<double>()); // / (double)nr_values_to_use;
                        avg_m_vec.first /= (double) nr_values_to_use;
                        avg_m_vec.second /= (double) nr_values_to_use;
                        cout << "avg_m_vec = (" << avg_m_vec.first << ", " << avg_m_vec.second << ")" << endl;

                        double avg_m = sqrt(avg_m_vec.first * avg_m_vec.first + avg_m_vec.second * avg_m_vec.second);
                        cout << "avg_m = " << avg_m << endl;
                        cout << "m_vec variance: " << U_L_variance << endl;
                        cout << "U_L: " << U_L << endl;
                        cout << "rel_stddev_total = " << rel_stddev_total << endl;

                        if (rel_stddev_total < max_error) {
                            cout << "U_L = " << U_L << " +- " << rel_stddev_total * U_L << endl;

                            // we set the system to be equilibrated
                            cout << "The system equilibrated, the equilibration lastet to t = " << t << endl;
                            sys.set_equilibration(t);
                            // write cumulant average
                            append_parameter(filepath, "U_L", U_L);
                            append_parameter(filepath, "m", avg_m);
                            append_parameter(filepath, "m2", m_L2);
                            append_parameter(filepath, "m4", m_L4);

                            // write error
                            append_parameter(filepath, "U_L_error", rel_stddev_total);
                            // write autocorrelation time
                            // append_parameter(filepath, "autocorrelation_time_U_L", autocorr_time);
                            // once we did this, we dont want to do that a second time?
                            equilibrated = true;
                        }
                    }
                }
            }
            timepoint += write_interval;
        }
    }
};

template <class system, class State>
class m_equilibration_observer_adaptive_autocorr: public obsver<system, State>{
    // does the same stuff as the normal equilibration observer but adapts the write density
    typedef obsver<system, State> obsver;
    using obsver::ofile;
    using obsver::open_stream;
    using obsver::open_app_stream;
    using obsver::close_stream;
    double write_interval = 1;
    double timepoint = 0;
    vector<double> m_vec{};
    vector<double> times{};
    int min_m_nr = 500;
    int m_nr_gpu = 2000;                 // nr of cum values at which the calculation switches to gpu
    double write_density = 1.0 / 100.0;      // we increase the write density to get a more useful value for the autocorrelation function?
    double min_write_density = 1.0 / 10000.0;   // I think 500 was already a bit slow for my taste
    double dt = 0.01;
    double dt_half;
    double equil_cutoff = 0.1;              // since the equilibration might influce the mean of m_vec a lot we cut a certain portion of m_vec values
    double max_error= 0.001;
    int m_nr = 0;                         // current number in the averageing process
    bool equilibrated = false;                      // for the usecase of the quench with dynamic equilibration
    double eval_factor = 10000.0;
    fs::path filepath;
    double max_moving_factor = 0.005;

public:
    m_equilibration_observer_adaptive_autocorr(int min_m_nr) : min_m_nr(min_m_nr) {}

    m_equilibration_observer_adaptive_autocorr(int min_m_nr, double write_density) : min_m_nr(min_m_nr), write_density(write_density) {}

    m_equilibration_observer_adaptive_autocorr(int min_cum_nr, double write_density, double equil_cutoff):
            min_m_nr(min_cum_nr), write_density(write_density), equil_cutoff(equil_cutoff) {}

    void init(fs::path path, map<Parameter, double>& paras, const system &sys) override {
        cout << "initializing m equilibration observer" << endl;
        int run_nr = (int)paras[Parameter::run_nr];
        max_error = paras[Parameter::equil_error];
        if(paras[Parameter::moving_factor] != 0.0) {
            // check that we do not set the moving factor to be zero
            max_moving_factor = paras[Parameter::moving_factor];
        }
        if (paras[Parameter::min_cum_nr]){
            min_m_nr = (int)paras[Parameter::min_cum_nr];
        }
        timepoint = paras[Parameter::start_time];
        dt = paras[Parameter::dt];
        equilibrated = false;
        // we also need to reset m_vec and times, dont we?
        m_vec = vector<double>{};
        times = vector<double>{};
        close_stream();

        bool pick_up = (paras[Parameter::random_init]  == -1.0);
        if(!pick_up) {
            filepath = path / (obsver::construct_filename(run_nr) + ".mag");
            open_stream(filepath);
            ofile << "t;m" << endl;
        } else {
            path += ".mag";
            filepath = path;
            readMagFromFile(path, m_vec, times);
            write_interval = times[1] - times[0];
            write_density = dt / write_interval;    // We should adapt the write interval of the file that we read
            cout << "picked up write density of " << write_density << "  corresponding to a write interval of " << write_interval;
            // the other observers that dont read the values in actually have a problem with defining the timepoint?
            // problem is that you set the timepoint here to the back of m_vec, but it might be that the .csv
            // file has not written that many states because it somehow didnt finish...
            timepoint += write_interval;
            cout << "Pickup, setting timepoint of " << this->get_name()  << " to " << timepoint << endl;
            open_app_stream(path);
        }
        cout << this->get_name() << " init called" << endl;


        // I think the starting write interval should be every one hundred steps
        dt_half = dt / 2.0;
        write_interval = dt / write_density;
        cout << "I am sure the write interval is " << write_interval << endl;
    }

    string get_name() override {
        return "m equilibration observer";
    }

    void operator()(system &sys, const State &x , double t ) override {
        if(t > timepoint - dt_half) {
            // advancing the cum nr
            m_nr++;
            // we calculate the cumulant and write it down
            thrust::host_vector<double> m_val_vec = sys.calc_m_vec(x);
            // with this observer we could think about writing at the end of the simulation, or at least rewriting
            // the file if we changed the stepsize
            // How do we do this, we cannot save every m value... I think we will try for now?
            ofile << t << ";";
            for (double m_val : m_val_vec) {
                ofile << m_val << ",";
            }
            ofile << endl;
            // add the cumulant and the times to the vectors to keep track
            // m_vec.push_back(m_val);
            m_vec.insert(m_vec.end(), m_val_vec.begin(), m_val_vec.end());
            times.push_back(t);
            // now we want to see if we have to adjust the write interval
            // we have to make sure that we got 5 fresh cum values
            if(!equilibrated) {
                int nr_m_measurements = m_vec.size() / m_val_vec.size();      // check how often we measured
                int nr_m_values = m_vec.size();
                // cout << "nr_m_values: " << nr_m_values << endl;
                // now use the last avg_nr of cum values to calculate a mean m_vec
                // use transform reduce?
                // does it work like this? m_vec.end() - avg_nr is the n-th last value in the vector?
                if(nr_m_measurements >= min_m_nr){
                    int eval_every_n_steps = (int)(eval_factor * write_density) + 1;   // I think this + 1 makes everything nicer? Or does it?
                    if(m_vec.size() % eval_every_n_steps == 0) {
                        // We dont need to calculate the error and stuff everytime we write down a m_vec
                        // IF the density is small, we should write down more often
                        // the equilibration phase might influence the mean a lot, should we cut of the first x% of the values?
                        int min_ind = (int)(equil_cutoff * nr_m_values);
                        int min_ind_measurements = (int)(equil_cutoff * nr_m_measurements);
                        int nr_values_to_use = nr_m_values - min_ind;
                        int nr_measurements_to_use = nr_m_measurements - min_ind_measurements;
                        // we calculate the Binder cumulant according to the ergodic hyptheses from all the ms we extracted
                        // TODO if you have time you can do this on GPU but it shouldnt be to useful
                        double* m = &m_vec[min_ind];
                        // cout << "min ind " << min_ind << endl;
                        // cout << "nr_values_to_use " << nr_values_to_use << endl;
                        // cout << "m_vector:" << endl;
                        // print_vector(m_vec);
                        // cout << "m array with nr_values_to use :" << nr_values_to_use  << endl;
                        // print_array(m, nr_values_to_use);
                        // cout << endl << endl;
                        // so... we need an new array, or we should just from beforehand use a vector of vectors
                        int nr_subsystems = m_val_vec.size();
                        vector<vector<double>> m_pow(nr_subsystems, std::vector<double>(nr_measurements_to_use, 0.0));
                        for(int t_ind = 0; t_ind < nr_measurements_to_use; t_ind++){
                            for(int sys_ind = 0; sys_ind < nr_subsystems; sys_ind++) {
                                // OMG this is going to be so slow
                                m_pow[sys_ind][t_ind] = pow(m[sys_ind  + t_ind * nr_subsystems], 2);
                            }
                        }
                        // cout << "m2" << endl;
                        // print_vector_of_vectors(m_pow);
                        // sleep(10);

                        vector<double> m2_avgs(nr_subsystems, 0.0);     // this will save the averages for all systems
                        vector<double> m2_errors(nr_subsystems, 0.0);   // this will save the errors for all systems     
                        
                        for(int sys_ind = 0; sys_ind < nr_subsystems; sys_ind++) {
                            // OMG this is going to be so slow
                            double* m2_arr = &m_pow[sys_ind][0];
                            // this should be all the m2s for all timepoints for the current system, for this I want to calculate the
                            // autocorrelation time now
                            double autocorr_time = get_autocorrtime_fft(m2_arr, nr_measurements_to_use, write_interval);
                            // cout << "autocorr_time m2 " << sys_ind << " = " << autocorr_time << endl;
                            // okay but we still need the mean aswell as the stddev of the distribution as I call it
                            double m_L2 = std::reduce(m2_arr, m2_arr + nr_measurements_to_use,
                                                                0.0, // initial value for the reduction (sum)
                                                                std::plus<double>()) / (double) nr_measurements_to_use;
                            double m_L2_err = std::transform_reduce(m2_arr, m2_arr + nr_measurements_to_use,
                                                                    0.0, // initial value for the reduction (sum)
                                                                    std::plus<double>(), // transformation (square)
                                                                    [&m_L2](double m2_val) { return pow((m2_val - m_L2), 2) ; });
                            double m2_variance = 2 * autocorr_time / ((nr_measurements_to_use) * write_interval) * m_L2_err;
                            // and these should now be saved
                            m2_avgs[sys_ind] = m_L2;
                            m2_errors[sys_ind] = m2_variance;
                        }
                        
                        // okay we need to do this whole stuff also for m4
                        for(int t_ind = 0; t_ind < nr_measurements_to_use; t_ind++){
                            for(int sys_ind = 0; sys_ind < nr_subsystems; sys_ind++) {
                                // we can use here the already calculated powers of m
                                m_pow[sys_ind][t_ind] = pow(m_pow[sys_ind][t_ind], 2);
                            }
                        }

                        // print_vector_of_vectors(m_pow);
                        // sleep(10);

                        vector<double> m4_avgs(nr_subsystems, 0.0);     // this will save the averages for all systems
                        vector<double> m4_errors(nr_subsystems, 0.0);   // this will save the errors for all systems     

                        for(int sys_ind = 0; sys_ind < nr_subsystems; sys_ind++) {
                            // OMG this is going to be so slow
                            double* m4_arr = &m_pow[sys_ind][0];
                            // this should be all the m2s for all timepoints for the current system, for this I want to calculate the
                            // autocorrelation time now
                            double autocorr_time = get_autocorrtime_fft(m4_arr, nr_measurements_to_use, write_interval);

                            // okay but we still need the mean aswell as the stddev of the distribution as I call it
                            double m_L4 = std::reduce(m4_arr, m4_arr + nr_measurements_to_use,
                                                      0.0, // initial value for the reduction (sum)
                                                      std::plus<double>()) / (double) nr_measurements_to_use;
                            double m_L4_err = std::transform_reduce(m4_arr, m4_arr + nr_measurements_to_use,
                                                                    0.0, // initial value for the reduction (sum)
                                                                    std::plus<double>(), // transformation (square)
                                                                    [&m_L4](double m4_val) { return pow((m4_val - m_L4), 2) ; });
                            double m4_variance = 2 * autocorr_time / ((nr_measurements_to_use) * write_interval) * m_L4_err;
                            // and these should now be saved
                            m4_avgs[sys_ind] = m_L4;
                            m4_errors[sys_ind] = m4_variance;
                        }
                        // cout << "m2_avgs vector:" << endl;
                        // print_vector(m2_avgs);
                        // cout << "m2_error vector:" << endl;
                        // print_vector(m2_errors);
//
                        // cout << "m4_avgs vector:" << endl;
                        // print_vector(m4_avgs);
                        // cout << "m4_error vector:" << endl;
                        // print_vector(m4_errors);
                        // the averages are not much now
                        double m_L2 = std::reduce(m2_avgs.begin(), m2_avgs.end(),
                                                  0.0, // initial value for the reduction (sum)
                                                  std::plus<double>()) / (double) nr_subsystems;
                        double m_L4 = std::reduce(m4_avgs.begin(), m4_avgs.end(),
                                                  0.0, // initial value for the reduction (sum)
                                                  std::plus<double>()) / (double) nr_subsystems;
                        // the errrors are standard error propagation
                        double m_L2_err = sqrt(transform_reduce(m2_errors.begin(), m2_errors.end(),
                                                                0.0, // initial value for the reduction (sum)
                                                                std::plus<double>(), // transformation (square)
                                                                [](double m2_err) { return pow(m2_err, 2);}) / (double)nr_subsystems);
                        double m_L4_err = sqrt(transform_reduce(m4_errors.begin(), m4_errors.end(),
                                                                0.0, // initial value for the reduction (sum)
                                                                std::plus<double>(), // transformation (square)
                                                                [](double m4_err) { return pow(m4_err, 2);}) / (double)nr_subsystems);

                        double U_L = m_L4 / (m_L2 * m_L2);
                        double U_L_error = sqrt(pow(1 / m_L2 / m_L2 * m_L4_err, 2) + pow(2 * m_L4 / pow(m_L2, 3) * m_L2_err, 2));

                        double U_L_variance = U_L_error * U_L_error;
                        double rel_stddev_total = U_L_error / U_L;

/*
                        std::vector<double> diff_total(nr_m_values - min_ind);

                        std::transform(m_vec.begin() + min_ind, m_vec.end(), diff_total.begin(), [U_L](double x) { return x - U_L; });

                        double dist_var = std::inner_product(diff_total.begin(), diff_total.end(), diff_total.begin(), 0.0) / (nr_m_values - min_ind);
                        // this is the variance of the distribution (in wrong because the values are correlated), so we
                        // need the autocorrelation time to adjust for this. The function works with arrays
                        double* U_L_arr = &m_vec[min_ind];
                        double autocorr_time;

                        if(nr_m_values - min_ind > m_nr_gpu) {
                            autocorr_time = get_autocorrtime_fft(U_L_arr, nr_m_values - min_ind, write_interval);
                        } else {
                            autocorr_time = get_autocorrtime(U_L_arr, nr_m_values - min_ind, write_interval);  // actually ds is just the write interval? which should be 1 or something like this
                        }
                        double U_L_variance = 2 * autocorr_time / ((nr_m_values - min_ind) * write_interval) * dist_var;
                        double rel_stddev_total = sqrt(U_L_variance) / U_L;

                        cout << "Distriubtion variance: " << dist_var << endl;
                        cout << "autocorrelation time: " << autocorr_time << endl;
                        */
                        double avg_m = std::reduce(m, m + nr_values_to_use) / (double)nr_values_to_use;
                        cout << "avg_m = " << avg_m << endl;
                        cout << "m_vec variance: " << U_L_variance << endl;
                        cout << "U_L: " << U_L << endl;
                        cout << "rel_stddev_total = " << rel_stddev_total << endl;

                        // sleep(10);
                        // Okay we found out that this is not the real error, so we need a function that calculates the error
                        // or first a function that calculates the integrated autocorrleation time
                        // We indded have two auto correlation times for the two directions. Could they be individually sized?
                        // the question is now if we want to extract U_L and its error if we are equilibrated
                        // we definitely should? But how and where? somehow into the parameter file?
                        if(rel_stddev_total < max_error) {
                            // I think in this case I want to see wether there is still a significant upwards or
                            // downards trend to ensure the system is equilibrated
                            // the question is what recent means, the last 20% of values? It shouldnt be to much
                            // that it wouldnt be recent anymore and it should be so few that statistical fluctuations
                            // would average out
                            // to be honest, the trend stuff is probably the prestuff we should do before even considering
                            // to calculate the error? Error calculation after all is pretty expensive
                            double moving_factor = getMovingFactor(nr_m_values, avg_m, min_ind);

                            cout << "U_L = " << U_L << " +- " << rel_stddev_total * U_L << endl;
                            cout << "MOVING FACTOR = " << moving_factor << endl;

                            if(moving_factor < max_moving_factor) {
                                // we set the system to be equilibrated
                                cout << "The system equilibrated, the equilibration lastet to t = " << t << endl;
                                sys.set_equilibration(t);
                                // write cumulant average
                                append_parameter(filepath, "U_L", U_L);
                                append_parameter(filepath, "m", avg_m);
                                append_parameter(filepath, "m2", m_L2);
                                append_parameter(filepath, "m4", m_L4);

                                // write error
                                append_parameter(filepath, "U_L_error", rel_stddev_total);
                                // write autocorrelation time
                                // append_parameter(filepath, "autocorrelation_time_U_L", autocorr_time);
                                // Also write the moving factor
                                append_parameter(filepath, "moving_factor_U_L", moving_factor);
                                // once we did this, we dont want to do that a second time?
                                equilibrated = true;
                            }
                        }
                        if(!equilibrated) { // TODO
                            // adapt_write_interval(nr_m_values, autocorr_time);
                        }
                    }
                }
            }
            timepoint += write_interval;
        }
    }

    double getMovingFactorOld(int nr_cum_values, long double avg_U_L) {
        int recent_ind = (int) (nr_cum_values * 0.9);
        double* recent_U_L = &m_vec[recent_ind];
        // the idea is to check if the differences between the values show a trend to be positive
        // or negative. But I think that comes done to calculating the deviation of the recent
        // average from the total average
        // small differences will
        // Idea would be to compare this deviation to the deviation per step
        // If we are in a high temperature state, the deviation per step will be relatively larger
        // If should not be dependent on the stepsize
        // first of all compute the avergage absolute difference per step
        int recent_size = nr_cum_values - recent_ind;
        double recent_avg_U_L = reduce(recent_U_L, recent_U_L + recent_size) / (double)recent_size;
        double mean_abs_delta_U = meanAbsDifference(recent_U_L, recent_size);
        // and we need the recent mean
        double moving_factor = fabs(avg_U_L - recent_avg_U_L) / ((double)recent_size * mean_abs_delta_U);

        cout << "recent_avg_U_L = " << recent_avg_U_L << endl;
        cout << "|U_L_recent - m_vec| = " << avg_U_L - recent_avg_U_L << endl;
        cout << "mean abs delta = " << mean_abs_delta_U << endl;
        return moving_factor;
    }

    double getMovingFactor(int nr_values, long double avg, int min_ind) {
        int recent_ind = (int) (nr_values * 0.8);
        // the idea is to check if the differences between the values show a trend to be positive
        // or negative. But I think that comes done to calculating the deviation of the recent
        // average from the total average
        // small differences will
        // Idea would be to compare this deviation to the deviation per step
        // If we are in a high temperature state, the deviation per step will be relatively larger
        // If should not be dependent on the stepsize
        // first of all compute the avergage absolute difference per step
        int recent_size = nr_values - recent_ind;
        double var_start = reduce(m_vec.begin() + min_ind, m_vec.begin() + recent_ind) / (double)(recent_ind - min_ind);
        double var_end = avg;

        double* recent_U_L = &m_vec[recent_ind];      // we need it still for the mean delta?
        double mean_abs_delta_U = meanAbsDifference(recent_U_L, recent_size);
        // and we need the recent mean
        double moving_factor = fabs(var_end - var_start) / ((double)recent_size * mean_abs_delta_U);

        cout << "var_start = " << var_start << endl;
        cout << "|var_end - var_start| = " << avg - var_start << endl;
        cout << "mean abs delta = " << mean_abs_delta_U << endl;
        return moving_factor;
    }

    void adapt_write_interval(int nr_cum_values, double autocorr_time, int values_per_autocorr_time=20) {
        // we want to adapt the write density on the autocorrelation time. But everytime we do this
        // we have to delete the values that have smaller spacing between them (the intermediate values)
        // It should be some kind of switch statement or something like this?
        // We want to write approximately 10 cumulant values during one autocorrelation time? Or maybe 100 ?
        // I think 10 should be fine
        double autocorr_time_floored = pow(10.0, floor(log10(autocorr_time)));      // we floor the autocorrelation time to the next power of 10, so 3000 would get rounded to 1000 and 101 to 100
        // we do this as we do not want to continuously change our stepping but in discrete steps
        double new_write_interval = autocorr_time_floored / (double)values_per_autocorr_time;    // if autocorr is 1000, the write interval is 100, using 100 writes per autocorr time now because ten is not satisfying somehow
        if(new_write_interval > write_interval) {
            // in this case we have to drop values accordingly to new write interval
            int keep_every_nth_val = (int)(new_write_interval / write_interval); // they should yield a glatt integer since we only work with powers of ten write interval might have been 1 and now it is 10 so we keep every 10th value
            int new_nr_cum_values = nr_cum_values / keep_every_nth_val;         // if we have 101 values with a spacing of 1, we want to keep 101, 91, ..., 1, making up for 11 values
            if (nr_cum_values % keep_every_nth_val != 0) {
                // if this is the case we can extract one more value
                new_nr_cum_values++;
            }
/*            cout << "the old number of cum values is " << nr_cum_values << endl;
            cout << "the new number of cum values is " << new_nr_cum_values << endl;*/
            // if we have 100 values with a spacing of 1, those would be the indices 0 to 99. We want to keep 99 and 89 (distance of 10 ds), 79, ..., 9 making up 10 values
            // if we have 101 values with a spacing of 1, those would have the indices 0 to 100. We would want to keep 100, 90, ..., 10, 0 so eleven values? Right now we would be missing out on 0, which would not be two bad...
            vector<double> new_U_L(new_nr_cum_values);
            vector<double> new_times(new_nr_cum_values);
            for(int i = 0; i < new_nr_cum_values; i++) {
                // TODO is the plus one correct? If you get memory errors look here?
                // okay so we work from behind because instead of keeping 0 - 100 we want to keep 1 - 101. But the order has to stay the same for the correlation function
                // cout << "i = " << i << "   nr_cum_values - keep_every_nth_val * i - 1 = " << nr_cum_values - keep_every_nth_val * i - 1 << "    m_vec = " << m_vec[nr_cum_values - keep_every_nth_val * i - 1] << endl;
                new_U_L[new_nr_cum_values - i - 1] = m_vec[nr_cum_values - keep_every_nth_val * i - 1];       // - one because in a vector of length 5 the last index is 4
                new_times[new_nr_cum_values - i - 1] = times[nr_cum_values - keep_every_nth_val * i - 1]; // if i = new_nr_cum_values - 1 = nr_cum_values / keep_ever_nth_val + 1 - 1 = nr_cum_values / keep .. so we are accessing times[-1] but why I dont get it
            }
            m_vec = new_U_L;
            times = new_times;
            cout << "ADAPTED NEW WRITE INTERVAL: write_interval_old = " << write_interval;
            write_density = max(dt / new_write_interval, min_write_density);
            write_interval =  dt / write_density;
            cout  << "  write_interval_new = " << new_write_interval << endl;
            // If we make the write_interval that much slower we also adapt the minimum cum nr
            min_m_nr = new_nr_cum_values;      // I think this is better, it should be some kind of initial minimum number of cumulant values
            cout << "new min cum nr is  min_m_nr = " << new_nr_cum_values << endl;
            // we decided to rewrite the file, it will save space and make the after simulation validation easier
            rewrite_file();
            // as soon as the autocorrelation time takes over the autocorrelation time should decide when to stop
            // this is all right?
        }
    }

    void rewrite_file() {
        int new_nr_cum_values = m_vec.size();
        close_stream();
        open_stream(filepath);
        ofile << "t,m_vec" << endl;
        for(int i = 0; i < new_nr_cum_values; i++) {
            ofile << times[i] << "," << m_vec[i] << endl;
        }
    }

};

template <class system, class State>
class corr_observer : public obsver<system, State>{
    // Okay the observing pattern could be totally wild, i probably somehow have to initialize the observer
    // outside of the simulation class We definetely need its own constructor here
    typedef obsver<system, State> obsver;
    using obsver::ofile;
    using obsver::open_stream;
    using obsver::close_stream;
    int nr_values;
    double write_interval = 1;
    double timepoint = 0;
public:
    corr_observer(int nr_values) : nr_values(nr_values) {
    }
    void init(fs::path folderpath, map<Parameter, double>& paras, const system &sys) override {
        int run_nr = (int)paras[Parameter::run_nr];
        timepoint = 0.0;

        close_stream();
        open_stream(folderpath / (obsver::construct_filename(run_nr) + ".corr"));
        cout << "corr observer init called" << endl;
        ofile << "t,xix,xiy" << endl;

        // I think this will have less performance impact than an if statement catching the first observer operation
        // Make sure to use this observer only with systems that have a get_end_T method
        double end_t = paras[end_time];
        write_interval = end_t / (double)nr_values;
    }

    string get_name() override {
        return "corr observer";
    }

    void operator()(system &sys, const State &x , double t ) override {
        if(t > timepoint) {
            double xix, xiy;
            sys.calc_xi(x, xix, xiy);
            ofile << t << "," << xix << "," << xiy << endl;
            timepoint += write_interval;
        }
    }
};

template <class system, class State>
class corr_equilibration_observer: public obsver<system, State>{
    // This observer will be a mixture of the cum equilibration observer and the new density observer
    // It will calculate xi with a fixed density during equilibration phase and then the specified number
    // of xi's during the quench.
    // The density during the equilibration does not have to change, or should it? I don't think so
    // Problem could be that if the equilibration takes very long, we calculate a whole lot of xi values during the
    // equilibratioon which slows down our simulation. An adaptive density doesn't look to good when plotting and has
    // the weird artefact that we write the correlation length at different times for different simulations
    // I think we will choose a mediocre density and start with this.
    typedef obsver<system, State> obsver;
    using obsver::ofile;
    using obsver::open_stream;
    using obsver::open_app_stream;
    using obsver::close_stream;
    double write_interval = 1;
    double timepoint = 0;
    vector<double> xix{};
    vector<double> xiy{};
    vector<double> times{};
    int min_corr_nr = 5000;
    int val_nr_gpu = 1000;                 // nr of cum values at which the calculation switches to gpu
    double dt = 0.01;
    double dt_half = dt / 2.0;
    double equil_cutoff = 0.1;              // since the equilibration might influce the mean of xi
    double density = 1.0 / 100.0;            // standard density writes once every 100 steps
    double error_factor = 10000.0;
    // TODO we have to judge whether this is large or not. The thing is the error is good for low temperature states to
    //  judge whether we are equilibrated but bad for high temperature states since we have large deviations
    // for high temperature states we would usually need a larger number of systems to judge the equilibration
    // TODO could there be a way of combining the temperature with the error? A higher temperature would allow
    // a larger error on xi and still judge it to be equilibrated. But this is again kind of handwavy, how large temperatures
    // would result in how large leeway?
    // also the error has the weird property that large systems with small nubmers of subsystems equilibrate later.
    // For the Tc-Binder cumulant calculation this error is very suitable since we are looking for a precise m_vec value anyway
    // here we are just looking to judge the system to be in thermal equilibrium
    // Would there be other ways of assessing the equilibration of our system? The energy? Will also fluctate, but maybe
    // not as strongly as the correlation length?
    // The good thing about the energy also is that the fluctuations scale with the number of lattice sites so large systems wont be
    // in disadvantage
    // Okay I think we will implement this for now and judge afterwards how good it works
    // the question remains if the relaxation of the energy also means that the correlation length is relaxed, I for my
    // case don't think so.
    // If we want to extract the correlation length during the run we have to cut the zero impuls, I think this has
    // to be adapted in the systems.
    double max_error= 0.05;
    bool equilibrated = false;                      // for the usecase of the quench with dynamic equilibration
    int nr_values;                  // nr of values that I want to be written down during the quench
    double quench_t;                // quench time, important to calculate the write interval during the quench
    size_t Lx;                      // sizes of the system. Important to validate if the extraction was meaningful
    size_t Ly;
    double xi_cap = 0.2;            // We cap the xi size to be at most 0.2 of the corresponding system size, otherwise the extraction is unreliable and
    // will distort the error estimation. If the total run was meaningful will then be judged afterwards.
    // (If the xi is now close to this threshold this probably means that our system was evolving in a state with larger correlation length than we can meaningful extract)
    // Okay short recap, we will now check evertime we calculate a xi if it is smaller than xi_cap * L and if it is not, we will cut it there
    // Afterwards in the automatic suite, we will calculate the avg xi and if this is smaller then for example xi_cap / 2, then we will accept the measurement
    // we will probably slightly influence the correlation lengths to be a bit smaller than the "real" value, since we do not include the statistical fluctuations where xi becomes larger than 0.2 L
    // The goal is that we hold those occasions very small so that the influence is negligeable
public:
    corr_equilibration_observer(int nr_values) : nr_values(nr_values) {
    }
    corr_equilibration_observer(int nr_values, int min_corr_nr, double density, double equil_cutoff):
    nr_values(nr_values), min_corr_nr(min_corr_nr), density(density), equil_cutoff(equil_cutoff) {}
    void init(fs::path path, map<Parameter, double>& paras, const system &sys) override {
        int run_nr = (int)paras[Parameter::run_nr];
        max_error = paras[Parameter::equil_error];
        equil_cutoff = paras[Parameter::equil_cutoff];
        min_corr_nr = (int)paras[Parameter::min_corr_nr];
        // the subsystems sizes for the xi cutting
        Lx = (size_t)paras[Parameter::subsystem_Lx];
        Ly = (size_t)paras[Parameter::subsystem_Ly];
        // the timepoint is only zero if we do not memory initialize
        timepoint = 0.0;
        equilibrated = false;
        // we also need to reset m_vec and times, dont we?
        xix = vector<double>{};
        xiy = vector<double>{};
        times = vector<double>{};
        close_stream();
        bool pick_up = (paras[Parameter::random_init]  == -1.0);
        if(!pick_up) {
            open_stream(path / (obsver::construct_filename(run_nr) + ".corr"));
            // We only want to do this if we are creating a new file?
            ofile << "t,xix,xiy" << endl;
        } else {
            // If we have the memory initialization we also want to load the correlation length values into cache
            // I thinkt we should do that before we decide to open this file also as output stream
            path += ".corr";
            cout << "Reading correlation lengths from file" << endl;
            readXiFromFile(path, xix, xiy, times);
            timepoint = times.back();
            cout << "successful!" << endl;
            open_app_stream(path);
        }
        cout << this->get_name() << " init called" << endl;


        // I think the starting write interval should be every one hundred steps
        dt = paras[Parameter::dt];
        dt_half = dt / 2.0;
        quench_t = sys.get_quench_time();
        write_interval = (1.0 / density) * dt;
    }

    string get_name() override {
        return "corr equilibration observer";
    }

    void operator()(system &sys, const State &x , double t ) override {
        if(t > timepoint - dt_half) {
            double xix_val, xiy_val;
            Singleton_timer::set_startpoint("Xi Calculation");
            sys.calc_xi(x, xix_val, xiy_val);
            Singleton_timer::set_endpoint("Xi Calculation");
            // The logic of the cutting should take place here or only for calculating the mean? I think we also cut it
            // but we might change that later
            xix_val = min(xix_val, xi_cap * Lx);
            xiy_val = min(xiy_val, xi_cap * Ly);
            ofile << t << "," << xix_val << "," << xiy_val << endl;
            // add the correlation lengths and the times to the vectors to keep track
            // we actually only need to keep track if we did not decide already that we are equilibrated?
            if(!equilibrated) {
                xix.push_back(xix_val);
                xiy.push_back(xiy_val);
                times.push_back(t);
                int nr_xi_values = xix.size();      // check how many corr values we already have, the number of xiy values equals the number of xix values
                // we lack the complete logic to change the stepsize since we said we use a constant density 
                if(nr_xi_values > min_corr_nr){
                    int error_every_n_steps = max((int)(error_factor * density), 1);
                    if(xix.size() % error_every_n_steps == 0) {
                        // if we reached the minimum number of values we check the error on the correlation lengths
                        int min_ind = (int)(equil_cutoff * nr_xi_values);
                        // we need to calculate the average aswell as the stddev for both directions
                        double avg_xix = accumulate(xix.begin() + min_ind, xix.end(), 0.0) / (double)(nr_xi_values - min_ind);
                        double avg_xiy = accumulate(xiy.begin() + min_ind, xiy.end(), 0.0) / (double)(nr_xi_values - min_ind);

                        std::vector<double> diff_xix_total(nr_xi_values - min_ind);
                        std::vector<double> diff_xiy_total(nr_xi_values - min_ind);

                        std::transform(xix.begin() + min_ind, xix.end(), diff_xix_total.begin(), [avg_xix](double x) { return x - avg_xix; });
                        std::transform(xiy.begin() + min_ind, xiy.end(), diff_xiy_total.begin(), [avg_xiy](double x) { return x - avg_xiy; });

                        double dist_var_x = std::inner_product(diff_xix_total.begin(), diff_xix_total.end(), diff_xix_total.begin(), 0.0) / (nr_xi_values - min_ind);
                        double dist_var_y = std::inner_product(diff_xiy_total.begin(), diff_xiy_total.end(), diff_xiy_total.begin(), 0.0) / (nr_xi_values - min_ind);

                        // autocorrelation time, works also for the correlation length
                        double* xix_arr = &xix[min_ind];
                        double* xiy_arr = &xiy[min_ind];

                        double autocorr_time_x;
                        double autocorr_time_y;

                        if(nr_xi_values - min_ind > val_nr_gpu) {
                            Singleton_timer::set_startpoint("Autocorrelation time fft calculation");
                            autocorr_time_x = get_autocorrtime_fft(xix_arr, nr_xi_values - min_ind, write_interval);
                            autocorr_time_y = get_autocorrtime_fft(xiy_arr, nr_xi_values - min_ind, write_interval);
                            Singleton_timer::set_endpoint("Autocorrelation time fft calculation");
                        } else {
                            Singleton_timer::set_startpoint("Autocorrelation time calculation");
                            autocorr_time_x = get_autocorrtime(xix_arr, nr_xi_values - min_ind, write_interval);
                            autocorr_time_y = get_autocorrtime(xiy_arr, nr_xi_values - min_ind, write_interval);  // actually ds is just the write interval? which should be 1 or something like this
                            Singleton_timer::set_endpoint("Autocorrelation time calculation");
                        }

                        double xix_variance = 2 * autocorr_time_x / ((nr_xi_values-min_ind) * write_interval) * dist_var_x;
                        double xiy_variance = 2 * autocorr_time_y / ((nr_xi_values-min_ind) * write_interval) * dist_var_y;


                        double rel_stddev_xix_total = sqrt(xix_variance) / avg_xix;
                        double rel_stddev_xiy_total = sqrt(xiy_variance) / avg_xiy;

                        cout << "xix = " << avg_xix << endl;
                        cout << "xiy = " << avg_xiy << endl;
                        cout << "autocorrelation time x: " << autocorr_time_x << endl;
                        cout << "autocorrelation time y: " << autocorr_time_y << endl;
                        cout << "xix_variance = " << xix_variance << endl;
                        cout << "xiy_variance = " << xiy_variance << endl;
                        cout << "rel_stddev_total xix = " << rel_stddev_xix_total << endl;
                        cout << "rel_stddev_total xiy = " << rel_stddev_xiy_total << endl;
                        // I mean you are calculating the errors anyway, you could also write them down? But takes some time ofc, but it would make it easier to see how long a certain simulation will still take...


                        // I would say if the mean of the two relative standard deviations satisfies the condition we are fine
                        double rel_stddev_total = 0.5 * (rel_stddev_xix_total + rel_stddev_xiy_total);

                        if(rel_stddev_total < max_error) {
                            cout << "The system equilibrated, the equilibration lastet to t = " << t << endl;
                            cout << "xix = " << avg_xix << " +- " << rel_stddev_xix_total * avg_xix << endl;
                            cout << "xiy = " << avg_xiy << " +- " << rel_stddev_xiy_total * avg_xiy << endl;
                            // we set the system to be equilibrated
                            sys.set_equilibration(t);           // For relaxation simulations this means that the simulation ends
                            // once we did this, we dont want to do that a second time?
                            equilibrated = true;
                            // the writ interval will now be the one that we destined for the quench, we just change it once here
                            write_interval = quench_t / (double)nr_values;
                        }
                    }
                }
            }
            timepoint += write_interval;
        }
    }

};

template <class system, class State>
class corr_equilibration_observer_adaptive: public obsver<system, State>{
    // This observer will be a mixture of the cum equilibration observer and the new density observer
    // It will calculate xi with a fixed density during equilibration phase and then the specified number
    // of xi's during the quench.
    // The density during the equilibration does not have to change, or should it? I don't think so
    // Problem could be that if the equilibration takes very long, we calculate a whole lot of xi values during the
    // equilibratioon which slows down our simulation. An adaptive density doesn't look to good when plotting and has
    // the weird artefact that we write the correlation length at different times for different simulations
    // I think we will choose a mediocre density and start with this.
    typedef obsver<system, State> obsver;
    using obsver::ofile;
    using obsver::open_stream;
    using obsver::open_app_stream;
    using obsver::close_stream;
    double write_interval = 1;
    double timepoint = 0;
    vector<double> xix{};
    vector<double> xiy{};
    vector<double> times{};
    int min_corr_nr = 5000;
    int val_nr_gpu = 1000;                 // nr of cum values at which the calculation switches to gpu
    double dt = 0.01;
    double dt_half = dt / 2.0;
    double equil_cutoff = 0.1;              // since the equilibration might influce the mean of xi
    double density = 1.0 / 100.0;            // standard density writes once every 100 steps
    double min_density = 1.0 / 10000.0;
    double error_factor = 10000.0;
    int observed_direction = 0;             // 0: both, 1: x, 2:y, for the error
    // TODO we have to judge whether this is large or not. The thing is the error is good for low temperature states to
    //  judge whether we are equilibrated but bad for high temperature states since we have large deviations
    // for high temperature states we would usually need a larger number of systems to judge the equilibration
    // TODO could there be a way of combining the temperature with the error? A higher temperature would allow
    // a larger error on xi and still judge it to be equilibrated. But this is again kind of handwavy, how large temperatures
    // would result in how large leeway?
    // also the error has the weird property that large systems with small nubmers of subsystems equilibrate later.
    // For the Tc-Binder cumulant calculation this error is very suitable since we are looking for a precise m_vec value anyway
    // here we are just looking to judge the system to be in thermal equilibrium
    // Would there be other ways of assessing the equilibration of our system? The energy? Will also fluctate, but maybe
    // not as strongly as the correlation length?
    // The good thing about the energy also is that the fluctuations scale with the number of lattice sites so large systems wont be
    // in disadvantage
    // Okay I think we will implement this for now and judge afterwards how good it works
    // the question remains if the relaxation of the energy also means that the correlation length is relaxed, I for my
    // case don't think so.
    // If we want to extract the correlation length during the run we have to cut the zero impuls, I think this has
    // to be adapted in the systems.
    double max_error= 0.05;
    bool equilibrated = false;                      // for the usecase of the quench with dynamic equilibration
    int nr_values;                  // nr of values that I want to be written down during the quench
    double quench_t;                // quench time, important to calculate the write interval during the quench
    size_t Lx;                      // sizes of the system. Important to validate if the extraction was meaningful
    size_t Ly;
    double xi_cap = 0.2;            // We cap the xi size to be at most 0.2 of the corresponding system size, otherwise the extraction is unreliable and
    fs::path filepath;              // path this observer writes to
    double max_moving_factor = 0.005;
    // will distort the error estimation. If the total run was meaningful will then be judged afterwards.
    // (If the xi is now close to this threshold this probably means that our system was evolving in a state with larger correlation length than we can meaningful extract)
    // Okay short recap, we will now check evertime we calculate a xi if it is smaller than xi_cap * L and if it is not, we will cut it there
    // Afterwards in the automatic suite, we will calculate the avg xi and if this is smaller then for example xi_cap / 2, then we will accept the measurement
    // we will probably slightly influence the correlation lengths to be a bit smaller than the "real" value, since we do not include the statistical fluctuations where xi becomes larger than 0.2 L
    // The goal is that we hold those occasions very small so that the influence is negligeable
public:
    corr_equilibration_observer_adaptive(int nr_values) : nr_values(nr_values) {
    }
    corr_equilibration_observer_adaptive(int nr_values, int min_corr_nr, double density, double equil_cutoff):
            nr_values(nr_values), min_corr_nr(min_corr_nr), density(density), equil_cutoff(equil_cutoff) {}
    void init(fs::path path, map<Parameter, double>& paras, const system &sys) override {
        int run_nr = (int)paras[Parameter::run_nr];
        max_error = paras[Parameter::equil_error];
        equil_cutoff = paras[Parameter::equil_cutoff];
        min_corr_nr = (int)paras[Parameter::min_corr_nr];
        max_moving_factor = paras[Parameter::moving_factor];
        // the subsystems sizes for the xi cutting
        Lx = (size_t)paras[Parameter::subsystem_Lx];
        Ly = (size_t)paras[Parameter::subsystem_Ly];
        observed_direction = (int)paras[Parameter::observed_direction];
        // the timepoint is only zero if we do not memory initialize
        timepoint = 0.0;
        equilibrated = false;
        // we also need to reset m_vec and times, dont we?
        xix = vector<double>{};
        xiy = vector<double>{};
        times = vector<double>{};
        close_stream();
        bool pick_up = (paras[Parameter::random_init]  == -1.0);
        if(!pick_up) {
            filepath = path / (obsver::construct_filename(run_nr) + ".corr");
            open_stream(filepath);
            // We only want to do this if we are creating a new file?
            ofile << "t,xix,xiy" << endl;
        } else {
            // If we have the memory initialization we also want to load the correlation length values into cache
            // I thinkt we should do that before we decide to open this file also as output stream
            path += ".corr";
            filepath = path;
            cout << "Reading correlation lengths from file" << endl;
            readXiFromFile(path, xix, xiy, times);
            timepoint = paras[Parameter::start_time];
            cout << "successful!" << endl;
            // we should adapt the write interval from the file that we read
            write_interval = times[times.size() - 1] - times[times.size() - 2];
            density = max(dt / write_interval, min_density);
            open_app_stream(filepath);
        }
        cout << this->get_name() << " init called" << endl;


        // I think the starting write interval should be every one hundred steps
        dt = paras[Parameter::dt];
        dt_half = dt / 2.0;
        quench_t = sys.get_quench_time();
        write_interval = (1.0 / density) * dt;
        cout << "write interval of corr equilibration observer adaptive is " << write_interval << endl;
        cout << "and the next timepoint is " << timepoint << endl;
    }

    string get_name() override {
        return "corr equilibration observer adaptive";
    }

    void operator()(system &sys, const State &x , double t ) override {
        if(t > timepoint - dt_half) {
            cout << "t > timepoint.... " << t << " > " << timepoint << endl;
            double xix_val, xiy_val;
            Singleton_timer::set_startpoint("Xi Calculation");
            sys.calc_xi(x, xix_val, xiy_val);
            Singleton_timer::set_endpoint("Xi Calculation");
            // The logic of the cutting should take place here or only for calculating the mean? I think we also cut it
            // but we might change that later
            xix_val = min(xix_val, xi_cap * Lx);
            xiy_val = min(xiy_val, xi_cap * Ly);
            ofile << t << "," << xix_val << "," << xiy_val << endl;
            // add the correlation lengths and the times to the vectors to keep track
            // we actually only need to keep track if we did not decide already that we are equilibrated?
            if(!equilibrated) {
                xix.push_back(xix_val);
                xiy.push_back(xiy_val);
                times.push_back(t);
                int nr_xi_values = xix.size();      // check how many corr values we already have, the number of xiy values equals the number of xix values
                // we lack the complete logic to change the stepsize since we said we use a constant density
                if(nr_xi_values > min_corr_nr){
                    int error_every_n_steps = (int)(error_factor * density) + 1;
                    if(xix.size() % error_every_n_steps == 0) {
                        // if we reached the minimum number of values we check the error on the correlation lengths
                        int min_ind = (int)(equil_cutoff * nr_xi_values);
                        // we need to calculate the average aswell as the stddev for both directions
                        double avg_xix = accumulate(xix.begin() + min_ind, xix.end(), 0.0) / (double)(nr_xi_values - min_ind);
                        double avg_xiy = accumulate(xiy.begin() + min_ind, xiy.end(), 0.0) / (double)(nr_xi_values - min_ind);

                        std::vector<double> diff_xix_total(nr_xi_values - min_ind);
                        std::vector<double> diff_xiy_total(nr_xi_values - min_ind);

                        std::transform(xix.begin() + min_ind, xix.end(), diff_xix_total.begin(), [avg_xix](double x) { return x - avg_xix; });
                        std::transform(xiy.begin() + min_ind, xiy.end(), diff_xiy_total.begin(), [avg_xiy](double x) { return x - avg_xiy; });

                        double dist_var_x = std::inner_product(diff_xix_total.begin(), diff_xix_total.end(), diff_xix_total.begin(), 0.0) / (nr_xi_values - min_ind);
                        double dist_var_y = std::inner_product(diff_xiy_total.begin(), diff_xiy_total.end(), diff_xiy_total.begin(), 0.0) / (nr_xi_values - min_ind);

                        // autocorrelation time, works also for the correlation length
                        double* xix_arr = &xix[min_ind];
                        double* xiy_arr = &xiy[min_ind];

                        double autocorr_time_x;
                        double autocorr_time_y;

                        if(nr_xi_values - min_ind > val_nr_gpu) {
                            Singleton_timer::set_startpoint("Autocorrelation time fft calculation");
                            autocorr_time_x = get_autocorrtime_fft(xix_arr, nr_xi_values - min_ind, write_interval);
                            autocorr_time_y = get_autocorrtime_fft(xiy_arr, nr_xi_values - min_ind, write_interval);
                            Singleton_timer::set_endpoint("Autocorrelation time fft calculation");
                        } else {
                            Singleton_timer::set_startpoint("Autocorrelation time calculation");
                            autocorr_time_x = get_autocorrtime(xix_arr, nr_xi_values - min_ind, write_interval);
                            autocorr_time_y = get_autocorrtime(xiy_arr, nr_xi_values - min_ind, write_interval);  // actually ds is just the write interval? which should be 1 or something like this
                            Singleton_timer::set_endpoint("Autocorrelation time calculation");
                        }

                        double xix_variance = 2 * autocorr_time_x / ((nr_xi_values-min_ind) * write_interval) * dist_var_x;
                        double xiy_variance = 2 * autocorr_time_y / ((nr_xi_values-min_ind) * write_interval) * dist_var_y;


                        double rel_stddev_xix_total = sqrt(xix_variance) / avg_xix;
                        double rel_stddev_xiy_total = sqrt(xiy_variance) / avg_xiy;

                        cout << "xix = " << avg_xix << endl;
                        cout << "xiy = " << avg_xiy << endl;
                        cout << "autocorrelation time x: " << autocorr_time_x << endl;
                        cout << "autocorrelation time y: " << autocorr_time_y << endl;
                        cout << "xix_variance = " << xix_variance << endl;
                        cout << "xiy_variance = " << xiy_variance << endl;
                        cout << "rel_stddev_total xix = " << rel_stddev_xix_total << endl;
                        cout << "rel_stddev_total xiy = " << rel_stddev_xiy_total << endl;
                        // I mean you are calculating the errors anyway, you could also write them down? But takes some time ofc, but it would make it easier to see how long a certain simulation will still take...


                        // I would say if the mean of the two relative standard deviations satisfies the condition we are fine
                        double rel_stddev_total;
                        if(observed_direction == 0) {
                            rel_stddev_total = 0.5 * (rel_stddev_xix_total + rel_stddev_xiy_total);
                        } else if(observed_direction == 1) {
                            rel_stddev_total = rel_stddev_xix_total;
                        } else if(observed_direction == 2) {
                            rel_stddev_total = rel_stddev_xiy_total;
                        }
                        cout << "observed direciton is " << observed_direction << " meaning rel_stddev_total = " << rel_stddev_total << endl;
                        if(rel_stddev_total < max_error) {
                            cout << "The system equilibrated, the equilibration lastet to t = " << t << endl;
                            cout << "xix = " << avg_xix << " +- " << rel_stddev_xix_total * avg_xix << endl;
                            cout << "xiy = " << avg_xiy << " +- " << rel_stddev_xiy_total * avg_xiy << endl;

                            // we also implement the moving factor here
                            double moving_factor_xix = getMovingFactor(nr_xi_values, min_ind, xix, avg_xix);
                            double moving_factor_xiy = getMovingFactor(nr_xi_values, min_ind, xiy, avg_xiy);
                            // lets just take the maximum...
                            double moving_factor = max(moving_factor_xix, moving_factor_xiy);
                            cout << "MOVING FACTOR = " << moving_factor << endl;

                            if(moving_factor < max_moving_factor) {
                                // we set the system to be equilibrated
                                sys.set_equilibration(t);           // For relaxation simulations this means that the simulation ends

                                // once we did this, we dont want to do that a second time?
                                equilibrated = true;
                                // the writ interval will now be the one that we destined for the quench, we just change it once here
                                write_interval = quench_t / (double)nr_values;
                            }
                        } else {
                            int minutes_remaining = get_remaining_minutes();
                            cout << "minutes remaining " << minutes_remaining << endl;
                            if(minutes_remaining < 15) {
                                // if less than 15 minutes are remaining until the jobs cancelled, we equilibrate
                                // so that the states get written and stuff
                                sys.set_equilibration(t);           // For relaxation simulations this means that the simulation ends
                                equilibrated = true;
                            }
                        }
                        if(!equilibrated) {
                            double autocorr_time = 0.5 * (autocorr_time_x + autocorr_time_y);
                            adapt_write_interval(nr_xi_values, autocorr_time);
                        } else {
                            append_parameter(filepath, "xix", avg_xix);
                            // write error
                            append_parameter(filepath, "xix_error", rel_stddev_xix_total);
                            // write autocorrelation time
                            append_parameter(filepath, "autocorrelation_time_xix", autocorr_time_x);
                            // Also write the moving factor
                            append_parameter(filepath, "moving_factor_xix", moving_factor);
                            // and the same stuff in the y direction
                            append_parameter(filepath, "xiy", avg_xiy);
                            append_parameter(filepath, "xiy_error", rel_stddev_xiy_total);
                            append_parameter(filepath, "autocorrelation_time_xiy", autocorr_time_y);
                            append_parameter(filepath, "moving_factor_xiy", moving_factor);
                        }
                    }
                }
            }
            timepoint += write_interval;
        }
    }

    void adapt_write_interval(int nr_cum_values, double autocorr_time, int values_per_autocorr_time=20) {
        // we want to adapt the write density on the autocorrelation time. But everytime we do this
        // we have to delete the values that have smaller spacing between them (the intermediate values)
        // It should be some kind of switch statement or something like this?
        // We want to write approximately 10 cumulant values during one autocorrelation time? Or maybe 100 ?
        // I think 10 should be fine
        double autocorr_time_floored = pow(10.0, floor(log10(autocorr_time)));      // we floor the autocorrelation time to the next power of 10, so 3000 would get rounded to 1000 and 101 to 100
        // we do this as we do not want to continuously change our stepping but in discrete steps
        double new_write_interval = autocorr_time_floored / (double)values_per_autocorr_time;    // if autocorr is 1000, the write interval is 100, using 100 writes per autocorr time now because ten is not satisfying somehow
        if(new_write_interval > write_interval) {
            // in this case we have to drop values accordingly to new write interval
            int keep_every_nth_val = (int)(new_write_interval / write_interval); // they should yield a glatt integer since we only work with powers of ten write interval might have been 1 and now it is 10 so we keep every 10th value
            int new_nr_cum_values = nr_cum_values / keep_every_nth_val;         // if we have 101 values with a spacing of 1, we want to keep 101, 91, ..., 1, making up for 11 values
            if (nr_cum_values % keep_every_nth_val != 0) {
                // if this is the case we can extract one more value
                new_nr_cum_values++;
            }
/*            cout << "the old number of cum values is " << nr_cum_values << endl;
            cout << "the new number of cum values is " << new_nr_cum_values << endl;*/
            // if we have 100 values with a spacing of 1, those would be the indices 0 to 99. We want to keep 99 and 89 (distance of 10 ds), 79, ..., 9 making up 10 values
            // if we have 101 values with a spacing of 1, those would have the indices 0 to 100. We would want to keep 100, 90, ..., 10, 0 so eleven values? Right now we would be missing out on 0, which would not be two bad...
            vector<double> new_xix(new_nr_cum_values);
            vector<double> new_xiy(new_nr_cum_values);
            vector<double> new_times(new_nr_cum_values);
            for(int i = 0; i < new_nr_cum_values; i++) {
                // TODO is the plus one correct? If you get memory errors look here?
                // okay so we work from behind because instead of keeping 0 - 100 we want to keep 1 - 101. But the order has to stay the same for the correlation function
                // cout << "i = " << i << "   nr_cum_values - keep_every_nth_val * i - 1 = " << nr_cum_values - keep_every_nth_val * i - 1 << "    m_vec = " << m_vec[nr_cum_values - keep_every_nth_val * i - 1] << endl;
                new_xix[new_nr_cum_values - i - 1] = xix[nr_cum_values - keep_every_nth_val * i - 1];       // - one because in a vector of length 5 the last index is 4
                new_xiy[new_nr_cum_values - i - 1] = xiy[nr_cum_values - keep_every_nth_val * i - 1];
                new_times[new_nr_cum_values - i - 1] = times[nr_cum_values - keep_every_nth_val * i - 1]; // if i = new_nr_cum_values - 1 = nr_cum_values / keep_ever_nth_val + 1 - 1 = nr_cum_values / keep .. so we are accessing times[-1] but why I dont get it
            }
            xix = new_xix;
            xiy = new_xiy;
            times = new_times;
            cout << "ADAPTED NEW WRITE INTERVAL: write_interval_old = " << write_interval;
            density = max(dt / new_write_interval, min_density);
            write_interval = dt / density;
            cout << "  write_interval_new = " << new_write_interval << endl;
            min_corr_nr = new_nr_cum_values;       // should behave as the other observer
            cout << "new minimum corr number is min_corr_nr = " << min_corr_nr << endl;
            // this is all right?
            // we decided to rewrite the file, it will save space and make the after simulation validation easier
            rewrite_file();
        }
    }

    void rewrite_file() {
        int new_nr_xi_values = xix.size();
        close_stream();
        open_stream(filepath);
        ofile << "t,xix,xiy" << endl;
        for(int i = 0; i < new_nr_xi_values; i++) {
            ofile << times[i] << "," << xix[i] << "," << xiy[i] << endl;
        }
    }

};

template <class system, class State>
class quench_combined_corr_observer: public obsver<system, State>{
    // This observer will be a mixture of the cum equilibration observer and the new density observer
    // It will calculate xi with a fixed density during equilibration phase and then the specified number
    // of xi's during the quench.
    // The density during the equilibration does not have to change, or should it? I don't think so
    // Problem could be that if the equilibration takes very long, we calculate a whole lot of xi values during the
    // equilibratioon which slows down our simulation. An adaptive density doesn't look to good when plotting and has
    // the weird artefact that we write the correlation length at different times for different simulations
    // I think we will choose a mediocre density and start with this.
    typedef obsver<system, State> obsver;
    using obsver::ofile;
    using obsver::open_stream;
    using obsver::open_app_stream;
    using obsver::close_stream;
    double write_interval = 1;
    double timepoint = 0;
    vector<double> xix{};
    vector<double> xiy{};
    vector<double> times{};
    int min_corr_nr = 5000;
    int val_nr_gpu = 1000;                 // nr of cum values at which the calculation switches to gpu
    double dt = 0.01;
    double dt_half = dt / 2.0;
    double equil_cutoff = 0.1;              // since the equilibration might influce the mean of xi
    double density = 1.0 / 100.0;            // standard density writes once every 100 steps
    double min_density = 1.0 / 10000.0;
    int observed_direction = 0;             // 0: both, 1: x, 2:y, for the error
    double max_error= 0.05;
    bool equilibrated = false;                      // for the usecase of the quench with dynamic equilibration
    int nr_values;                  // nr of values that I want to be written down during the quench
    double quench_t;                // quench time, important to calculate the write interval during the quench
    size_t Lx;                      // sizes of the system. Important to validate if the extraction was meaningful
    size_t Ly;
    double xi_cap = 0.2;            // We cap the xi size to be at most 0.2 of the corresponding system size, otherwise the extraction is unreliable and
    fs::path filepath;              // path this observer writes to
    fs::path filepath_ft;              // path this observer writes to
    double max_moving_factor = 0.02;
    int equil_counter = 0;      // This variable will count how many timepoints we already wrote down during equilibration
    int equil_nr = 100;          // nr ob observations during equilibration
    ofstream ft_file;

public:
    quench_combined_corr_observer(int nr_values) : nr_values(nr_values) {
    }
    quench_combined_corr_observer(int nr_values, int min_corr_nr, double density, double equil_cutoff):
            nr_values(nr_values), min_corr_nr(min_corr_nr), density(density), equil_cutoff(equil_cutoff) {}
    void init(fs::path path, map<Parameter, double>& paras, const system &sys) override {
        int run_nr = (int)paras[Parameter::run_nr];
        max_error = paras[Parameter::equil_error];
        if (paras[Parameter::equil_cutoff]) {
            equil_cutoff = paras[Parameter::equil_cutoff];
        }
        cout << "setting equil_cutoff to " <<  equil_cutoff << endl;
        min_corr_nr = (int)paras[Parameter::min_corr_nr];
        if (paras[Parameter::moving_factor]) {
            max_moving_factor = paras[Parameter::moving_factor];
        }

        if (paras[Parameter::corr_write_density]) {
            density = paras[Parameter::corr_write_density];
        }
        // the subsystems sizes for the xi cutting
        Lx = (size_t)paras[Parameter::subsystem_Lx];
        Ly = (size_t)paras[Parameter::subsystem_Ly];
        observed_direction = (int)paras[Parameter::observed_direction];
        // the timepoint is only zero if we do not memory initialize
        timepoint = 0.0;
        equilibrated = false;
        // we also need to reset m_vec and times, dont we?
        xix = vector<double>{};
        xiy = vector<double>{};
        times = vector<double>{};
        close_stream();
        ft_file.close();
        bool pick_up = (paras[Parameter::random_init]  == -1.0);
        if(!pick_up) {
            filepath = path / (obsver::construct_filename(run_nr) + ".corr");
            filepath_ft = path / (obsver::construct_filename(run_nr) + ".ft");
            open_stream(filepath);
            ft_file.open(filepath_ft);
            // We only want to do this if we are creating a new file?
            ofile << "t,xix,xiy" << endl;
            ft_file << "t;ft_k;ft_l" << endl;
        } else {
            // If we have the memory initialization we also want to load the correlation length values into cache
            // I thinkt we should do that before we decide to open this file also as output stream
            path += ".corr";
            filepath = path;
            filepath_ft = path.replace_extension(".ft");

            cout << "Reading correlation lengths from file" << endl;
            readXiFromFile(path, xix, xiy, times);
            timepoint = times.back();
            cout << "successful!" << endl;
            // we should adapt the write interval from the file that we read
            write_interval = times[1] - times[0];
            density = dt / write_interval;
            open_app_stream(filepath);
            ft_file.open(filepath_ft, ios::app);
        }
        cout << this->get_name() << " init called" << endl;


        // I think the starting write interval should be every one hundred steps
        dt = paras[Parameter::dt];
        dt_half = dt / 2.0;
        quench_t = sys.get_quench_time();
        cout << "the quench time is quench_t = " << quench_t << endl;
        write_interval = (1.0 / density) * dt;
    }

    string get_name() override {
        return "combined observer";
    }

    void operator()(system &sys, const State &x , double t ) override {
        // How exactly do I want this to behave? It should definitely if we are quenching write a defined number of
        // of points
        // If we are not quenching, it may vary
        if(t > timepoint - dt_half) {
            double xix_val, xiy_val;
            Singleton_timer::set_startpoint("Xi Calculation");
            thrust::host_vector<double> ft_squared_k;
            thrust::host_vector<double> ft_squared_l;
            sys.get_ft_vectors(x, ft_squared_k, ft_squared_l);
            sys.fit_ft(ft_squared_k, ft_squared_l, xix_val, xiy_val);
            Singleton_timer::set_endpoint("Xi Calculation");
            // The logic of the cutting should take place here or only for calculating the mean? I think we also cut it
            // but we might change that later
            xix_val = min(xix_val, xi_cap * Lx);
            xiy_val = min(xiy_val, xi_cap * Ly);
            ofile << t << "," << xix_val << "," << xiy_val << endl;

            ft_file << t << ";";
            for(int i = 0; i < Lx; i++) {
                if(i == 0) {
                    ft_file << ft_squared_k[i];
                } else {
                    ft_file << "," << ft_squared_k[i];
                }
            }
            ft_file << ";";
            for(int j = 0; j < Ly; j++) {
                if(j == 0) {
                    ft_file << ft_squared_l[j];
                } else {
                    ft_file << "," << ft_squared_l[j];
                }
            }
            ft_file << endl;
            equil_counter++;

            double end_quench_t = sys.get_end_quench_time();
            double s_eq_t = end_quench_t - quench_t;       // i think this returns infinity if we are not equilibrated?
            // shouldnt it be as easy as this?
            if (t > s_eq_t) {
                if (t > end_quench_t) {
                    write_interval = (1.0 / density) * dt;
                } else {
                    write_interval = quench_t / ((double)nr_values);
                }
            }
            else {
                xix.push_back(xix_val);
                xiy.push_back(xiy_val);
                times.push_back(t);
                int nr_xi_values = xix.size();      // check how many corr values we already have, the number of xiy values equals the number of xix values
                // we lack the complete logic to change the stepsize since we said we use a constant density
                if(nr_xi_values > min_corr_nr){
                    int error_every_n_steps = 10;
                    if(xix.size() % error_every_n_steps == 0) {
                        int min_ind;
                        double avg_xix;
                        double avg_xiy;
                        double autocorr_time_x;
                        double autocorr_time_y;
                        double xix_variance;
                        double xiy_variance;
                        double rel_stddev_xix_total;
                        double rel_stddev_xiy_total;
                        evaluate(nr_xi_values, min_ind, avg_xix, avg_xiy, autocorr_time_x, autocorr_time_y,
                                 xix_variance, xiy_variance,
                                 rel_stddev_xix_total, rel_stddev_xiy_total);


                        cout << "xix = " << avg_xix << endl;
                        cout << "xiy = " << avg_xiy << endl;
                        cout << "autocorrelation time x: " << autocorr_time_x << endl;
                        cout << "autocorrelation time y: " << autocorr_time_y << endl;
                        cout << "xix_variance = " << xix_variance << endl;
                        cout << "xiy_variance = " << xiy_variance << endl;
                        cout << "rel_stddev_total xix = " << rel_stddev_xix_total << endl;
                        cout << "rel_stddev_total xiy = " << rel_stddev_xiy_total << endl;
                        // I mean you are calculating the errors anyway, you could also write them down? But takes some time ofc, but it would make it easier to see how long a certain simulation will still take...


                        // I would say if the mean of the two relative standard deviations satisfies the condition we are fine
                        double rel_stddev_total;
                        if(observed_direction == 0) {
                            rel_stddev_total = 0.5 * (rel_stddev_xix_total + rel_stddev_xiy_total);
                        } else if(observed_direction == 1) {
                            rel_stddev_total = rel_stddev_xix_total;
                        } else if(observed_direction == 2) {
                            rel_stddev_total = rel_stddev_xiy_total;
                        }
                        if(rel_stddev_total < max_error) {
                            cout << "xix = " << avg_xix << " +- " << rel_stddev_xix_total * avg_xix << endl;
                            cout << "xiy = " << avg_xiy << " +- " << rel_stddev_xiy_total * avg_xiy << endl;
                            cout << "min ind for moving factor " << min_ind << endl;
                            // we also implement the moving factor here
                            double moving_factor_xix = getMovingFactor(nr_xi_values, min_ind, xix, avg_xix);
                            double moving_factor_xiy = getMovingFactor(nr_xi_values, min_ind, xiy, avg_xiy);
                            // lets just take the maximum...
                            double moving_factor = max(moving_factor_xix, moving_factor_xiy);

                            if(moving_factor < max_moving_factor) {
                                // we set the system to be equilibrated
                                cout << "The system equilibrated, the equilibration lastet to t = " << t << endl;

                                sys.set_equilibration(t);           // For relaxation simulations this means that the simulation ends

                                append_parameter(filepath, "xix", avg_xix);
                                // write error
                                append_parameter(filepath, "xix_error", rel_stddev_xix_total);
                                // write autocorrelation time
                                append_parameter(filepath, "autocorrelation_time_xix", autocorr_time_x);
                                // Also write the moving factor
                                append_parameter(filepath, "moving_factor_xix", moving_factor);
                                // and the same stuff in the y direction
                                append_parameter(filepath, "xiy", avg_xiy);
                                append_parameter(filepath, "xiy_error", rel_stddev_xiy_total);
                                append_parameter(filepath, "autocorrelation_time_xiy", autocorr_time_y);
                                append_parameter(filepath, "moving_factor_xiy", moving_factor);
                                // once we did this, we dont want to do that a second time?
                                equilibrated = true;
                                // the writ interval will now be the one that we destined for the quench, we just change it once here
                                write_interval = quench_t / (double)nr_values;
                            }
                            cout << "MOVING FACTOR = " << moving_factor << " > " << max_moving_factor  << endl;
                        }
                        cout << "observed direciton is " << observed_direction << " meaning rel_stddev_total = " << rel_stddev_total;
                        cout << " > " << max_error << endl;

/*                        if(!equilibrated) {
                            cout << "equil_counter = " << equil_counter << endl;
                            if(equil_counter >= 2 * equil_nr) {
                                adapt_write_interval();
                                equil_counter = equil_nr;
                            }
                        }*/
                    }
                }
            }
            timepoint += write_interval;
        }
    }

    void evaluate(int nr_xi_values, int &min_ind, double &avg_xix, double &avg_xiy, double &autocorr_time_x,
                  double &autocorr_time_y, double &xix_variance, double &xiy_variance, double &rel_stddev_xix_total,
                  double &rel_stddev_xiy_total) {
        min_ind = (int)(equil_cutoff * nr_xi_values);
        // we need to calculate the average aswell as the stddev for both directions
        avg_xix = accumulate(xix.begin() + min_ind, xix.end(), 0.0) / (double)(nr_xi_values - min_ind);
        avg_xiy = accumulate(xiy.begin() + min_ind, xiy.end(), 0.0) / (double)(nr_xi_values - min_ind);

        std::vector<double> diff_xix_total(nr_xi_values - min_ind);
        std::vector<double> diff_xiy_total(nr_xi_values - min_ind);

        std::transform(xix.begin() + min_ind, xix.end(), diff_xix_total.begin(), [avg_xix](double x) { return x - avg_xix; });
        std::transform(xiy.begin() + min_ind, xiy.end(), diff_xiy_total.begin(), [avg_xiy](double x) { return x - avg_xiy; });

        double dist_var_x = std::inner_product(diff_xix_total.begin(), diff_xix_total.end(), diff_xix_total.begin(), 0.0) / (nr_xi_values - min_ind);
        double dist_var_y = std::inner_product(diff_xiy_total.begin(), diff_xiy_total.end(), diff_xiy_total.begin(), 0.0) / (nr_xi_values - min_ind);

        // autocorrelation time, works also for the correlation length
        double* xix_arr = &xix[min_ind];
        double* xiy_arr = &xiy[min_ind];

        if(nr_xi_values - min_ind > val_nr_gpu) {
            Singleton_timer::set_startpoint("Autocorrelation time fft calculation");
            autocorr_time_x = get_autocorrtime_fft(xix_arr, nr_xi_values - min_ind, write_interval);
            autocorr_time_y = get_autocorrtime_fft(xiy_arr, nr_xi_values - min_ind, write_interval);
            Singleton_timer::set_endpoint("Autocorrelation time fft calculation");
        } else {
            Singleton_timer::set_startpoint("Autocorrelation time calculation");
            autocorr_time_x = get_autocorrtime(xix_arr, nr_xi_values - min_ind, write_interval);
            autocorr_time_y = get_autocorrtime(xiy_arr, nr_xi_values - min_ind, write_interval);  // actually ds is just the write interval? which should be 1 or something like this
            Singleton_timer::set_endpoint("Autocorrelation time calculation");
        }

        xix_variance = 2 * autocorr_time_x / ((nr_xi_values-min_ind) * write_interval) * dist_var_x;
        xiy_variance = 2 * autocorr_time_y / ((nr_xi_values-min_ind) * write_interval) * dist_var_y;

        rel_stddev_xix_total = sqrt(xix_variance) / avg_xix;
        rel_stddev_xiy_total = sqrt(xiy_variance) / avg_xiy;
    }

    void adapt_write_interval() {
        int nr_values = xix.size(); // this should be 200 but who knows amk
        double new_write_interval = 2 * write_interval;    // if autocorr is 1000, the write interval is 100, using 100 writes per autocorr time now because ten is not satisfying somehow
            // in this case we have to drop values accordingly to new write interval
        int keep_every_nth_val = 2; // they should yield a glatt integer since we only work with powers of ten write interval might have been 1 and now it is 10 so we keep every 10th value
        int new_nr_cum_values = nr_values / keep_every_nth_val;         // if we have 101 values with a spacing of 1, we want to keep 101, 91, ..., 1, making up for 11 values
        if (nr_values % keep_every_nth_val != 0) {
            // if this is the case we can extract one more value
            new_nr_cum_values++;
        }
/*            cout << "the old number of cum values is " << nr_values << endl;
        cout << "the new number of cum values is " << new_nr_cum_values << endl;*/
        // if we have 100 values with a spacing of 1, those would be the indices 0 to 99. We want to keep 99 and 89 (distance of 10 ds), 79, ..., 9 making up 10 values
        // if we have 101 values with a spacing of 1, those would have the indices 0 to 100. We would want to keep 100, 90, ..., 10, 0 so eleven values? Right now we would be missing out on 0, which would not be two bad...
        vector<double> new_xix(new_nr_cum_values);
        vector<double> new_xiy(new_nr_cum_values);
        vector<double> new_times(new_nr_cum_values);
        for(int i = 0; i < new_nr_cum_values; i++) {
            // TODO is the plus one correct? If you get memory errors look here?
            // okay so we work from behind because instead of keeping 0 - 100 we want to keep 1 - 101. But the order has to stay the same for the correlation function
            // cout << "i = " << i << "   nr_values - keep_every_nth_val * i - 1 = " << nr_values - keep_every_nth_val * i - 1 << "    m_vec = " << m_vec[nr_values - keep_every_nth_val * i - 1] << endl;
            new_xix[i] = xix[keep_every_nth_val * i];       // - one because in a vector of length 5 the last index is 4
            new_xiy[i] = xiy[keep_every_nth_val * i];
            new_times[i] = times[keep_every_nth_val * i]; // if i = new_nr_cum_values - 1 = nr_values / keep_ever_nth_val + 1 - 1 = nr_values / keep .. so we are accessing times[-1] but why I dont get it
        }
        xix = new_xix;
        xiy = new_xiy;
        times = new_times;
        cout << "ADAPTED NEW WRITE INTERVAL: write_interval_old = " << write_interval;
        density = max(dt / new_write_interval, min_density);
        write_interval = dt / density;
        cout << "  write_interval_new = " << new_write_interval << endl;
        min_corr_nr = new_nr_cum_values;       // should behave as the other observer
        cout << "new minimum corr number is min_corr_nr = " << min_corr_nr << endl;
        // this is all right?
        // we decided to rewrite the file, it will save space and make the after simulation validation easier
        rewrite_xi_file();
        rewrite_ft_file();
    }

    void rewrite_xi_file() {
        int new_nr_xi_values = xix.size();
        close_stream();
        open_stream(filepath);
        ofile << "t,xix,xiy" << endl;
        for(int i = 0; i < new_nr_xi_values; i++) {
            ofile << times[i] << "," << xix[i] << "," << xiy[i] << endl;
        }
    }

    void rewrite_ft_file() {
        cout << "REWRITING FT FILE?" << endl;
        ft_file.close();
        // we know have to read everything?
        omitHalfTimes(filepath_ft);
        ft_file.open(filepath_ft, ios::app);
    }
};

template <class system, class State>
class ft_observer : public obsver<system, State>{
    // Okay the observing pattern could be totally wild, i probably somehow have to initialize the observer
    // outside of the simulation class We definetely need its own constructor here
    typedef obsver<system, State> obsver;
    using obsver::ofile;
    using obsver::open_stream;
    using obsver::open_app_stream;
    using obsver::close_stream;
    int nr_values;
    double write_interval = 1;
    double timepoint = 0;
    size_t Lx;
    size_t Ly;
public:
    ft_observer(int nr_values) : nr_values(nr_values) {
    }
    void init(fs::path path, map<Parameter, double>& paras, const system &sys) override {
        int run_nr = (int)paras[Parameter::run_nr];
        timepoint = 0.0;

        close_stream();
        bool pick_up = (paras[Parameter::random_init]  == -1.0);
        fs::path folderpath;
        if(!pick_up) {
            open_stream(path / (obsver::construct_filename(run_nr) + ".ft"));
            ofile << "t;ft_k;ft_l" << endl;
        } else {
            path += ".ft";
            open_app_stream(path);
        }
        cout << "ft observer init called" << endl;

        // I think this will have less performance impact than an if statement catching the first observer operation
        // Make sure to use this observer only with systems that have a get_end_T method
        double end_t = paras[end_time];
        Lx = paras[subsystem_Lx];
        Ly = paras[subsystem_Ly];
        write_interval = end_t / (double)nr_values;
    }

    string get_name() override {
        return "ft observer";
    }

    void operator()(system &sys, const State &x , double t ) override {
        if(t > timepoint) {
            double *ft_k, *ft_l;
            ft_k = new double[Lx];
            ft_l = new double[Ly];

            sys.calc_ft(x, ft_k, ft_l);
            ofile << t << ";";
            for(int i = 0; i < Lx; i++) {
                if(i == 0) {
                    ofile << ft_k[i];
                } else {
                    ofile << "," << ft_k[i];
                }
            }
            ofile << ";";
            for(int j = 0; j < Ly; j++) {
                if(j == 0) {
                    ofile << ft_l[j];
                } else {
                    ofile << "," << ft_l[j];
                }
            }
            ofile << endl;
            timepoint += write_interval;
        }
    }
};

template <class system, class State>
class ft_observer_density : public obsver<system, State>{
    // Okay the observing pattern could be totally wild, i probably somehow have to initialize the observer
    // outside of the simulation class We definetely need its own constructor here
    typedef obsver<system, State> obsver;
    using obsver::ofile;
    using obsver::open_stream;
    using obsver::open_app_stream;
    using obsver::close_stream;
    int nr_values;
    double write_interval;
    double write_density;
    double timepoint = 0;
    double dt;
    size_t Lx;
    size_t Ly;
public:
    ft_observer_density() {
    }
    void init(fs::path path, map<Parameter, double>& paras, const system &sys) override {
        int run_nr = (int)paras[Parameter::run_nr];
        write_density = paras[Parameter::ft_write_density];
        dt = paras[Parameter::dt];
        timepoint = 0.0;

        close_stream();
        bool pick_up = (paras[Parameter::random_init]  == -1.0);
        fs::path folderpath;
        if(!pick_up) {
            open_stream(path / (obsver::construct_filename(run_nr) + ".ft"));
            ofile << "t;ft_k;ft_l" << endl;
        } else {
            path += ".ft";
            open_app_stream(path);
        }
        cout << "ft observer init called" << endl;

        // I think this will have less performance impact than an if statement catching the first observer operation
        // Make sure to use this observer only with systems that have a get_end_T method
        double end_t = paras[end_time];
        Lx = paras[subsystem_Lx];
        Ly = paras[subsystem_Ly];
        write_interval = dt / write_density;
    }

    string get_name() override {
        return "ft density observer";
    }

    void operator()(system &sys, const State &x , double t ) override {
        if(t > timepoint - 0.5 * dt) {
            double *ft_k, *ft_l;
            ft_k = new double[Lx];
            ft_l = new double[Ly];

            sys.calc_ft(x, ft_k, ft_l);
            ofile << t << ";";
            for(int i = 0; i < Lx; i++) {
                if(i == 0) {
                    ofile << ft_k[i];
                } else {
                    ofile << "," << ft_k[i];
                }
            }
            ofile << ";";
            for(int j = 0; j < Ly; j++) {
                if(j == 0) {
                    ofile << ft_l[j];
                } else {
                    ofile << "," << ft_l[j];
                }
            }
            ofile << endl;
            timepoint += write_interval;
        }
    }
};

template <class system, class State>
class quench_ft_observer : public obsver<system, State>{
    // Okay the observing pattern could be totally wild, i probably somehow have to initialize the observer
    // outside of the simulation class We definetely need its own constructor here
    typedef obsver<system, State> obsver;
    using obsver::ofile;
    using obsver::open_stream;
    using obsver::open_app_stream;
    using obsver::close_stream;
    int nr_values;
    double write_interval = 1;
    double timepoint = 0;
    double end_quench_t;
    double quench_t;
    double s_eq_t;
    size_t Lx;
    size_t Ly;
    int equil_counter = 0;      // This variable will count how many timepoints we already wrote down during equilibration
    int equil_nr = 100;          // nr ob observations during equilibration
public:
    quench_ft_observer(int nr_values) : nr_values(nr_values) {
    }
    void init(fs::path path, map<Parameter, double>& paras, const system &sys) override {
        int run_nr = (int)paras[Parameter::run_nr];
        timepoint = 0.0;

        close_stream();
        bool pick_up = (paras[Parameter::random_init]  == -1.0);
        if(pick_up) {
            open_stream(path / (obsver::construct_filename(run_nr) + ".ft"));
            ofile << "t;ft_k;ft_l" << endl;
        } else {
            path += ".ft";
            open_app_stream(path);
        }
        cout << "quench ft observer init called" << endl;

        // I think this will have less performance impact than an if statement catching the first observer operation
        // Make sure to use this observer only with systems that have a get_end_T method
        double end_t = paras[end_time];
        Lx = paras[subsystem_Lx];
        Ly = paras[subsystem_Ly];

        end_quench_t = sys.get_end_quench_time();
        quench_t = sys.get_quench_time();
        s_eq_t = end_quench_t - quench_t;
        cout << "end_quench_t " << end_quench_t << "  quench_t " << quench_t << "  s_eq_t " << s_eq_t << endl;
        // the starting write interval is chosen so that 10% of the nr_values lie in the starting equilibration
        write_interval = s_eq_t / ((double)nr_values * 0.1);
    }

    string get_name() override {
        return "quench ft observer";
    }

    void operator()(system &sys, const State &x , double t ) override {
        if(t > timepoint) {
            double *ft_k, *ft_l;
            ft_k = new double[Lx];
            ft_l = new double[Ly];

            sys.calc_ft(x, ft_k, ft_l);
            ofile << t << ";";
            for(int i = 0; i < Lx; i++) {
                if(i == 0) {
                    ofile << ft_k[i];
                } else {
                    ofile << "," << ft_k[i];
                }
            }
            ofile << ";";
            for(int j = 0; j < Ly; j++) {
                if(j == 0) {
                    ofile << ft_l[j];
                } else {
                    ofile << "," << ft_l[j];
                }
            }
            ofile << endl;

            // here we adjust the write interval depending on where we are?
            if (timepoint < s_eq_t) {
                write_interval = s_eq_t / ((double)nr_values * 0.1);
                timepoint += write_interval;
                timepoint = min(timepoint, s_eq_t);
                equil_counter++;
                if(equil_counter > equil_nr) {
                    // if this is the case we want to rewrite the file to limit the data generation. Also we will increase
                    // the write interval since it seems that it will take some more time
                    // I hope that doesnt clash with the recordings of the corr obs. To be honest we should directly combine them
                }

            } else if ((end_quench_t > timepoint) & (timepoint >= s_eq_t)) {
                // case that we are during the quench, we want 80 % of values to lie here
                write_interval = quench_t / ((double)nr_values * 0.8);
                timepoint += write_interval;
                timepoint = min(timepoint, end_quench_t);
            } else if (timepoint >= end_quench_t) {
                // we are in equilibration phase after quench, again 10% of values
                write_interval = s_eq_t / ((double)nr_values * 0.1);
                timepoint += write_interval;
            }
            cout << "write interval = " << write_interval << endl;
        }
    }
};

template <class system, class State>
class quench_equilibration_observer : public standard_observer<system, State> {
    // this observer is supposed to write down some meshs even for the dynamically equilibrating quench
    double quench_t;        // the time that the quench takes
    bool startpoint = true; // whether we are at the starting point whcih we want to write down
    double dt_half;         // we adopt this variable to f*ing have better timepoints?
    double density = 1.0/100.0;  // should have the same density as the corr equilibration observer
public:
    typedef standard_observer<system, State> standard_observer;
    using standard_observer::write_interval;
    using standard_observer::nr_values;
    using standard_observer::timepoint;
    using standard_observer::ofile;

    quench_equilibration_observer(int nr_values): standard_observer(nr_values) {
        // should this be the observer that only writes a few values for the equilibrium process?
        // I am not sure at the moment, I think i am just going to write the simplest one for now, so
        // just equidistant saving values.
        // We try to initialize the observer outside the class since it probably needs
        // multiple unique-to-observer parameters
        // problem is it still might depend on system parameters. and i don't want another if for every step
        // even though i am not really sure if it impacts the performance AT ALL or if the compiler knows better
    }
    void init(fs::path folderpath, map<Parameter, double>& paras, const system &sys) override {
        cout <<  this->get_name() << " init is called" << endl;
        // the write interval is kind of tricky, since we don't know when the corr_observer will
        // determine that we are equilibrated, atm the corr equilibration observer observes every 100 steps and this is
        // not variable, for starters we can do that here to.
        double dt = paras[Parameter::dt];
        dt_half = dt/2.0;
        if(paras[corr_write_density]) {
            density = paras[corr_write_density];
        }
        write_interval = dt / density;      // We probably should make this density variable
        startpoint = true;              // reset the startpoint variable
        // set the quench_t
        quench_t = sys.get_quench_time();

        cout << "Write interval is " << write_interval << endl;
        standard_observer::init(folderpath, paras, sys);
    }

    string get_name() override {
        return "quench equilibration observer";
    }

    void operator()(system &sys, const State &x , double t ) override {
        // the observing will also be a bit different than before, but not much
        if(t > timepoint - dt_half) {
            // we definetely dont write every timepoint
            // we check if we se the s_eq_t of the quench system to a time that is smaller or equal to the current time
            // we have to get the equilibration time with a detour
            double eq_t = sys.get_end_quench_time() - quench_t;     // this should be the time that the quench ends, so eq_t + quench_t - the time that the quench takes, so quench_t
            // I guess I want to know the  now?
            double h = sys.get_h();
            //cout << "Current h = " << h << endl;
            double T = sys.get_cur_T();
            // cout << "Current T = " << T << endl;
            // cout << "h/T = " << h / T << endl;
            // write
            if(startpoint){
                this->write_state(sys, x, t);
                startpoint = false;
            } else if(t > sys.get_end_t() - dt_half) {
                this->write_state(sys, x, t);
            }
            else if((t > eq_t - dt_half) & (t < sys.get_end_quench_time() + dt_half)) {
                // If the current time is larger or equal the time we need to equilibrate, we write down and change the write interval
                this->write_state(sys, x, t);

                // advance timeoint
                write_interval = quench_t / nr_values; // assume the quench takes t = 100 and we want to write two values, we want to write them in the middle and in the end
                // so we set the write interval to 50 and look at t=50 and t = 100. We have to make sure that we catch the last value
            }
            timepoint += write_interval;
        }
    }

};

template <class system, class State>
class density_quench_ft_observer : public obsver<system, State>{
    // Okay the observing pattern could be totally wild, i probably somehow have to initialize the observer
    // outside of the simulation class We definetely need its own constructor here
    typedef obsver<system, State> obsver;
    using obsver::ofile;
    using obsver::open_stream;
    using obsver::open_app_stream;
    using obsver::close_stream;
    int nr_values;
    double write_interval = 1;
    double timepoint = 0;
    double end_quench_t;
    double quench_t;
    double s_eq_t;
    double dt;
    double standard_density = 1.0/100.0;    // the standard density will be one value every 100 steps, but adjustable?
    size_t Lx;
    size_t Ly;
public:
    density_quench_ft_observer(int nr_values) : nr_values(nr_values) {
    }
    void init(fs::path path, map<Parameter, double>& paras, const system &sys) override {
        int run_nr = (int)paras[Parameter::run_nr];
        timepoint = 0.0;

        close_stream();
        bool pick_up = (paras[Parameter::random_init]  == -1.0);
        if(!pick_up) {
            open_stream(path / (obsver::construct_filename(run_nr) + ".ft"));
            ofile << "t;ft_k;ft_l" << endl;
        } else {
            path += ".ft";
            open_app_stream(path);
        }
        cout << "density quench ft observer init called" << endl;

        // I think this will have less performance impact than an if statement catching the first observer operation
        // Make sure to use this observer only with systems that have a get_end_T method
        double end_t = paras[end_time];
        Lx = paras[subsystem_Lx];
        Ly = paras[subsystem_Ly];

        // the plan is to have nr_values ft during the quench and the same density during equilibration.
        // no wait this is again dumb for short and long quenches.
        // damn
        dt = paras[Parameter::dt];
        quench_t = sys.get_quench_time();
        write_interval = (1.0 / standard_density) * dt;
    }

    string get_name() override {
        return "density quench ft observer";
    }

    void operator()(system &sys, const State &x , double t ) override {
        if(t > timepoint) {
            double *ft_k, *ft_l;
            ft_k = new double[Lx];
            ft_l = new double[Ly];

            sys.calc_ft(x, ft_k, ft_l);
            ofile << t << ";";
            for(int i = 0; i < Lx; i++) {
                if(i == 0) {
                    ofile << ft_k[i];
                } else {
                    ofile << "," << ft_k[i];
                }
            }
            ofile << ";";
            for(int j = 0; j < Ly; j++) {
                if(j == 0) {
                    ofile << ft_l[j];
                } else {
                    ofile << "," << ft_l[j];
                }
            }
            ofile << endl;

            // here we adjust the write interval depending on where we are?
            double s_eq_t = sys.get_end_quench_time() - quench_t;
            // shouldnt it be as easy as this?
            if (t > s_eq_t) {
                write_interval = quench_t / ((double)nr_values);
            }
            timepoint += write_interval;
            // was that all, should I test my abomination?
        }
    }
};
#endif //CUDAPROJECT_OBSERVERS_CUH
