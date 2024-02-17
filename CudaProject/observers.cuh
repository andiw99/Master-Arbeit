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

// We need those observers that do stuff like writing the values to a file
class observer {
public:
    ofstream ofile;
    template<class State, class System>
    void operator()(const System sys, const State &x , double t ) {
        cout << "Base Observer gets called" << endl;
    }

    observer() {
    }

    template<class State, class System>
    void write(System &sys, const State &x , double t ) {

    }

    void open_stream(fs::path filepath) {
        ofile.open(filepath);
    }

    void close_stream(fs::path filepath) {
        ofile.close();
    }
    template <class System>
    void init(fs::path folderpath, map<Parameter, double>& paras, const System &sys) {
        cout << "dummy init is called" << endl;
    }

    string get_name() {
        return "base observer";
    }
};

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
    using obsver::close_stream;
    standard_observer(int nr_values): nr_values(nr_values) {
        // this one has just the logic that we write nr_values equidistant points, applies to quench aswell as to
        // Relaxation. The init is a bit different since we need to find out the write interval
        // It also writes the run parameters to a file
    }
    virtual void init(fs::path folderpath, map<Parameter, double>& paras, const system &sys)  {
        // I think this will have less performance impact than an if statement catching the first observer operation
        // Make sure to use this observer only with systems that have a get_end_T method

        // open the file to write the info to, in this case it will be just run_nr.csv
        timepoint = 0.0;
        // I think we will add the run number to the paras of the run
        int run_nr = (int)paras[Parameter::run_nr];
        // we just write the parameters first
        close_stream();
        string filename = obsver::construct_filename(run_nr);
        open_stream(folderpath / (filename + ".txt"));
        write_parameters(ofile, paras);
        ofile << "system" << "," << sys.get_name() << endl;
        // dont forget to close!;
        close_stream();

        open_stream(folderpath / (filename + ".csv"));
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
        cout << "equilibration obs init is called, write_interval is  " << endl;
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

/*template <class system, class State>
class runtime_observer : public obsver<system, State> {
    // supposed to write the time the run took
    // I am thinking about adding a 'conclude'-function to the observers since i dont really know how
    // to efficiently read and write the runtime otherwise
    // or we just don't care about performance since we can turn this observer off if we do real runs...
    double end_t;
public:
    typedef obsver<system, State> obsver;
    using obsver::ofile;
    using obsver::open_stream;
    using obsver::close_stream;
    fs::path path = "";
    void operator()(system &sys, const State &x , double t) override {
        if(t >= end_t) {
            int duration = sys.timer.get_elapsed_time();
            ofile << duration << ",";
        }
    }

    void init(fs::path folderpath, map<Parameter, double>& paras, const system &sys) override {
        // this observer is supposed to work with a finite t
        end_t = paras[end_time];
        // we open the runtimes file if we are doing a new measurment, so if the folderpath chagnes
        // and only if the file isnt open already
        cout << "stream is open: " << ofile.is_open() << endl;
        if(path != folderpath || !ofile.is_open()) {
            path = folderpath;
            close_stream();
            cout << "opening new stream" << endl;
            open_stream(folderpath/ "runtimes");
        }
    }

    ~runtime_observer() {
        cout << "deleting runtime observer " << endl;
        close_stream();
    }

    string get_name() override {
        return "runtime observer";
    }
};*/

template <class system, class State>
class NER_observer : public obsver<system, State>{
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
    NER_observer(int nr_values) : nr_values(nr_values) {
    }
    void init(fs::path folderpath, map<Parameter, double>& paras, const system &sys) override {
        int run_nr = (int)paras[Parameter::run_nr];
        timepoint = 0.0;

        close_stream();
        open_stream(folderpath / (to_string(run_nr) + "-ner"));
        cout << "NER init is called" << endl;
        ofile << "t,m,f_mm,f_me" << endl;

        // I think this will have less performance impact than an if statement catching the first observer operation
        // Make sure to use this observer only with systems that have a get_end_T method
        double end_t = paras[end_time];
        write_interval = end_t / (double)nr_values;
    }

    string get_name() override {
        return "NER observer";
    }

    void operator()(system &sys, const State &x , double t ) override {
        if(t > timepoint) {
            ofile << t << "," << sys.calc_m(x) << "," << sys.calc_f_mm(x) << ",";
            ofile << sys.calc_f_me(x) << endl;
            timepoint += write_interval;
        }
    }

    };

template <class system, class State>
class cum_observer : public obsver<system, State>{
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
    cum_observer(int nr_values) : nr_values(nr_values) {
    }
    void init(fs::path folderpath, map<Parameter, double>& paras, const system &sys) override {
        int run_nr = (int)paras[Parameter::run_nr];
        timepoint = 0.0;

        close_stream();
        open_stream(folderpath / (obsver::construct_filename(run_nr) + ".cum"));
        cout << this->get_name() << " init called" << endl;
        ofile << "t,U_L" << endl;

        // I think this will have less performance impact than an if statement catching the first observer operation
        // Make sure to use this observer only with systems that have a get_end_T method
        double end_t = paras[end_time];
        write_interval = end_t / (double)nr_values;
    }

    string get_name() override {
        return "cum observer";
    }

    void operator()(system &sys, const State &x , double t ) override {
        if(t > timepoint) {
            ofile << t << "," << sys.calc_binder(x) << endl;
            timepoint += write_interval;
        }
    }

};

template <class system, class State>
class cum_equilibration_observer: public obsver<system, State>{
    // Okay the observing pattern could be totally wild, i probably somehow have to initialize the observer
    // outside of the simulation class We definetely need its own constructor here
    typedef obsver<system, State> obsver;
    using obsver::ofile;
    using obsver::open_stream;
    using obsver::close_stream;
    double write_interval = 1;
    double timepoint = 0;
    vector<double> U_L{};
    vector<double> times{};
    int min_cum_nr = 500;
    int cum_nr_gpu = 25000;                 // nr of cum values at which the calculation switches to gpu
    double write_density = 1.0 / 100.0;      // we increase the write density to get a more useful value for the autocorrelation function?
    double dt = 0.01;
    double dt_half;
    double equil_cutoff = 0.1;              // since the equilibration might influce the mean of U_L a lot we cut a certain portion of U_L values
    double max_error= 0.001;
    int cum_nr = 0;                         // current number in the averageing process
    bool equilibrated = false;                      // for the usecase of the quench with dynamic equilibration
public:
    cum_equilibration_observer(int min_cum_nr) : min_cum_nr(min_cum_nr) {
    }

    cum_equilibration_observer(int min_cum_nr, double write_density) : min_cum_nr(min_cum_nr), write_density(write_density) {
    }

    cum_equilibration_observer(int min_cum_nr, double write_density, double equil_cutoff) : min_cum_nr(min_cum_nr), write_density(write_density), equil_cutoff(equil_cutoff) {
    }

    void init(fs::path folderpath, map<Parameter, double>& paras, const system &sys) override {
        int run_nr = (int)paras[Parameter::run_nr];
        max_error = paras[Parameter::equil_error];
        if (paras[Parameter::min_cum_nr]){
            min_cum_nr = (int)paras[Parameter::min_cum_nr];
        }
        timepoint = 0.0;
        equilibrated = false;
        // we also need to reset U_L and times, dont we?
        U_L = vector<double>{};
        times = vector<double>{};
        close_stream();
        open_stream(folderpath / (obsver::construct_filename(run_nr) + ".cum"));
        cout << this->get_name() << " init called" << endl;
        ofile << "t,U_L" << endl;

        // I think the starting write interval should be every one hundred steps
        dt = paras[Parameter::dt];
        dt_half = dt / 2.0;
        write_interval = dt / write_density;
    }

    string get_name() override {
        return "cum equilibration observer";
    }

    void operator()(system &sys, const State &x , double t ) override {
        if(t > timepoint - dt_half) {
            // advancing the cum nr
            cum_nr++;
            // we calculate the cumulant and write it down
            double cum = sys.calc_binder(x);
            ofile << t << "," << cum << endl;
            // add the cumulant and the times to the vectors to keep track
            U_L.push_back(cum);
            times.push_back(t);
            // now we want to see if we have to adjust the write interval
            // we have to make sure that we got 5 fresh cum values
            if(!equilibrated) {
                int nr_cum_values = U_L.size();      // check how many cum values we already have
                // now use the last avg_nr of cum values to calculate a mean U_L
                // use transform reduce?
                // does it work like this? U_L.end() - avg_nr is the n-th last value in the vector?
                if(nr_cum_values >= min_cum_nr){
                    int error_every_n_steps = (int)(10000.0 * write_density);
                    if(U_L.size() % error_every_n_steps == 0) {
                        // We dont need to calculate the error and stuff everytime we write down a U_L
                        // IF the density is small, we should write down more often
                        // the equilibration phase might influence the mean a lot, should we cut of the first x% of the values?
                        int min_ind = (int)(equil_cutoff * nr_cum_values);
                        // again calculate mean and stddev. We want the standarddeviation of the mean value this time?
                        double avg_U_L = accumulate(U_L.begin() + min_ind, U_L.end(), 0.0) / (double)(nr_cum_values - min_ind);

                        std::vector<double> diff_total(nr_cum_values - min_ind);

                        std::transform(U_L.begin() + min_ind, U_L.end(), diff_total.begin(), [avg_U_L](double x) { return x - avg_U_L; });

                        double dist_var = std::inner_product(diff_total.begin(), diff_total.end(), diff_total.begin(), 0.0) / (nr_cum_values-min_ind);
                        // this is the variance of the distribution (in wrong because the values are correlated), so we
                        // need the autocorrelation time to adjust for this. The function works with arrays
                        double* U_L_arr = &U_L[min_ind];
                        double autocorr_time;

                        if(nr_cum_values - min_ind > cum_nr_gpu) {
                            autocorr_time = get_autocorrtime_gpu(U_L_arr, nr_cum_values - min_ind, write_interval);
                        } else {
                            autocorr_time = get_autocorrtime(U_L_arr, nr_cum_values - min_ind, write_interval);  // actually ds is just the write interval? which should be 1 or something like this
                        }
                        // (nr_cum_values-min_ind) * the write interval is the total time of the part of the simulation that we are considering
                        // sq sum total should be the thing that we called the standard deviation of the distribution
                        double U_L_variance = 2 * autocorr_time / ((nr_cum_values-min_ind) * write_interval) * dist_var;
                        double rel_stddev_total = sqrt(U_L_variance) / avg_U_L;

                        cout << "autocorrelation time: " << autocorr_time << endl;
                        cout << "U_L variance: " << U_L_variance << endl;
                        cout << "rel_stddev_total = " << rel_stddev_total << endl;


                        // Okay we found out that this is not the real error, so we need a function that calculates the error
                        // or first a function that calculates the integrated autocorrleation time
                        // We indded have two auto correlation times for the two directions. Could they be individually sized?
                        // the question is now if we want to extract avg_U_L and its error if we are equilibrated
                        // we definitely should? But how and where? somehow into the parameter file?
                        if(rel_stddev_total < max_error) {
                            cout << "The system equilibrated, the equilibration lastet to t = " << t << endl;
                            cout << "U_L = " << avg_U_L << " +- " << rel_stddev_total * avg_U_L << endl;
                            // we set the system to be equilibrated
                            sys.set_equilibration(t);
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
class old_cum_equilibration_observer: public obsver<system, State>{
    // Okay the observing pattern could be totally wild, i probably somehow have to initialize the observer
    // outside of the simulation class We definetely need its own constructor here
    typedef obsver<system, State> obsver;
    using obsver::ofile;
    using obsver::open_stream;
    using obsver::close_stream;
    double write_interval = 1;
    double timepoint = 0;
    vector<double> U_L{};
    vector<double> times{};
    int min_cum_nr = 500;
    int avg_nr = 5;     // after avg_nr of U_L calculations we check if we should adjust the stepsize
    double min_density = 1.0 / 1000.0;      // the minimum density of U_L calculations will be once in 1000 steps
    double max_density = 1.0 / 10.0;        // the maximum density of U_L calculations will be once in 10 steps
    double min_stepsize_div = 0.001;         // the stddeviation of avg_nr U_Ls should be at least min_stepsize_div of mean_UL
    double max_stepsize_div = 0.003;         // the stddeviation of avg_nr U_Ls should be at most max_stepsize_div of mean_UL
    double dt = 0.01;
    double equil_cutoff = 0.1;              // since the equilibration might influce the mean of U_L a lot we cut a certain portion of U_L values
    double max_error= 0.001;
    int cum_nr = 0;                         // current number in the averageing process
    bool equilibrated = false;                      // for the usecase of the quench with dynamic equilibration
public:

    void init(fs::path folderpath, map<Parameter, double>& paras, const system &sys) override {
        int run_nr = (int)paras[Parameter::run_nr];
        max_error = paras[Parameter::equil_error];
        if (paras[Parameter::min_cum_nr]){
            min_cum_nr = (int)paras[Parameter::min_cum_nr];
        }
        timepoint = 0.0;
        equilibrated = false;
        // we also need to reset U_L and times, dont we?
        U_L = vector<double>{};
        times = vector<double>{};
        close_stream();
        open_stream(folderpath / (obsver::construct_filename(run_nr) + ".cum"));
        cout << this->get_name() << " init called" << endl;
        ofile << "t,U_L" << endl;

        // I think the starting write interval should be every one hundred steps
        dt = paras[Parameter::dt];
        write_interval = 100 * dt;
    }

    string get_name() override {
        return "cum equilibration observer";
    }

    void operator()(system &sys, const State &x , double t ) override {
        if(t > timepoint) {
            // advancing the cum nr
            cum_nr++;
            // we calculate the cumulant and write it down
            double cum = sys.calc_binder(x);
            ofile << t << "," << cum << endl;
            // add the cumulant and the times to the vectors to keep track
            U_L.push_back(cum);
            times.push_back(t);
            // now we want to see if we have to adjust the write interval
            // we have to make sure that we got 5 fresh cum values
            if((cum_nr >= avg_nr) & !equilibrated) {
                cum_nr = 0;     // reset the cum nr
                int nr_cum_values = U_L.size();      // check how many cum values we already have
                // now use the last avg_nr of cum values to calculate a mean U_L
                // use transform reduce?
                // does it work like this? U_L.end() - avg_nr is the n-th last value in the vector?
                double mean_U_L = accumulate(U_L.end() - avg_nr, U_L.end(), 0.0) / (double)avg_nr;
                // calculate the stddev
                std::vector<double> diff(avg_nr);
                std::transform(U_L.end() - avg_nr, U_L.end(), diff.begin(), [mean_U_L](double x) { return x - mean_U_L; });
                double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);

                // we will use the relative stddeve
                double rel_stddev = sqrt(sq_sum / (double)avg_nr) / mean_U_L;
                // if the relative stddev is smaller than min_stepsize_div, we will increase the write interval
                if(rel_stddev < min_stepsize_div) {
                    cout << "increasing stepsize since rel_stddev = " << rel_stddev;
                    write_interval *= 2;
                    cout << "new stepsize = " << write_interval << endl;
                } else if (rel_stddev > max_stepsize_div) {
                    cout << "decreasing stepsize since rel_stddev = " << rel_stddev;
                    // if it is larger, we will reduce it
                    write_interval /= 2;
                    cout << "new stepsize = " << write_interval << endl;
                }
                // the write interval should be in bounds of min and max density
                write_interval = min((1.0 / min_density) * dt, write_interval);
                write_interval = max((1.0 / max_density) * dt, write_interval);

                // okay so we adjusted the write interval
                // now we still want to check if we are equilibrated
                // the simulation runs at least so long that there are min_cum_nr valuse
                if(nr_cum_values > min_cum_nr){
                    // the equilibration phase might influence the mean a lot, should we cut of the first x% of the values?
                    int min_ind = (int)(equil_cutoff * nr_cum_values);
                    // again calculate mean and stddev. We want the standarddeviation of the mean value this time?
                    double avg_U_L = accumulate(U_L.begin() + min_ind, U_L.end(), 0.0) / (double)(nr_cum_values - min_ind);
                    std::vector<double> diff_total(nr_cum_values - min_ind);
                    std::transform(U_L.begin() + min_ind, U_L.end(), diff_total.begin(), [avg_U_L](double x) { return x - avg_U_L; });
                    double sq_sum_total = std::inner_product(diff_total.begin(), diff_total.end(), diff_total.begin(), 0.0);
                    // this is the variance of the distribution (in wrong because the values are correlated), so we
                    // need the autocorrelation time to adjust for this. The function works with arrays
                    double* U_L_arr = &U_L[min_ind];
                    double autocorr_time = get_autocorrtime(U_L_arr, nr_cum_values - min_ind, write_interval);  // actually ds is just the write interval? which should be 1 or something like this
                    // (nr_cum_values-min_ind) * the write interval is the total time of the part of the simulation that we are considering
                    // sq sum total should be the thing that we called the standard deviation of the distribution
                    double U_L_variance = 2 * autocorr_time / ((nr_cum_values-min_ind) * write_interval) * sq_sum_total;
                    double rel_stddev_total = sqrt(U_L_variance) / avg_U_L;

                    // Okay we found out that this is not the real error, so we need a function that calculates the error
                    // or first a function that calculates the integrated autocorrleation time
                    // We indded have two auto correlation times for the two directions. Could they be individually sized?

                    cout << "rel_stddev_total = " << rel_stddev_total << endl;
                    // the question is now if we want to extract avg_U_L and its error if we are equilibrated
                    // we definitely should? But how and where? somehow into the parameter file?
                    if(rel_stddev_total < max_error) {
                        cout << "The system equilibrated, the equilibration lastet to t = " << t << endl;
                        cout << "U_L = " << avg_U_L << " +- " << rel_stddev_total * avg_U_L << endl;
                        // we set the system to be equilibrated
                        sys.set_equilibration(t);
                        // once we did this, we dont want to do that a second time?
                        equilibrated = true;
                        // actually we do not normally need any cumulant values in the quench case so we set
                        // the write interval to be just the maximum interval
                        write_interval = (1.0 / min_density) * dt;
                    }
                }
            }
            timepoint += write_interval;
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
    using obsver::close_stream;
    double write_interval = 1;
    double timepoint = 0;
    vector<double> xix{};
    vector<double> xiy{};
    vector<double> times{};
    int min_corr_nr = 5000;
    int val_nr_gpu = 25000;                 // nr of cum values at which the calculation switches to gpu
    double dt = 0.01;
    double dt_half = dt / 2.0;
    double equil_cutoff = 0.1;              // since the equilibration might influce the mean of xi
    double density = 1.0 / 10.0;            // standard density writes once every 100 steps
    double error_factor = 1000.0;
    // TODO we have to judge whether this is large or not. The thing is the error is good for low temperature states to
    //  judge whether we are equilibrated but bad for high temperature states since we have large deviations
    // for high temperature states we would usually need a larger number of systems to judge the equilibration
    // TODO could there be a way of combining the temperature with the error? A higher temperature would allow
    // a larger error on xi and still judge it to be equilibrated. But this is again kind of handwavy, how large temperatures
    // would result in how large leeway?
    // also the error has the weird property that large systems with small nubmers of subsystems equilibrate later.
    // For the Tc-Binder cumulant calculation this error is very suitable since we are looking for a precise U_L value anyway
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
    void init(fs::path folderpath, map<Parameter, double>& paras, const system &sys) override {
        int run_nr = (int)paras[Parameter::run_nr];
        max_error = paras[Parameter::equil_error];
        equil_cutoff = paras[Parameter::equil_cutoff];
        // the subsystems sizes for the xi cutting
        Lx = (size_t)paras[Parameter::subsystem_Lx];
        Ly = (size_t)paras[Parameter::subsystem_Ly];
        timepoint = 0.0;
        equilibrated = false;
        // we also need to reset U_L and times, dont we?
        xix = vector<double>{};
        xiy = vector<double>{};
        times = vector<double>{};
        close_stream();
        open_stream(folderpath / (obsver::construct_filename(run_nr) + ".corr"));
        cout << this->get_name() << " init called" << endl;
        ofile << "t,xix,xiy" << endl;

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
            sys.calc_xi(x, xix_val, xiy_val);
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
                    int error_every_n_steps = (int)(error_factor * density);
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

                        if(nr_xi_values - min_ind > 0) {
                            autocorr_time_x = get_autocorrtime_gpu(xix_arr, nr_xi_values - min_ind, write_interval);
                            autocorr_time_y = get_autocorrtime_gpu(xiy_arr, nr_xi_values - min_ind, write_interval);
                        } else {
                            autocorr_time_x = get_autocorrtime(xix_arr, nr_xi_values - min_ind, write_interval);
                            autocorr_time_y = get_autocorrtime(xiy_arr, nr_xi_values - min_ind, write_interval);  // actually ds is just the write interval? which should be 1 or something like this
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
class ft_observer : public obsver<system, State>{
    // Okay the observing pattern could be totally wild, i probably somehow have to initialize the observer
    // outside of the simulation class We definetely need its own constructor here
    typedef obsver<system, State> obsver;
    using obsver::ofile;
    using obsver::open_stream;
    using obsver::close_stream;
    int nr_values;
    double write_interval = 1;
    double timepoint = 0;
    size_t Lx;
    size_t Ly;
public:
    ft_observer(int nr_values) : nr_values(nr_values) {
    }
    void init(fs::path folderpath, map<Parameter, double>& paras, const system &sys) override {
        int run_nr = (int)paras[Parameter::run_nr];
        timepoint = 0.0;

        close_stream();
        open_stream(folderpath / (obsver::construct_filename(run_nr) + ".ft"));
        cout << "ft observer init called" << endl;
        ofile << "t;ft_k;ft_l" << endl;

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
class quench_ft_observer : public obsver<system, State>{
    // Okay the observing pattern could be totally wild, i probably somehow have to initialize the observer
    // outside of the simulation class We definetely need its own constructor here
    typedef obsver<system, State> obsver;
    using obsver::ofile;
    using obsver::open_stream;
    using obsver::close_stream;
    int nr_values;
    double write_interval = 1;
    double timepoint = 0;
    double end_quench_t;
    double quench_t;
    double s_eq_t;
    size_t Lx;
    size_t Ly;
public:
    quench_ft_observer(int nr_values) : nr_values(nr_values) {
    }
    void init(fs::path folderpath, map<Parameter, double>& paras, const system &sys) override {
        int run_nr = (int)paras[Parameter::run_nr];
        timepoint = 0.0;

        close_stream();
        open_stream(folderpath / (obsver::construct_filename(run_nr) + ".ft"));
        cout << "quench ft observer init called" << endl;
        ofile << "t;ft_k;ft_l" << endl;

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
        write_interval = 100 * dt;      // We probably should make this density variable
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
            // write
            if(startpoint){
                this->write_state(sys, x, t);
                startpoint = false;
            } else if(t > eq_t - dt_half) {
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
    void init(fs::path folderpath, map<Parameter, double>& paras, const system &sys) override {
        int run_nr = (int)paras[Parameter::run_nr];
        timepoint = 0.0;

        close_stream();
        open_stream(folderpath / (obsver::construct_filename(run_nr) + ".ft"));
        cout << "density quench ft observer init called" << endl;
        ofile << "t;ft_k;ft_l" << endl;

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

/*class quench_observer : public observer {
    // Okay the observing pattern could be totally wild, i probably somehow have to initialize the observer
    // outside of the simulation class We definetely need its own constructor here
    double write_interval = 1;
    int nr_values;
    double timepoint = 0;
public:
    quench_observer(int nr_values): nr_values(nr_values) {
        // should this be the observer that only writes a few values for the equilibrium process?
        // I am not sure at the moment, I think i am just going to write the simplest one for now, so
        // just equidistant saving values.
        // We try to initialize the observer outside the class since it probably needs
        // multiple unique-to-observer parameters
        // problem is it still might depend on system parameters. and i don't want another if for every step
        // even though i am not really sure if it impacts the performance AT ALL or if the compiler knows better
    }

    template <class System, class State>
    void operator()(const System &sys, const State &x , double t )  {
        if(t > timepoint) {
            // write
            // Doing some refactoring, not putting t in every row and Temperature directly after T
            double T = sys.get_cur_T();
            int lat_dim = sys.get_lat_dim();
            ofile << t << "," << T << ",";
            for(int i = 0; i < lat_dim * lat_dim; i++) {
                ofile << x[i] << ",";
                // cout << x[i] << endl;

            }
            // Zeilenumbruch
            ofile << "\n";
            // advance timeoint
            timepoint += write_interval;
        }
    }
    template <class System>
    void init(fs::path folderpath, map<string, double>& paras, const System &sys)  {
        // I think this will have less performance impact than an if statement catching the first observer operation
        // Make sure to use this observer only with systems that have a get_end_T method
        cout << "Quench obs init is called" << endl;
        double end_T = sys.get_end_t();
        write_interval = end_T / (double)nr_values;
        timepoint = 0.0;

        // open the file to write the info to, in this case it will be just run_nr.csv
        // I think we will add the run number to the paras of the run
        int run_nr = (int)paras["run_nr"];
        open_stream(folderpath / (to_string(run_nr) + ".csv"));
    }

    void init(fs::path folderpath, map<string, double>& paras)  {
        cout << "call this mf" << endl;
    }

    string get_name() {
        return "quench observer";
    }
};*/
// Observer for lattic on bath specific for gpu_bath system?
class bath_observer : public observer {
    // this one has a filestream
    ofstream& file;
    // and a write_every
    const size_t write_every;
    // and a count? Is that all efficient?
    size_t count;

public:
    bath_observer(ofstream& out, size_t write_every = 1000) : file(out), write_every(write_every) {
        count = 1;
    }

/*    bath_observer(size_t write_every = 1000) : write_every(write_every) {
    }*/

    template<class State, class System>
    void operator()(System &sys, const State &x , double t ) {
        // so i think if statements are not to costly on c++
        if (count % write_every == 0 || count == 1) {
            // now we write the state to the file
            // i call a function here that would not be available for all Systems, not to clean and a bit fishy
            // but i guess it works for now
            size_t lat_dim = sys.get_lattice_dim();
            double T = sys.get_cur_T();

            // writing, but only q for the moment
            file << "t, " << t << ",";
            for(int i = 0; i < lat_dim * lat_dim; i++) {
                file << x[i] << ",";
                // cout << x[i] << endl;

            }
            // for the last we write T? I mean thats how i did it, right?
            // Temperatur saven
            file << T;

            // Zeilenumbruch
            file << "\n";
        }
        count++;
    }
    template<class State, class System>
    void write(System &sys, const State &x , double t ) {
        size_t lat_dim = sys.get_lattice_dim();
        double T = sys.get_cur_T();

        file << "t, " << t << ",";
        for(int i = 0; i < lat_dim * lat_dim; i++) {
            file << x[i] << ",";
            // cout << x[i] << endl;

        }
        // for the last we write T? I mean thats how i did it, right?
        // Temperatur saven
        file << T;

        // Zeilenumbruch
        file << "\n";
    }

    template<class State, class System>
    void writev2(System &sys, const State &x , double t ) {
        // Doing some refactoring, not putting t in every row and Temperature directly after T
        size_t lat_dim = sys.get_lattice_dim();
        double T = sys.get_cur_T();

        file << t << "," << T << ",";
        for(int i = 0; i < lat_dim * lat_dim; i++) {
            file << x[i] << ",";
            // cout << x[i] << endl;

        }
        // Zeilenumbruch
        file << "\n";
    }
};


/*class quench_observer : public bath_observer {
    ofstream &xi_out;
public:
    quench_observer(ofstream& sys_out, ofstream& xi_out) : bath_observer(sys_out), xi_out(xi_out){
        // write the header
        xi_out << "t, xi" << endl;
    }

    template<class State, class System>
    void write_xi(System &sys, const State &x, double t) {
        // we need to calculate the correlation function and the correlation length
        // correlation function is easy but length is hard, need to do fit
        // the thing is, x is a gpu vector and to copy it to host to use my old corr_func might be slow
        // rewriting in thrust will take some time tho...
        // probably for now best to use slow variant since we probably switch to fftw sometime anyway?
        // so we use thrust copy to do this
        size_t lat_dim = sys.get_lattice_dim();
        size_t n = lat_dim * lat_dim;
        size_t nr_dists = lat_dim / 2 + 1;
        vector<double> vals(n);         // init standard vector with size of x_values which is x.size()/2

        // copy from x to f
        thrust::copy(x.begin(), x.begin() + n, vals.begin());      // okay, should be copied, we have one d f
        vector<vector<double>> f = oneD_to_twoD(vals);
        // ready to calc corr func, need vectors to save it
        Eigen::VectorXd C_x = Eigen::VectorXd::Zero((long)nr_dists);
        Eigen::VectorXd C_y = Eigen::VectorXd::Zero((long)nr_dists);

        calc_corr(f, C_x, C_y);
        // average them for now
        Eigen::VectorXd C = average_two(C_x, C_y);
        // vector with distances?
        Eigen::VectorXd dists = Eigen::VectorXd::LinSpaced((long)nr_dists, 0, (int)nr_dists - 1);
        // construct matrix out of the two
        Eigen::MatrixXd C_x_vals(nr_dists, 2);
        C_x_vals = construct_matrix(C, dists);
        // now we need to do the fit
        // init starting values
        Eigen::VectorXd params(2);
        params << 1.0, 5.0;
        NumericalExpDecayDiffFunctor functor(C_x_vals);
        Eigen::LevenbergMarquardt<NumericalExpDecayDiffFunctor> lm(functor);
        Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(params);
        std::cout << "status: " << status << std::endl;

        //std::cout << "info: " << lm.info() << std::endl;

        std::cout << "params that minimizes the function: " << std::endl << params << std::endl;

        // now printing to file
        // TODO i am now done and now i realize that it isnt even that smart to fit for every system...
        // TODO we should actually average the Correlation function and fit afterwards
        // but good that you now know how to fit in c++...
        // params(1) is the correlation length that we want, so we write this to the file
        xi_out << t << "," << params(1) << endl;
    }
};*/


// Observer for lattic on bath specific for gpu_bath system?
class chain_observer : public observer {
    // this one has a filestream
    ofstream& file;
    // and a write_every
    const size_t write_every;
    // and a count? Is that all efficient?
    size_t count;

public:
    chain_observer(ofstream& out, size_t write_every = 1000) : file(out), write_every(write_every) {
        count = 1;
    }

    template<class State, class System>
    void operator()(const System sys, const State &x , double t ) {
        // so i think if statements are not to costly on c++
        if (count % write_every == 0 || count == 1) {
            // now we write the state to the file
            // i call a function here that would not be available for all Systems, not to clean and a bit fishy
            // but i guess it works for now
            size_t lat_dim = sys.get_lattice_dim();
            double T = sys.get_cur_T();

            // writing, but only q for the moment
            file << "t, " << t << ",";
            for(int i = 0; i < lat_dim; i++) {
                file << x[i] << ",";
                // cout << x[i] << endl;

            }
            // for the last we write T? I mean thats how i did it, right?
            // Temperatur saven
            file << T;

            // Zeilenumbruch
            file << "\n";
        }
        count++;
    }

    template<class State, class System>
    void write(const System sys, const State &x , double t ) {
        size_t lat_dim = sys.get_lattice_dim();
        double T = sys.get_cur_T();

        file << "t, " << t << ",";
        for(int i = 0; i < lat_dim; i++) {
            file << x[i] << ",";
            // cout << x[i] << endl;

        }
        // for the last we write T? I mean thats how i did it, right?
        // Temperatur saven
        file << T;

        // Zeilenumbruch
        file << "\n";
    }

};


#endif //CUDAPROJECT_OBSERVERS_CUH
