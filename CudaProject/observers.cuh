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
};

template <class system, class State>
class standard_observer : public obsver<system, State> {
    double timepoint = 0;
public:
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
        open_stream(folderpath / construct_filename(run_nr) + ".txt"));
        write_parameters(ofile, paras);
        // dont forget to close!;
        close_stream();

        open_stream(folderpath / construct_filename(run_nr) + ".csv"));
    }
    void operator()(system &sys, const State &x , double t ) override {
        if(t > timepoint) {
            // write
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
            // advance timeoint
            timepoint += write_interval;
        }
    }

    string get_name() override {
        return "standard observer";
    }

    string construct_filename(int run_nr) {
        return to_string(run_nr) + "-" + get_current_time() + "-" + boost::asio::ip::host_name();
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
};

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