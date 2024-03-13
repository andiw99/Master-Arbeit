//
// Created by weitze73 on 08.06.23.
//
#include "../main.cuh"
#include "../systems-cuda.cuh"
#include "parameters.cuh"

namespace fs = std::filesystem;

int rep;

template <  size_t lattice_dim,
            class state_type, class alg, class oper,
            template<class, class, class, class, class> class stepper,
            template<size_t lat_dim> class system,
            class obsver>
void single_quench(map<string, double> paras, fs::path &save_dir, size_t seed = 0) {
    // I want to run this for:
    // different tau (most important)
    // different starting Ts (less important)
    // in general different values. Can there be an architecture that runs them dynamically?
    // lets just assume everything we need to know is written in parameters
    // we need to initialize the starting values
    // we somehow need to make sure that we reach end_T within t, it depends on tau
    // only some important values:
    size_t n = lattice_dim * lattice_dim;
    size_t N = (size_t)paras["N"];
    double dt_max = paras["dt"];
    // we initialize the system
    system<lattice_dim> sys(
                paras["starting_T"],
                paras["end_T"],
                paras["eta"],
                paras["alpha"],
                paras["beta"],
                paras["J"],
                paras["tau"],
                seed,
                paras["t_eq"]
            );
    // it does not make sense to quench for a certain time
    cout << "Starting Simulation for a " << lattice_dim << " by " << lattice_dim << " lattice from " <<
         paras["starting_T"] << " to " << paras["end_T"] << " for a time of " << sys.get_end_t() << endl;
    cout << "Quenching takes a time of " << sys.get_quench_time() << endl;
    // we should define a equilibrate_t for start and end
    // we dont need a general t now since the time of the simulation is determined by the
    // two equilibrate times and the quench time (which is defined by tau and start and endpoint)
    // we now need to worry about where to save the file
    create_dir(save_dir);
    // and the names of the files
    fs::path name = save_dir / to_string(rep);
    // initialize the output stream
    cout << "Trying to save to " << name << endl;
    ofstream file, parafile;
    file.open(name.replace_extension(".csv"));
    parafile.open(name.replace_extension(".txt"));
    // init the observer
    obsver obs(file);
    // init the stepper
    stepper<state_type, alg, oper, double, double> Stepper(n * N, paras["K"], paras["tol"]);
    // initial state, we this time want it to have approximately the energy of the system in thermal equilibrium
    // or at least above?
    double x0 = 0;  // double well potential makes it hard to use it to init a state of certain energy
    double p0 = sqrt(M_PI * paras["starting_T"]);       // energy impuls relation makes it easy to init
    // init empty vector                                       // a state that has a certain energy expectation value
    state_type x(N * n);                                       // since we now the expectation value of |p|
    // fill it
    chrono::milliseconds ms = chrono::duration_cast<chrono::milliseconds >(
            chrono::system_clock::now().time_since_epoch()
    );
    fill_init_values<state_type, lattice_dim * lattice_dim>(x, (float) x0, (float) p0, ms.count() % 10000);
    // we should check for the energy?
    // this is a transform reduce operation, but is it worth implementing it?
    // maybe in the long run...
    cout << "Initial kinetic energy: E_kin = " << sys.calc_kinetic_energy(x) << " vs theoretical thermic energy "
    << (double)n * paras["starting_T"] << endl;
    // now we can do the stepping?
    double end_t = sys.get_end_t();
    double s_eq_t = paras["t_eq"];
    double e_eq_t = paras["t_eq"];
    double quench_time = sys.get_quench_time();
    double end_quench_t = sys.get_end_quench_time();
    double write_interval_xi = end_t / paras["nr_xis"];
    double write_interval_sys = quench_time / (paras["nr_save_values"] - 5);
    double t = 0, xi_timepoint = 0, sys_timepoint = s_eq_t;
    // we need to know when to write, we probably want more xis than whole system data?
    // We might want to make the observing a bit more intelligent. We want one saved value
    // for the initial state, then two during start equilibrating, including the timepoint s_eq_t
    // then 11 of the quench, and then again two for the end equilibration

    obs.writev2(sys, x, t);
    Stepper.step_until(s_eq_t/2, sys, x, dt_max, t);
    obs.writev2(sys, x, t);
    Stepper.step_until(s_eq_t, sys, x, dt_max, t);
    while(t < end_quench_t) {
        // we step until sys_timepoint, write and increment for the next sys_timepoint
        Stepper.step_until(sys_timepoint, sys, x, dt_max, t);
        obs.writev2(sys, x, t);
        sys_timepoint += write_interval_sys;
    }
    cout << end_quench_t + e_eq_t/2 << endl;
    Stepper.step_until(end_quench_t + e_eq_t/2, sys, x, dt_max, t);
    obs.writev2(sys, x, t);
    Stepper.step_until(end_t, sys, x, dt_max, t);
    obs.writev2(sys, x, t);

    // are we done now? no write parameters
    write_parameters(parafile, paras);
}


template <  size_t lattice_dim,
        class state_type, class alg, class oper,
        template<class, class, class, class, class> class stepper,
        template<size_t lat_dim> class system,
        class obsver>
void repeat(map<string, double> parameters, int runs, size_t seed = 0, fs::path dir_path="") {
    // seed is the seed for the random numbers so that we can have different random numbers per run
    if(runs == 0) {
        return;
    }

    cout << runs << " runs left" << endl;

    single_quench<
                                lattice_dim,
                                state_type, alg, oper,
                                stepper,
                                system,
                                obsver
                                >(parameters, dir_path, seed);
    // increase the repetition number
    rep++;
    // how to get the number of steps that were done? let single calc routine return it?
    // or put it also into repeat
    repeat<lattice_dim,
            state_type, alg, oper,
            stepper,
            system,
            obsver
            >(parameters, runs - 1, 100000000 * rep, dir_path);    // TODO better seed
}

void scan() {
    // typedefs here?
    typedef thrust::device_vector<double> state_type;
    typedef thrust_algebra algebra;
    typedef thrust_operations operations;
    typedef euler_combined<state_type, algebra, operations> stepper;
    // I don't think we need another observer, we just do what we did the whole time, but when we calculate the
    // correlation function, we do it for every time we saved
    // Advantages: We have the data which generated the correlation lengths
    // Disadvatages: We have to save them and possibly large data structure.
    typedef bath_observer observer;


    fs::path root = quench_root;

    map<string, double> paras = quench_paras;

    // so i think we want to space the tau logarithmically, meaning we either
    // give in the minimum and maximum value and have to recalc here
    // or we just enter the factors and have to think while making parameters
    // how long can a Quench take? 10000s should be an upper bound, when quenching for 100K this means
    // 10000 / 100 = tau_max = 100
    // New considerations
    // Since we seemingly have to average over a lot of runs, 1000s should be an upper bound, when quenching for
    // 1K this meens tau_max = 1000
    // logspace: tau = 2 ^ tau_factor, so max_tau_factor =  7 (10)
    // maybe if we later know our system better we can go to quench for less temp, like 10K, 20K
    // we just have to make sure we are outside the freezeout time?
    const vector<double> taus = logspace(paras["min_tau_factor"],
                                      paras["max_tau_factor"],
                                      (int)paras["nr_taus"] + 1);
    cout << "Defined taus = " << endl;
    print_vector(taus);

    for(double tau : taus){
        rep = 0;
        paras["tau"] = tau;
        string dirpath = root / to_string(tau);
        cout << "parameters:" << endl;
        printMap(paras);
        // call repeat with right templates
        repeat< lattice_dim,
                state_type, algebra, operations,
                euler_combined,
                gpu_bath,
                observer
                >(paras, (int)paras["repeat"], 0, dirpath);
    }

}



int main() {
    scan();
    return 0;
}