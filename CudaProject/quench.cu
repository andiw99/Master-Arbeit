//
// Created by weitze73 on 08.06.23.
//
#include "main.cuh"
#include "systems.cuh"
#include "parameters.cuh"

namespace fs = std::filesystem;

int rep;

template <  size_t lattice_dim,
            class state_type, class alg, class oper,
            template<class statetype, class algebra, class operations> class stepper,
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
         paras["starting_T"] << " to " << paras["end_T"] << " for a time of " << paras["t"] << endl;
    // we should define a equilibrate_t for start and end
    // we dont need a general t now since the time of the simulation is determined by the
    // two equilibrate times and the quench time (which is defined by tau and start and endpoint)
    // we now need to worry about where to save the file
    create_dir(save_dir);
    // and the names of the files
    fs::path name = save_dir / to_string(rep);
    // initialize the output streams
    ofstream file, parafile, corrfile;
    file.open(name / ".csv");
    parafile.open(name / ".txt");
    corrfile.open(name / "-corr.txt");
    // init the observer
    obsver obs(file, corrfile);
    // init the stepper
    stepper<state_type, alg, oper> Stepper(n * N, paras["K"], paras["tol"]);
    // initial state, we this time want it to have approximately the energy of the system in thermal equilibrium
    double x0 = 0;  // double well potential makes it hard to use it to init a state of certain energy
    double p0 = sqrt(M_PI * paras["starting_T"] / 4);       // energy impuls relation makes it easy to init
    // init empty vector                                       // a state that has a certain energy expectation value
    state_type x(N * n);                                       // since we now the expectation value of |p|
    // fill it
    fill_init_values<state_type, lattice_dim * lattice_dim>(x, (float) x0, (float) p0);
    // we should check for the energy?
    // this is a transform reduce operation, but is it worth implementing it?
    // maybe in the long run...
    cout << "Initial kinetic energy: E_kin = " << sys.calc_kinetic_energy(x) << " vs theoretical thermic energy "
    << (double)n * paras["starting_T"] << endl;
    // now we can do the stepping?
    double t = 0, xi_timepoint = 0, sys_timepoint = 0;
    double end_t = sys.get_end_t();
    double write_interval_xi = end_t / paras["nr_xis"];
    double write_interval_sys = end_t / paras["nr_save_values"];
    // we need to know when to write, we probably want more xis than whole system data?
    while(t < end_t) {
        // write first
        if(t > xi_timepoint) {
            obs.write_xi(sys, x, t);
            xi_timepoint += write_interval_xi;
            // since we usually want more xi values, we only check to write the system every
            // time we also check to write xi
            if(t > sys_timepoint) {
                obs.write(sys, x, t);
                sys_timepoint += write_interval_xi;
            }
        }

        // now we do the step
        Stepper.do_step(sys, x, dt_max, t);
    }
}

template <  size_t lattice_dim,
        class state_type, class alg, class oper,
        template<class statetype, class algebra, class operations> class stepper,
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
            >(parameters, runs - 1, seed, dir_path);
}

void scan() {
    const size_t lattice_dim = 128;
    // typedefs here?
    typedef thrust::device_vector<double> state_type;
    typedef thrust_algebra algebra;
    typedef thrust_operations operations;
    typedef euler_combined<state_type, algebra, operations> stepper;
    typedef quench_observer observer;


    fs::path root = quench_root;

    map<string, double> paras = quench_paras;

    // so i think we want to space the tau logarithmically, meaning we either
    // give in the minimum and maximum value and have to recalc here
    // or we just enter the factors and have to think while making parameters
    const vector<double> taus = linspace(paras["min_tau_factor"],
                                      paras["max_tau_facotr"],
                                      (int)paras["nr_taus"] + 1);
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
    return 0;
}