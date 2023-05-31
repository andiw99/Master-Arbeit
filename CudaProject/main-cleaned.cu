//
// Created by andi on 14.05.23.
//
#include "main.cuh"
#include "systems.cuh"
#include "parameters.cuh"
#include <numeric>
#include <chrono>
#include <map>


template <size_t lattice_dim>
int single_calc_routine(map<string, double> parameters, long seed = 0, string system="default", string save_dir = "") {
    // We try out the code for the brownian motion i would say
    // But we cannot use our old class system I think because there the whole system is already on a lattice
//
// Created by andi on 21.04.23.
    // But we can quickly write another system i guess
    const int       steps = (int)parameters["steps"];
    double    dt = parameters["dt"];
    const double    T = parameters["T"];
    const double    J = parameters["J"];
    const double    alpha = parameters["alpha"];
    const double    beta = parameters["beta"];
    const double    tau = parameters["tau"];
    const double    eta = parameters["eta"];
    const size_t    N = (size_t) parameters["N"];
    const int       nr_save_values = (int)parameters["nr_save_values"];
    const double x0 = parameters["x0"];
    const double p0 = parameters["p0"];
    const           size_t n = lattice_dim * lattice_dim;

    size_t write_every = steps / nr_save_values;
    cout << "Starting Simulation for a " << lattice_dim << " by " << lattice_dim << " lattice for " << steps << " steps." << endl;

    // last time i didnt have to specify the dimensionality i think (in terms of (x, p) )
    const double D = T / eta;
    double theo_msd = 2 * D * dt * steps;
    double mu = 0;
    double msd = 0;
    // file stuff
    string dir_name;
    if(save_dir.empty()) {
        string storage_root = "../../../Generated content/Default/";
        dir_name = create_tree_name(eta, T, dt, n, alpha, beta, J, tau,
                                    storage_root);
    } else {
        dir_name = save_dir;
        // create directory?
        create_dir(dir_name);
    }
    string name = dir_name + "/" + to_string(seed/steps);
    ofstream file;
    ofstream parafile;
    file.open(name + ".csv");
    parafile.open(name + ".txt");

    bath_observer Obs(file, write_every);

    typedef thrust::device_vector<double> gpu_state_type;
    euler_mayurama_stepper<gpu_state_type, thrust_algebra, thrust_operations > gpu_stepper(N * n);

    // init and print initial state

    gpu_state_type x(N * n, x0);
    // set the impulses to be zero
    thrust::fill(x.begin() + n, x.begin() + N * n, p0);
    // okay we overwrite this here
    fill_init_values<n>(x, (float) x0, (float) p0);
    for (int i = 0; i < n; i++) {
        mu += x[i];
        msd += x[i] * x[i];
    }
/*    cout << "Initial values:" << endl;
    cout << mu / (n) << endl;
    cout << msd / (n) << endl;*/


    auto start = chrono::high_resolution_clock::now();
    double t = 0;


    // i THINK we will do the ugly way here until we understand polymorphism
    cout << "creating System " << system << endl;
    if (system == "constant") {
        constant_bath<lattice_dim> gpu_system(T, eta, alpha, beta, J, seed);
        for( size_t i=0 ; i<steps ; ++i ) {
            // we need small stepsizes at the beginning to guarantee stability
            // but after some time, we can increase the stepsize
            // there should be a t that is equal to t_relax?
            double dt_max = parameters["dt_max"];
            int t_relax = (int)parameters["t_relax"];
            if(t >= t_relax) {
                dt = dt_max;
            }
            gpu_stepper.do_step(gpu_system, x, dt, t);
            Obs(gpu_system, x, t);
            t += dt;
        }
    }
    else if(system == "quadratic_chain") {
        cout << "creating quadratic chain obj" << endl;
        quadratic_chain<lattice_dim> gpu_system(T, eta, J, seed);
        for( size_t i=0 ; i<steps ; ++i ) {
            gpu_stepper.do_step(gpu_system, x, dt, t);
            Obs(gpu_system, x, t);
            t += dt;
        }
        // calculate the energy and print it
        {
            double E = gpu_system.calc_energy(x);
            double Ekin = gpu_system.calc_kinetic_energy(x);
            double Epot = gpu_system.calc_potential_energy(x);
            double d2 = gpu_system.calc_total_squared_dist(x);
            cout << "Energy of the System: " << E << endl;
            cout << "Theoretical Energy: " << (double) lattice_dim * T << endl;
            cout << "kinetic Energy of the System: " << Ekin << endl;
            cout << "potential Energy of the System: " << Epot << endl;
            cout << "total squared dist: " << d2 << endl;
            cout << "theoretical total squared dist: " << (double) lattice_dim * T / J;
        }
    } else {
        gpu_bath<lattice_dim> gpu_system(T, eta, alpha, beta, J, tau, seed);
        for( size_t i=0 ; i<steps ; ++i ) {
            gpu_stepper.do_step(gpu_system, x, dt, t);
            Obs(gpu_system, x, t);
            t += dt;
        }
    }

    write_parameters(parafile, eta, T, dt, n, alpha, beta, J, tau);

    // print stuff
    {
        file.close();
        auto end = chrono::high_resolution_clock::now();
        auto duration = chrono::duration_cast<std::chrono::milliseconds>(end - start);
        cout << "execution took " << duration.count() << "ms, meaning " <<
             duration.count() * 1000/(steps) << "ms per 1000 steps." << endl;
        cout << "for a " << lattice_dim << " by " << lattice_dim << " lattice." << endl;
        // print this shit
        // TODO we could use this reduction stuff to compute the moments
        mu = 0;
        msd = 0;
        for(int i = 0; i < n; i++) {
            mu += x[i];
            msd += x[i] * x[i];
        }
        mu /= n;
        msd /= n;

        cout << "mu = " << mu << endl;
        cout << "msd = " << msd << "   theo value msd = " << theo_msd << endl;
        return steps;
    }
}

template <size_t lattice_dim>
int single_calc_routine(long seed = 0, string system="default", string save_dir="") {

    map<string, double> paras = single_calc_standard;

    return single_calc_routine<lattice_dim>(paras, seed, system, save_dir);
}

template <size_t lattice_dim>
void repeat(int runs, long seed = 0, string system="default", string save_dir="") {
    // seed is the seed for the random numbers so that we can have different random numbers per run
    if(runs == 0) {
        return;
    }

    cout << runs << " runs left" << endl;

    int steps = single_calc_routine<lattice_dim>(seed, system, save_dir);
    // how to get the number of steps that were done? let single calc routine return it?
    // or put it also into repeat
    repeat<lattice_dim>(runs - 1, seed + steps, system, save_dir);

}

template <size_t lattice_dim>
void repeat(map<string, double> parameters, int runs, long seed = 0, string system="default", string dir_path="") {
    // seed is the seed for the random numbers so that we can have different random numbers per run
    if(runs == 0) {
        return;
    }

    cout << runs << " runs left" << endl;

    int steps = single_calc_routine<lattice_dim>(parameters, seed, system, dir_path);
    // how to get the number of steps that were done? let single calc routine return it?
    // or put it also into repeat
    repeat<lattice_dim>(parameters, runs - 1, seed + steps, system, dir_path);

}

void scan_temps_routine(const int steps_val = 0, const int end_t_val = 0, const string& root_val = "", double mu = 0, double sigma=1) {

    string root = tempscan_root;
    root = (root_val.empty()) ? root : root_val;
    create_dir(root);

    const size_t lattice_dim = 100;


    double t_relax = temp_scan_standard["t_relax"];
    int end_t = (end_t_val == 0) ? (int)temp_scan_standard["end_t"] : end_t_val;

    map<string, double> paras = temp_scan_standard;

    int steps = (int)(t_relax / temp_scan_standard["dt_start"] + (end_t - t_relax) / temp_scan_standard["dt_start"]);

    paras["steps"] = steps;
    paras["dt"] = temp_scan_standard["dt_start"];

    const vector<double> T = linspace(temp_scan_standard["min_temp"],
                                      temp_scan_standard["max_temp"], (int)temp_scan_standard["nr_temps"] + 1);

    print_vector(T);

    cout << "Starting Simulation for a " << lattice_dim << " by " << lattice_dim << " lattice for " << steps << " steps." << endl;
    cout << "Stepsize is dt = " << paras["dt"] << endl;



    // keep track at which
    int count = 0;
    for(double temp : T) {
        paras["T"] = temp;
        // for every T we need to initalize a new system, but first i think we maybe should check our system?
        // for every T we need another directory
        string dirpath = root + "/" + to_string(temp);
        // we could repeat?
        // we need a new observer with a file name
        cout << "Running repeat with following parameters:" << endl;
        printMap(paras);
        repeat<lattice_dim>(paras, (int)temp_scan_standard["repeat_nr"], 0, "constant", dirpath);
        count++;
    }
}


void convergence_check_oscillatorchain() {
    // We can reuse the temp scan here, we just do a scan with also multiple lattice sizes and stepsizes
    // lattice sizes have to be clear at compile time, so we have to run it multiple times for the different sizes
    // but we can scan for temperature and stepsize
    // it was pretty dampened for t = 12, i would say t_end will be 15
    int t_end = 15;
    // steps von 0.02 - 0.0001
    // dt  = t_end / steps -> steps= t_end / dt
    string root = "../../../Generated content/Convergence Check MSD/";
    vector<int> steps = linspace((int)(t_end / 0.02), (int)(t_end / (0.001)), 3);
    for(int step : steps) {
        cout << step << endl;
        scan_temps_routine(step, t_end, root, 1.0, 0.0);
    }
}


void quadratic_chain_routine() {
    const size_t N = 10000;
    int runs = 10;
    string save_dir = "../../Generated content/quadratic chain/";
    // if we here just use the standard repeat function, we have to set the parameters in the standard single
    // calc function
    // TODO this is all not optimal, I need to rewrite this soon
    repeat<N>(runs, 0, "quadratic_chain", save_dir);
}


int main() {
    scan_temps_routine();
    // convergence_check_oscillatorchain();
    // quadratic_chain_routine();


    // single_calc_routine();
    // add to msd and mu
    return 0;
}
