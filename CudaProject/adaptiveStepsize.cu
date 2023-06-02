//
// Created by andi on 02.06.23.
//
#include "main.cuh"
#include "systems.cuh"
#include "parameters.cuh"

template <size_t lattice_dim>
int adaptive_routine(map<string, double> parameters, long seed = 0, string system="default", string save_dir = "", int count=0) {
    // We try out the code for the brownian motion i would say
    // But we cannot use our old class system I think because there the whole system is already on a lattice
//
// Created by andi on 21.04.23.
    // But we can quickly write another system i guess
    double    dt_max = parameters["dt_max"];
    const double    T = parameters["T"];
    const double    J = parameters["J"];
    const double    alpha = parameters["alpha"];
    const double    beta = parameters["beta"];
    const double    tau = parameters["tau"];
    const double    eta = parameters["eta"];
    const size_t    N = (size_t) parameters["N"];
    const int       nr_save_values = (int)parameters["nr_save_values"];
    // for the adaptive routine we have two new parameters, K_start and tol
    const int K = (int) parameters["K"];
    double tol = parameters["tol"];
    double end_t = parameters["end_t"];


    const double x0 = parameters["x0"];
    const double p0 = parameters["p0"];
    const           size_t n = lattice_dim * lattice_dim;

    double write_interval = parameters["end_t"] / nr_save_values;
    cout << "Starting Simulation for a " << lattice_dim << " by " << lattice_dim << " lattice until" << end_t << endl;

    // last time i didnt have to specify the dimensionality i think (in terms of (x, p) )
    const double D = T / eta;
    double mu = 0;
    double msd = 0;
    // file stuff
    string dir_name;
    if(save_dir.empty()) {
        string storage_root = "../../../Generated content/Default/";
        dir_name = create_tree_name(eta, T, dt_max, n, alpha, beta, J, tau,
                                    storage_root);
    } else {
        dir_name = save_dir;
        // create directory?
        create_dir(dir_name);
    }
    string name = dir_name + "/" + to_string(count);
    ofstream file;
    ofstream parafile;
    file.open(name + ".csv");
    parafile.open(name + ".txt");

    bath_observer Obs(file, write_interval);

    typedef thrust::device_vector<double> gpu_state_type;

    euler_simple_adaptive<gpu_state_type, thrust_algebra, thrust_operations > gpu_stepper(N * n, K, tol);

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
    double t = 0;
    double write_timepoint = 0;


    // i THINK we will do the ugly way here until we understand polymorphism
    cout << "creating System " << system << endl;
    if (system == "constant") {
        constant_bath<lattice_dim> gpu_system(T, eta, alpha, beta, J, seed);
        while(t < end_t) {
            // we need small stepsizes at the beginning to guarantee stability
            // but after some time, we can increase the stepsize
            // there should be a t that is equal to t_relax?

            gpu_stepper.do_step(gpu_system, x, dt_max, t);
            if (t >= write_timepoint) {
                Obs.write(gpu_system, x, t);
                write_timepoint += write_interval;
                cout << "current k = " << gpu_stepper.get_k() << endl;
            }
        }
    }
    else if(system == "quadratic_chain") {
        cout << "creating quadratic chain obj" << endl;
        quadratic_chain<lattice_dim> gpu_system(T, eta, J, seed);
        while(t < end_t) {
            gpu_stepper.do_step(gpu_system, x, dt_max, t);
            Obs(gpu_system, x, t);
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
        while(t < end_t) {
            gpu_stepper.do_step(gpu_system, x, dt_max, t);
            Obs(gpu_system, x, t);
        }
    }

    write_parameters(parafile, eta, T, dt_max, n, alpha, beta, J, tau);

    // print stuff
    {
        file.close();


        mu = 0;
        msd = 0;
        for(int i = 0; i < n; i++) {
            mu += x[i];
            msd += x[i] * x[i];
        }
        mu /= n;
        msd /= n;

        cout << "mu = " << mu << endl;
        return (int)(end_t/dt_max);
    }
}

template <size_t lattice_dim>
void repeat(map<string, double> parameters, int runs, long seed = 0, string system="default", string dir_path="", int count=0) {
    // seed is the seed for the random numbers so that we can have different random numbers per run
    if(runs == 0) {
        return;
    }

    cout << runs << " runs left" << endl;

    int steps = adaptive_routine<lattice_dim>(parameters, seed, system, dir_path);
    // how to get the number of steps that were done? let single calc routine return it?
    // or put it also into repeat
    repeat<lattice_dim>(parameters, runs - 1, seed + steps, system, dir_path);

}


void simple_temps_scan() {
    // we always need to specify the lattice dim
    const size_t lattice_dim = 100;

    string root = adaptive_tempscan_root;

    map<string, double> paras = temp_scan_standard;
    // we do not use the fast forward here
    paras["dt"] = adaptive_temp_scan_standard["dt_start"];

    const vector<double> T = linspace(temp_scan_standard["min_temp"],
                                      temp_scan_standard["max_temp"], (int)temp_scan_standard["nr_temps"] + 1);
    // printing
    {
        print_vector(T);

        cout << "Starting Simulation for a " << lattice_dim << " by " << lattice_dim << " lattice unitl t = " << paras["end_t"] << endl;
        cout << "Stepsize is dt = " << paras["dt"] << endl;
    }

    // now cycling
    for(double temp : T) {
        paras["T"] = temp;
        // for every T we need to initalize a new system, but first i think we maybe should check our system?
        // for every T we need another directory
        string dirpath = root + "/" + to_string(temp);

        cout << "Running repeat with following parameters:" << endl;
        printMap(paras);
        repeat<lattice_dim>(paras, (int)temp_scan_standard["repeat_nr"], 0, "constant", dirpath);
    }
}

int main() {
    return 0;
}