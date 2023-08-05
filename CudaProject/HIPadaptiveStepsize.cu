//
// Created by andi on 02.06.23.
//
#include "main.cuh"
#include "systems.cuh"
#include "parameters.cuh"

template <template<class, class, class, class, class> class stepper, size_t lattice_dim>
int adaptive_routine(map<string, double> parameters, long seed = 0, string system="default",
                     string root = "", int count=0, double pre_T = -1.0) {
    // We try out the code for the brownian motion i would say
    // But we cannot use our old class system I think because there the whole system is already on a lattice
//
// Created by andi on 21.04.23.
    // But we can quickly write another system i guess
    double    dt_max = parameters["dt_max"];
    const double    T = parameters["T"];
    const double    J = parameters["J"];            // in anisotropic case Jx
    const double    Jy= parameters["Jy"];            // in anisotropic case Jy
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


    const double x0 = parameters["x0"] * sqrt(beta / 2.0);
    const double p0 = parameters["p0"] * sqrt(beta / 2.0);
    const           size_t n = lattice_dim * lattice_dim;

    double write_interval = parameters["end_t"] / nr_save_values;
    cout << "Starting Simulation for a " << lattice_dim << " by " << lattice_dim << " lattice until t = " << end_t << endl;

    // last time i didnt have to specify the dimensionality i think (in terms of (x, p) )
    const double D = T / eta;
    double mu = 0;
    double msd = 0;
    // file stuff
    string save_dir = root + "/" + to_string(parameters["T"]);
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

    stepper<gpu_state_type, thrust_algebra, thrust_operations, double, double> gpu_stepper(N * n, K, tol);

    // init and print initial state, we start in an equilibrium position, in the positive minimum
    gpu_state_type x(N * n, sqrt(beta / 2.0));

    // set the impulses to be zero
    thrust::fill(x.begin() + n, x.begin() + N * n, p0);
    // okay we overwrite this here
    if(parameters["random"] == 1.0) {
        // if random parameter is true we initialize high temperature random initial state
        chrono::milliseconds ms = chrono::duration_cast<chrono::milliseconds >(
                chrono::system_clock::now().time_since_epoch()
        );
        fill_init_values<gpu_state_type, n>(x, (float) x0, (float) p0, ms.count() % 10000);
    } else if (parameters["random"] == -1.0) {
        if (pre_T >= 0) {
            // else we need to read in the previous state
            string pre_dir_name = root + "/" + to_string(pre_T);
            string pre_name = pre_dir_name + "/" + to_string(count) + ".csv";
            ifstream pre_file = safe_read(pre_name, true);
            // and we are ready to read in the last file?
            double prev_T;
            double prev_t;
            vector<double> pre_lattice = readDoubleValuesAt(pre_file, -1, prev_T, prev_t);
            // we need to copy them into the gpu state type
            // just for loop?
            for(int i = 0; i < n; i++) {
                x[i] = pre_lattice[i];
            }
        } else if(pre_T < 0.0) {
            // if pre_T is smaller than zero that means that we didn't have a previous T so wi initialize random.
            // this is now code to check for runs that are already there
            cout << "checking for initial state in folder..." << endl;
            // listing the temp folders that are already inside
            vector<fs::path> temp_paths = list_dir_paths(root);
            // we check every folder name for the value and use the one that is closest to our actual temp
            string closest_T = findClosestDir(temp_paths, T);
            cout << "clostest folder to " << T << " already existing is " << closest_T << endl;
            if(closest_T != "None") {
                // now we list every csv file and take the one that is closest to i
                vector<fs::path> csv_files = list_csv_files(closest_T);
                print_vector(csv_files);
                string closest_i = findClosestStem(csv_files, count);
                cout << "clostest index to " << count << " already existing is " << closest_i << endl;
                if(closest_i != "None") {
                    fs::path pre_name = closest_i;
                    cout << "Trying to read " << pre_name << endl;
                    ifstream pre_file = safe_read(pre_name, true);
                    // and we are ready to read in the last file?
                    double prev_T;
                    double prev_t;
                    vector<double> pre_lattice = readDoubleValuesAt(pre_file, -1, prev_T, prev_t);
                    for(int i = 0; i < n; i++) {
                        x[i] = pre_lattice[i];
                    }
                }
            }
        }

    }

    for (int i = 0; i < n; i++) {
        mu += x[i];
        msd += x[i] * x[i];
    }
    double t = 0;
    double write_timepoint = 0;
    string obs_checkpoint = "Observing";



    // i THINK we will do the ugly way here until we understand polymorphism
    cout << "creating System " << system << endl;
    if (system == "constant") {
        checkpoint_timer obs_timer({obs_checkpoint});
        constant_bath<lattice_dim> gpu_system(T, eta, alpha, beta, J, seed);
        while(t < end_t) {
            // we need small stepsizes at the beginning to guarantee stability
            // but after some time, we can increase the stepsize
            // there should be a t that is equal to t_relax?

            gpu_stepper.do_step(gpu_system, x, dt_max, t);
            obs_timer.set_startpoint(obs_checkpoint);
            if (t >= write_timepoint) {
                Obs.writev2(gpu_system, x, t);
                write_timepoint += write_interval;
                cout << "current k = " << gpu_stepper.get_k() << " at t = " << t << endl;
            }
            obs_timer.set_endpoint(obs_checkpoint);
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
    } else if(system=="coulomb constant") {
        coulomb_constant<lattice_dim> gpu_system(T, eta, alpha, beta, J, seed);
        while(t < end_t) {
            // we need small stepsizes at the beginning to guarantee stability
            // but after some time, we can increase the stepsize
            // there should be a t that is equal to t_relax?

            gpu_stepper.do_step(gpu_system, x, dt_max, t);
            if (t >= write_timepoint) {
                Obs.writev2(gpu_system, x, t);
                write_timepoint += write_interval;
                cout << "current k = " << gpu_stepper.get_k() << " at t = " << t << endl;
            }
        }
    } else if(system=="anisotropic coulomb constant") {
        anisotropic_coulomb_constant<lattice_dim> ani_coulomb(T, eta, alpha, beta, J, Jy, seed);
        while(t < end_t) {
            // we need small stepsizes at the beginning to guarantee stability
            // but after some time, we can increase the stepsize
            // there should be a t that is equal to t_relax?

            gpu_stepper.do_step(ani_coulomb, x, dt_max, t);
            if (t >= write_timepoint) {
                Obs.writev2(ani_coulomb, x, t);
                write_timepoint += write_interval;
                cout << "current k = " << gpu_stepper.get_k() << " at t = " << t << endl;
            }
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

template <template<class, class, class, class, class> class stepper, size_t lattice_dim>
void repeat(map<string, double> parameters, int runs, long seed = 0, string system="default",
            string root="", int count=0, double pre_T = -1.0) {
    // seed is the seed for the random numbers so that we can have different random numbers per run
    if(runs == 0) {
        return;
    }

    cout << runs << " runs left" << endl;

    int steps = adaptive_routine<stepper, lattice_dim>(parameters, seed, system, root, count, pre_T);
    // how to get the number of steps that were done? let single calc routine return it?
    // or put it also into repeat
    repeat<stepper, lattice_dim>(parameters, runs - 1, seed + steps, system, root, count+1, pre_T);

}


void simple_temps_scan(string stepper = "adaptive", string system="constant") {
    // we always need to specify the lattice dim
    // const size_t* lattice_dim = &(size_t)adaptive_temp_scan_standard["lat_dim"];
    // const size_t lattice_dim = 100;
    string root = adaptive_tempscan_root;


    map<string, double> paras = adaptive_temp_scan_standard;
    // we do not use the fast forward here
    vector<double> T;
    if(paras["logspace"] == 1.0) {
        T = geomspace(paras["min_temp"],
                                          paras["max_temp"], (int)paras["nr_temps"] + 1);
    } else {
        T = linspace(paras["min_temp"],
                                          paras["max_temp"], (int)paras["nr_temps"] + 1);
    }
    // printing
    {
        print_vector(T);

        cout << "Starting Simulation for a " << lattice_dim << " by " << lattice_dim << " lattice unitl t = " << paras["end_t"] << endl;
        cout << "Stepsize is dt = " << paras["dt"] << endl;
    }

    // now cycling
    int i = 0;
    for(double temp : T) {
        paras["T"] = temp;
        // for every T we need to initalize a new system, but first i think we maybe should check our system?
        // for every T we need another directory

        double pre_T = -1.0;
        if(i > 0) {
            // if we already did a sim, eg if i >0, the previous temp can be obtained
            pre_T = T[i - 1];
        }

        cout << "Running repeat with following parameters:" << endl;
        // need a function here that takes the dirpath, looks if there are already files inside and
        // retunrs the highest number so that i can adjust the count
        int count = findHighestCSVNumber(root + "/" + to_string(temp)) + 1;
        cout << count << " Files already in Folder" << endl;
        int runs = (int)paras["repeat_nr"] - count;
        cout << "We fill the folder up to " << count << " + " << runs << " = " << count + runs << " realizations";
        printMap(paras);
        if(stepper == "adaptive") {
            repeat<euler_simple_adaptive, lattice_dim>(paras, runs, 0, system, root, count, pre_T);
        } else {
            repeat<euler_combined, lattice_dim>(paras, runs, 0, system, root, count, pre_T);
        }
        i++;
    }

    // print the root
    cout << root << endl;
}

int main() {
    simple_temps_scan("combined", "anisotropic coulomb constant");
    return 0;
}