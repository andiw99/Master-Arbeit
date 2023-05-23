//
// Created by andi on 14.05.23.
//
#include "main.cuh"
#include "systems.cuh"
#include <numeric>
#include <chrono>
#include <map>


struct rand_init_values
{
    double mu, sigma, ampl;

    __host__ __device__
    rand_init_values(double ampl, double mu = 0.0, double sigma = 1.0) : ampl(ampl), mu(mu), sigma(sigma) {};

    __host__ __device__
    float operator()(const unsigned int ind) const
    {
        thrust::default_random_engine rng;
        thrust::normal_distribution<double> dist(mu, sigma);
        rng.discard(ind);

        // dist * ampl zu returnen ist wie... aus dist mit std ampl zu ziehen: b * N(m, o) = N(m, b*o)

        return ampl * dist(rng);
    }
};

template <size_t n>
void fill_init_values(thrust::device_vector<double>& state, float x0, float p0, int run = 0, double mu=0, double sigma=1) {
    cout << "do i get called?  " << x0 << endl;
    thrust::counting_iterator<size_t> index_sequence_begin(run * state.size());
    // thrust::fill(theta.begin(), theta.begin() + n, 0);
    // n is system size
    // fill the starting positions
    thrust::transform(index_sequence_begin,
                      index_sequence_begin + n,
                      state.begin(),
                      rand_init_values(x0, mu, sigma));
    // fill starting impulses
    thrust::transform(index_sequence_begin + n,
                      index_sequence_begin + 2*n,
                      state.begin() + n,
                      rand_init_values(p0, mu, sigma));
}


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

    // initialize the system...
    // ugly af
    // TODO this won't work again probalby?
/*    System<lattice_dim>* gpu_system;
    if(system == "default") {
        gpu_system = new gpu_bath<lattice_dim>(T, eta, alpha, beta, J, tau, seed);
    } else if (system == "constant") {
        gpu_system = new constant_bath<lattice_dim>(T, eta, alpha, beta, J, seed);
    } else {
        throw runtime_error("invalid system name");
    }*/


    // gpu_bath(const double T, const double eta, const double alpha, const double beta, const double J, const double tau);
    // constant_bath<lattice_dim> gpu_system(T, eta, alpha, beta, J);
    // gpu_oscillator_chain<lattice_dim> gpu_system(T, eta, alpha);
    // this state type sets x and p to be x0, meaning 100 in our case.
    gpu_state_type x(N * n, x0);
    // set the impulses to be zero
    thrust::fill(x.begin() + n, x.begin() + N * n, p0);
    // okay we overwrite this here
    fill_init_values<n>(x, (float)x0, (float)p0);
    for(int i = 0; i < n; i++) {
        mu += x[i];
        msd += x[i] * x[i];
    }
    cout << "Initial values:" << endl;
    cout << mu / (n) << endl;
    cout << msd / (n) << endl;
    mu = 0;
    msd = 0;

    auto start = chrono::high_resolution_clock::now();
    double t = 0;


    // i THINK we will do the ugly way here until we understand polymorphism
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
    } else if(system == "quadratic chain") {
        quadratic_chain<lattice_dim> gpu_system(T, eta, J, seed);
        for( size_t i=0 ; i<steps ; ++i ) {
            gpu_stepper.do_step(gpu_system, x, dt, t);
            Obs(gpu_system, x, t);
            t += dt;
        }
        // calculate the energy and print it
        double E = gpu_system.calc_energy(x);
        cout << "Energy of the System: " << E << endl;
        cout << "Theoretical Energy: " << (double)N * T << endl;
    } else {
        gpu_bath<lattice_dim> gpu_system(T, eta, alpha, beta, J, tau, seed);
        for( size_t i=0 ; i<steps ; ++i ) {
            gpu_stepper.do_step(gpu_system, x, dt, t);
            Obs(gpu_system, x, t);
            t += dt;
        }
    }
/*    for( size_t i=0 ; i<steps ; ++i ) {
        gpu_stepper.do_step(*gpu_system, x, dt, t);
        Obs(*gpu_system, x, t);
        t += dt;
    }*/

    write_parameters(parafile, eta, T, dt, n, alpha, beta, J, tau);

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

template <size_t lattice_dim>
int single_calc_routine(long seed = 0, string system="default", string save_dir="") {
    map<string, double> paras;

    // Adding key-value pairs to the map
    paras["steps"] = (double) 100000;
    paras["dt"] = 0.001;
    paras["T"] = 30;
    paras["J"] = 10;
    paras["alpha"] = 5;
    paras["beta"] = 10;
    paras["tau"] = 10;
    paras["eta"] = 1.2;
    paras["N"] = 32;
    paras["nr_save_values"] = 250;
    paras["x0"] = 0;
    paras["p0"] = 0;

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

    string root = "../../Generated content/New Scan/";
    root = (root_val.empty()) ? root : root_val;
    // make sure root exists

    int             end_t = 2000;
    int             t_relax = 50;                     // approximate time where system is relaxed, probably math to approximate?
    double dt_max =  0.001;                             // max dt of 0.005 was to high
    double dt_start = 0.0001;
    double dt = dt_start;
    int steps = (int)(t_relax / dt_start + (end_t - t_relax) / dt_max);
    const int       nr_temps = 40;
    const double    J = 2;
    const double    alpha = 1;
    const double    beta = 10;
    const double    tau = 10;
    const double    eta = 0.2;
    const int       nr_save_values = 32;
    const size_t    lattice_dim = 100;
    const size_t    n = lattice_dim * lattice_dim;
    const size_t    N = 2;
    int             repeat_nr = 5;

    const vector<double> T = linspace(5.0, 7.0, nr_temps + 1);

    steps = (steps_val == 0) ? steps : steps_val;
    end_t = (end_t_val == 0) ? end_t : end_t_val;

    map<string, double> paras;
    paras["steps"] = steps;
    paras["end_t"] = end_t;
    paras["t_relax"] = t_relax;
    paras["dt"] = dt_start;
    paras["dt_max"] = dt_max;
    paras["J"] = J;
    paras["alpha"] = alpha;
    paras["beta"] = beta;
    paras["tau"] = tau;
    paras["eta"] = eta;
    paras["lattice_dim"] = lattice_dim;
    paras["N"] = N;
    paras["nr_save_values"] = nr_save_values;


    print_vector(T);


    cout << "Starting Simulation for a " << lattice_dim << " by " << lattice_dim << " lattice for " << steps << " steps." << endl;
    cout << "Stepsize is dt = " << dt << endl;
    const double x0 = 8.0;
    const double p0 = 8.0;
    paras["x0"] = x0;
    paras["p0"] = p0;

    create_dir(root);

    // stepper will be the same for all temps
    euler_mayurama_stepper<gpu_state_type, thrust_algebra, thrust_operations > gpu_stepper(N * n);

    // keep track at which
    int count = 0;
    for(double temp : T) {
        paras["T"] = temp;
        // for every T we need to initalize a new system, but first i think we maybe should check our system?
        // for every T we need another directory
        string dirpath = root + "/" + to_string(temp);
        // we could repeat?
        // we need a new observer with a file name
        repeat<lattice_dim>(paras, repeat_nr, 0, "constant", dirpath);
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



    // single_calc_routine();
    // add to msd and mu
    return 0;
}
