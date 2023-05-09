#include "main.cuh"
#include "systems.cuh"
#include <numeric>
#include <chrono>


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


int single_calc_routine(long seed = 0) {
    // We try out the code for the brownian motion i would say
    // But we cannot use our old class system I think because there the whole system is already on a lattice
//
// Created by andi on 21.04.23.
    // But we can quickly write another system i guess

    const int steps = 100000;
    const double dt = 0.001;
    const double T = 30;
    const double J = 50;
    const double alpha = 5;
    const double beta = 10;
    const double tau = 10;
    const double eta = 1.2;
    const int nr_save_values = 32;
    size_t write_every = steps / nr_save_values;
    const size_t lattice_dim = 250;
    // system size
    const size_t n = lattice_dim * lattice_dim;
    // DGLs per lattice site
    const size_t N = 2;

    cout << "Starting Simulation for a " << lattice_dim << " by " << lattice_dim << " lattice for " << steps << "steps." << endl;
    const double x0 = 8.0;
    const double p0 = 8.0;

    // last time i didnt have to specify the dimensionality i think (in terms of (x, p) )
    const double D = T / eta;
    double theo_msd = 2 * D * dt * steps;
    double mu = 0;
    double msd = 0;
    // file stuff
    string storage_root = "../../../Generated content/Constant Bath/";
    string dir_name = create_tree_name(eta, T, dt, n, alpha, beta, J, tau,
                                       storage_root);
    string name = dir_name + "/" + to_string(seed/steps);
    ofstream file;
    ofstream parafile;
    file.open(name + ".csv");
    parafile.open(name + ".txt");
/*    observer base_obs = observer();

    euler_mayurama_stepper<state_type, container_algebra, default_operations> stepper(2, &base_obs);
    brownian_particel system(eta2, T);

    // I guess we directly have to check whether the distribution parameters are still the same.
    const size_t runs = 10;



    for(size_t i = 0; i < runs; i++) {
        // init the inital values for every run
        state_type x(2, 0.0);
        for( size_t n=0 ; n<steps ; ++n ) {
            stepper.do_step(system, x, dt, n*dt);
            // cout << n*dt << " ";
            // cout << x[0] << " " << x[1] << endl;
        }
        // add to msd and mu
        mu += x[0];
        msd += x[0] * x[0];
    }
    // averaging
    mu /= runs;
    msd /= runs;

    cout << "mu = " << mu << endl;
    cout << "msd = " << msd << "   theo value msd = " << theo_msd << endl;*/



    // init observer

    bath_observer Obs(file, write_every);

    typedef thrust::device_vector<double> gpu_state_type;
    euler_mayurama_stepper<gpu_state_type, thrust_algebra, thrust_operations > gpu_stepper(N * n);

    // initialize the system...
    // gpu_bath(const double T, const double eta, const double alpha, const double beta, const double J, const double tau);
    gpu_bath<lattice_dim> gpu_system(T, eta, alpha, beta, J, tau, seed);
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
        cout << x[i] << endl;
    }
    cout << "Initial values:" << endl;
    cout << mu / (n) << endl;
    cout << msd / (n) << endl;
    mu = 0;
    exit(0);
    msd = 0;
    /*


    // We initialize a system of size 50000... what does that even mean?
    gpu_brownian_system gpu_system(eta, T, n);

    // so this state type only has N=2 values, which are set to 0. Maybe we need to initialize N * n values?
    gpu_state_type x(N * n, 0.0);
      */

    auto start = chrono::high_resolution_clock::now();
    double t = 0;

    for( size_t i=0 ; i<steps ; ++i ) {
        gpu_stepper.do_step(gpu_system, x, dt, t);
        // ... why don't we just write here? wouldn't that be faster?
        // i guess we do that
        Obs(gpu_system, x, t);
        t += dt;
        // cout << n*dt << " ";
        // cout << x[0] << " " << x[1] << endl;
    }

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


void scan_temps_routine(const int steps_val = 0, const int end_t_val = 0, const string& root_val = "", double mu = 0, double sigma=1) {

    string root = "../../Generated content/Relax Scan Underdamped Detailed/";
    root = (root_val.empty()) ? root : root_val;
    // make sure root exists

    int steps = 3000000;
    steps = (steps_val == 0) ? steps : steps_val;
    int end_t = 100;
    end_t = (end_t_val == 0) ? end_t : end_t_val;
    const double dt = (double)end_t / steps;
    // select Temps
    const int nr_temps = 75;
    const vector<double> T = linspace(7.5, 45.0, nr_temps + 1);
    print_vector(T);
    const double J = 50;
    const double alpha = 1;
    const double beta = 10;
    const double tau = 10;
    const double eta = 0.2;
    const int nr_save_values = 32;
    size_t write_every = steps / nr_save_values;
    const size_t lattice_dim = 100;
    // system size
    const size_t n = lattice_dim * lattice_dim;
    // DGLs per lattice site
    const size_t N = 2;

    cout << "Starting Simulation for a " << lattice_dim << " by " << lattice_dim << " lattice for " << steps << " steps." << endl;
    cout << "Stepsize is dt = " << dt << endl;
    const double x0 = 8.0;
    const double p0 = 8.0;

    create_dir(root);

    // stepper will be the same for all temps
    euler_mayurama_stepper<gpu_state_type, thrust_algebra, thrust_operations > gpu_stepper(N * n);

    // keep track at which
    int count = 0;
    for(double temp : T) {
        // for every T we need to initalize a new system, but first i think we maybe should check our system?
        // we need a new observer with a file name
        count++;
        ofstream file;
        file.open(root + to_string(temp) + to_string(dt) + to_string(lattice_dim) + ".csv");
        bath_observer Obs(file, write_every);
        // change the system here
        // gpu_oscillator_chain<lattice_dim> gpu_system(temp, eta, alpha);
        constant_bath<lattice_dim> gpu_system(temp, eta, alpha, beta, J);
        // init state with random numbers
        gpu_state_type x(N * n, 0.0);
        fill_init_values<n>(x, (float)x0, (float)p0, count, mu, sigma);
        // okay we overwrite this here for now and fill with constant values
        // thrust::fill(x.begin(), x.begin() + n, x0);
        // thrust::fill(x.begin() + n, x.begin() + N * n, p0);

        // check real quick that we have different random numbers in every run
        double mu = 0;
        double mu_p = 0;
        double msd = 0;
        double msd_p = 0;
        for(int i = 0; i < n; i++) {
            mu += x[i];
            mu_p += x[i + n];
            msd += x[i] * x[i];
            msd_p += x[i + n] * x[i + n];
        }
        mu /= n;
        mu_p /= n;
        msd /= n;
        msd_p /= n;

        cout << "Starting run " << count << " / " << T.size() << ".\n";
        cout << "Initial state was initialized with mu = " << mu << " and msd = " << msd << endl;
        cout << "And mu_p = " << mu_p << " and msd_p = " << msd_p << endl;
        cout << "Theoretical msd should be " << x0 * x0 << endl;


        double t = 0;

        for( size_t i=0 ; i<steps ; ++i ) {
            gpu_stepper.do_step(gpu_system, x, dt, t);
            // ... why don't we just write here? wouldn't that be faster?
            // i guess we do that
            Obs(gpu_system, x, t);
            t += dt;
            // cout << n*dt << " ";
            // cout << x[0] << " " << x[1] << endl;
        }


        ofstream parafile;
        parafile.open(root + to_string(temp) + to_string(dt) + to_string(lattice_dim)  + ".txt");
        write_parameters(parafile, eta, temp, dt, n, alpha, beta, J, tau);
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

void repeat(int runs, long seed = 0) {
    // seed is the seed for the random numbers so that we can have different random numbers per run
    if(runs == 0) {
        return;
    }

    int steps = single_calc_routine(seed);
    // how to get the number of steps that were done? let single calc routine return it?
    // or put it also into repeat
    repeat(runs - 1, seed + steps);

}


int main() {
    scan_temps_routine();
    // convergence_check_oscillatorchain();



    // single_calc_routine();
    // add to msd and mu
    return 0;
}
