#include "main.cuh"
#include "systems.cuh"


int main() {
    // We try out the code for the brownian motion i would say
    // But we cannot use our old class system I think because there the whole system is already on a lattice
//
// Created by andi on 21.04.23.
    // But we can quickly write another system i guess

    const int steps = 100000;
    const double dt = 0.005;
    const double T = 70;
    const double J = 50;
    const double alpha = 5;
    const double beta = 10;
    const double tau = 10;
    const double eta = 5;
    const int nr_save_values = 150;
    size_t write_every = steps / nr_save_values;
    const size_t lattice_dim = 250;
    // system size
    const size_t n = lattice_dim * lattice_dim;
    // DGLs per lattice site
    const size_t N = 2;
    const double x0 = 0.0;
    const double p0 = 0.0;

    // last time i didnt have to specify the dimensionality i think (in terms of (x, p) )
    const double D = T / eta;
    double theo_msd = 2 * D * dt * steps;
    double mu = 0;
    double msd = 0;
    // file stuff
    string storage_root = "../../../Generated content/DomainSize/";
    string dir_name = create_tree_name(eta, T, dt, n, alpha, beta, J, tau,
                                       storage_root);
    string name = dir_name + "/" + to_string(lattice_dim);
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
    gpu_bath<lattice_dim> gpu_system(T, eta, alpha, beta, J, tau);
    // gpu_oscillator_chain<lattice_dim> gpu_system(T, eta, alpha);
    // this state type sets x and p to be x0, meaning 100 in our case.
    gpu_state_type x(N * n, x0);
    // set the impulses to be zero
    thrust::fill(x.begin() + n, x.begin() + N * n, p0);
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
    for(int i = 0; i <= n; i++) {
        mu += x[i];
        msd += x[i] * x[i];
    }
    mu /= n;
    msd /= n;

    cout << "mu = " << mu << endl;
    cout << "msd = " << msd << "   theo value msd = " << theo_msd << endl;


    // add to msd and mu
    return 0;
}
