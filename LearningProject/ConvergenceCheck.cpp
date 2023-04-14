//
// Created by andi on 06.04.23.
//
#include "Helpfunctions and Classes.h"
#include "Systems.h"
#include "Solvers.h"

int main() {
    // we wont to check whether the stochastic euler method converges.
    // for brownian motion in harmonic oscillators, the
    double eta = 5;
    double T = 100; // 10
    double passed_time = 100;
    int steps = 1000000;
    double dt = passed_time/steps;
    int n = 20;
    double alpha = 1;   // 1
    double beta  = 5;
    double J = 10;
    double tau = 2;
    double starting_t = -10;
    int runs = 10;

    double msd_harmonic_osc = T / alpha;
    double mu_harmonic_osc = 0;

    string storage_root = "../../Generated content/";

    vector<vector<double>> x0 =
            vector<vector<double>>(n, vector<double>(n, 0));
    vector<vector<double>> v0 =
            vector<vector<double>>(n, vector<double>(n, 0));
    /*
    full_interaction_system specific_system =
            full_interaction_system(eta, T, n, x0,
                                    v0,alpha,beta, J,
                                    tau);
    */
     harmonic_system HarmonicSystem =
            harmonic_system(eta, T, n, x0, v0, alpha);
    solver suite =
            solver(dt, &HarmonicSystem, starting_t);
    string dir_name = create_directory(eta, T, dt, n, alpha, beta, J, tau,
                                       storage_root);
    // since dir_name is per construction
    // a directory we dont have to check if its a directory?
    string name = dir_name + "/" + to_string(count_files(dir_name));
    std::ofstream file;
    file.open (name + ".csv");
    cout << "Evaluation time t = " << dt * steps << endl;
    cout << "Step size dt = " << dt << endl;
    cout << "n = " << n << endl;
    cout << "nr of runs: " << runs << endl;
    double error_mu = 0;
    double error_msd = 0;
    for(int i = 0; i < runs; i++) {
        state_type x = suite.run(steps, file, false);
        double mu = suite.calc_mu(x);
        double msd = suite.calc_msd(x, 0);
        cout << "supposed value: mu = " << mu_harmonic_osc << "   actual value: mu = " << mu << endl;
        cout << "supposed value: msd = " << msd_harmonic_osc << "   actual value: msd = " << msd << endl;
        error_mu += abs(mu_harmonic_osc - mu);
        error_msd += abs(msd_harmonic_osc - msd);
    }
    error_mu /= runs;
    error_msd /= runs;
    cout << "Average deviation dmu = " << error_mu << endl;
    cout << "Average deviation dmsd = " << error_msd << endl;
    write_parameters(file, eta, T, dt, n, alpha, beta, J, tau);
    file.close();




    return 0;
}