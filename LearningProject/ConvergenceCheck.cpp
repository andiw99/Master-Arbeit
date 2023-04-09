//
// Created by andi on 06.04.23.
//
#include "Functions and Classes.h"

int main() {
    // we wont to check whether the stochastic euler method converges.
    // for brownian motion in harmonic oscillators, the
    double eta = 5;
    double T = 100; // 10
    double dt = 0.01;
    int steps = 100000;
    int n = 20;
    double alpha = 1;   // 1
    double beta  = 5;
    double J = 10;
    double tau = 2;
    double starting_t = -10;
    int runs = 10;

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
    for(int i = 0; i < runs; i++) {
        suite.run(steps, file, false);
    }
    write_parameters(file, eta, T, dt, n, alpha, beta, J, tau);
    file.close();




    return 0;
}