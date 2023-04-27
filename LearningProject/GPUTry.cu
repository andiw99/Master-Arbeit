//
// Created by andi on 21.04.23.
//

#include "GPUTry.h"
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>

    const int steps = 50000;
    const double dt = 0.001;
    const double eta = 5;
    const double T = 100;

    // H has storage for 4 integers
    thrust::host_vector<int> H(4);

    // initialize individual elements
    H[0] = 14;
    H[1] = 20;
    H[2] = 38;
    H[3] = 46;

    // H.size() returns the size of vector H
    std::cout << "H has size " << H.size() << std::endl;

    // print contents of H
    for(int i = 0; i < H.size(); i++)
        std::cout << "H[" << i << "] = " << H[i] << std::endl;

    // resize H
    H.resize(2);

    std::cout << "H now has size " << H.size() << std::endl;

    // Copy host_vector H to device_vector D
    thrust::device_vector<int> D = H;

    // elements of D can be modified
    D[0] = 99;
    D[1] = 88;

    // print contents of D
    for(int i = 0; i < D.size(); i++)
        std::cout << "D[" << i << "] = " << D[i] << std::endl;

    // last time i didnt have to specify the dimensionality i think (in terms of (x, p) )

    /*
     *

    const size_t N = 2;
    euler_mayurama_stepper<state_type, container_algebra, default_operations> stepper(N);
    brownian_particel system(eta, T);

    // I guess we directly have to check whether the distribution parameters are still the same.
    const size_t runs = 10000;
    double mu = 0;
    double msd = 0;
    double D = T / eta;
    double theo_msd = 2 * D * dt * steps;
    double theo_mu = 0;
    for(size_t i = 0; i < runs; i++) {
        // init the inital values for every run
        state_type x(N, 0.0);
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
    cout << "msd = " << msd << "   theo value msd = " << theo_msd << endl;
     */
    return 0;
}
