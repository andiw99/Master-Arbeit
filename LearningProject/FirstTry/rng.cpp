//
// Created by weitze73 on 13.03.23.
//

#include "rng.h"
#include <random>
#include <functional>
#include <array>
#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

// problem: with adaptive step size you don't know how many values for the noise term you need
// we probably just generate a very large amount...

/**
 * Generates gaussian noise term values with the default random engine generator of the <random> package/header
 *
 * @param n         size of the grid of Silicon dimers, to calc how many values are needed per step
 * @param steps     number of steps of the numerical solution of the ODEs
 * @param sigma     stddev of the gaussian distribution
 * @param theta     temperature of the bath
 * @param eta       damping term of langevin equation
 */
template <size_t n, size_t steps>
array<array<double, n*n>, steps> generate(double sigma, double theta, double eta) {
    default_random_engine generator;
    // TODO: mean is zero? Or is the mean theta * eta?
    // TODO: i guess the direction of the bump should not be biased, so let the mean be 0 for now
    normal_distribution<double> distribution(0, sigma);
    auto dice = bind(distribution, generator);
    /*  We probably create an array first which we will then save/write to a file. Array of the form
     *      Step    Index 00    01      02  ..  0n      10      11  ..  nn
     *       1            0.2   0.1     ....     more random numbers
     *       ..
     *       steps
     */
    // saveable in memory? 1000 doubles per Step, 10000 steps? 10 Million doubles 100 Million Bytes, 100 MB EZ!
    // init array
    // TODO i don't really get how to use dynamic memory
    /*
    array<array<double, n*n>, steps> * rngs;
    rngs = new array<array<double, n*n>, steps>;
    */
    array<array<double, n*n>, steps> rngs{};
    // nested loops
    for(int step=0; step < steps; step++) {
        for(int i = 0; i < n * n; i++) {
            rngs[step][i] = dice();
        }
    }
    return rngs;
}

/**
 * Writes the generated random numbers to a file
 * @tparam n        lattice size
 * @tparam steps    number of steps of the numerical solution of ODEs
 * @param rngs      generated random numbers
 */
template <size_t n, size_t steps>
void write_to_csv(array<array<double, n*n>, steps> rngs) {
    ofstream out("rngs.csv");
    for(auto& row : rngs) {
        for(auto col : row) {
            out << col << ", ";
        }
        out << "\n";
    }
}

/**
 * new version to generate random numbers
 * @param n         lattice size
 * @param steps     number of steps of the numerical solution of ODEs
 */
void generate_and_write(size_t n, size_t steps) {
    default_random_engine generator;
    normal_distribution<double> distribution(0, 1);
    auto dice = bind(distribution, generator);
    ofstream out("rngs.csv");
    for(int j = 0; j < steps; j++) {
        for(int i = 0; i < n*n; i++) {
            out << dice() << ", ";
        }
        out << "\n";
    }
}

int main() {
    // parameters
    const size_t n = 20;          // lattice size
    const size_t steps = 10000;     // preferably larger number than steps are needed to be done in numerical solution of ODEs
    // TODO hier workaround weil du irgendwie keine größere Matrix als 400 * 1000 im memory haben kannst
    const size_t block_size = 1000;
    const size_t iterations = steps / block_size;
    double sigma = 1, theta = 0.1, eta = 1.0;
    generate_and_write(n, steps);
    /*
    for(int iter = 0; iter < iterations; iter++) {
        array<array<double, n*n>, block_size> rngs = generate<n, block_size>(sigma, theta, eta);
        write_to_csv<n, block_size>(rngs);
    }
     */
    return 0;
}