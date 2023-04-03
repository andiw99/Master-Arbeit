//
// Created by andi on 28.03.23.
//

#include "RightAndDown.h"
//
// Created by andi on 27.03.23.
//

#include "main.h"
#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <random>
#include <fstream>

// using namespaces!
using namespace std;
using namespace boost::numeric::odeint;
// new state is composed of every
typedef vector<double> entry_type;
typedef  vector<vector<entry_type>> state_type;

class lattice_model {
// unser state_type ist einfach 2d Vektor, lass uns ublas matrix nehmen? Oder normalen vector?
    double eta;
    //static const auto dice = bind(normal_distribution<double>(0, 1), default_random_engine());
    // _Bind<normal_distribution(linear_congruential_engine<unsigned long, 16807, 0, 2147483647>)> dice;
    double T;
    double dt;
    static inline auto normal_dice = bind(normal_distribution<double>(0, 1), default_random_engine(12));
    double alpha;
    double beta;
    double J;
    vector<vector<double>> x0;
    vector<vector<double>> v0;
    int n;      // size of the lattice
    double tau;
    double tmin;       // starting time

    /*
    normal_distribution<double> normal_dist;
    default_random_engine generator;
    */
    void stochastic_euler_method(state_type &x, const state_type &dxdt, const state_type &theta) {
        // huh we dont have vector operations, maybe implement them
        // We use the state_type here, so we need to cycle through every lattice
        for(int i = 0; i < n; i++) {
            for(int j=0; j < n; j++) {
                x[i][j][0] = x[i][j][0] + dxdt[i][j][0] * dt + theta[i][j][0] * sqrt(dt);
                x[i][j][1] = x[i][j][1] + dxdt[i][j][1] * dt + theta[i][j][1] * sqrt(dt);
            }
        }
    }

    double dVsdq(double q, double t) {
        return alpha * (2 * q * q * q - beta * eps(t) * q);
    }

    double eps(double t) {
        // TODO normally with time scale, but is this important here?
        return min(t/tau, 1.0);
    }

    double dVidq(double qij, double qijp1, double qip1j) {
        return J * (sin(qij - qijp1) + sin(qij - qip1j));
    }

    double dVdq(double qij, double qijp1, double qip1j, double t) {
        return dVsdq(qij, t) + dVidq(qij, qijp1, qip1j);
    }

    void brownian_step(const state_type &x, state_type &dxdt, state_type &theta, double t){
        for(int i = 0; i < n; i++) {
            for(int j=0; j < n; j++) {
                dxdt[i][j][0] = x[i][j][1];
                // TODO Problem when calculating this are the border terms but here we can just use a modulo for now?
                dxdt[i][j][1] = (-eta) * x[i][j][1] - dVdq(x[i][j][0], x[i][(j + 1) % n][0], x[(i+1) % n][j][0], t);
                theta[i][j][0] = 0;
                theta[i][j][1] = sqrt(2 * eta * T) * normal_dice();
            }
        }

    }

public:
    void run(int steps, string name="Brownian Motion in Potential Data") {
        // runs brownian motion with stepsize dt for steps steps
        // initialize
        state_type x = state_type (n, vector<entry_type>(n, entry_type(2, 0)));
        // set initial values
        // naive setting because i don't know c++ to much atm and it only runs once when initializing
        for(int i = 0; i < n; i++) {
            for(int j=0; j < n; j++) {
                x[i][j][0] = x0[i][j];
                x[i][j][1] = v0[i][j];

            }
        }
        state_type dxdt = state_type (n, vector<entry_type>(n, entry_type(2, 0)));
        state_type theta = state_type (n, vector<entry_type>(n, entry_type(2, 0)));

        // set ininital time to starting time
        double t = tmin;

        std::ofstream myfile;
        myfile.open (name + ".csv");

        // average value
        double mu = 0;
        for(int k = 0; k < steps; k++) {
            brownian_step(x, dxdt, theta, t);
            stochastic_euler_method(x, dxdt, theta);

            // write to file
            // everything?
            myfile << t << ",";
            for(int i = 0; i < n; i++) {
                for(int j=0; j < n; j++) {
                    myfile << x[i][j][0] << ",";
                }
            }
            for(int i = 0; i < n; i++) {
                for(int j=0; j < n; j++) {
                    myfile << x[i][j][1] << ",";
                }
            }

            // Zeilenumbruch
            myfile << "\n";

            t += dt;
        }
        myfile.close();
        // since t is now a class variable, we h
    }
    lattice_model(double eta_val, double dt_val, double T_val,  int n_val,
                  vector<vector<double>> x0_val,
                  vector<vector<double>> v0_val,
                  double alpha_val, double beta_val, double J_val, double tau_val, double t_val = -5) {
        eta = eta_val;
        dt = dt_val;
        T = T_val;
        x0 = x0_val;
        v0 = v0_val;
        n = n_val;
        alpha = alpha_val;
        beta = beta_val;
        J = J_val;
        tau = tau_val;
        tmin = t_val;
    }
};

int main() {
    double eta = 5;
    double T = 0.01; // 10
    double dt = 0.01;
    int steps = 2000;
    int n = 20;
    double alpha = 1;   // 1
    double beta  = 10;
    double J = 10;
    double tau = 10;
    double starting_t = -10;

    vector<vector<double>> x0 = vector<vector<double>>(n, vector<double>(n, 0));
    vector<vector<double>> v0 = vector<vector<double>>(n, vector<double>(n, 0));

    lattice_model suite = lattice_model(eta, dt, T,n , x0, v0,
                                        alpha, beta, J, tau, starting_t);
    // number of runs
    int runs = 1;
    for(int i = 0; i < runs; i++) {
        string name = "Brownian Motion Two Side Interaction/" + to_string(i);
        suite.run(steps, name);
    }


    return 0;
}

