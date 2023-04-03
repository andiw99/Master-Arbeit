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
    double w;
    vector<vector<double>> x0;
    vector<vector<double>> v0;
    int n;      // size of the lattice

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

    double dVdx(double x, double t) {
        return w * w * x * x * x - 3 * min(t, 1.0) * x;
    }

    void brownian_step(const state_type &x, state_type &dxdt, state_type &theta, double t){
        for(int i = 0; i < n; i++) {
            for(int j=0; j < n; j++) {
                dxdt[i][j][0] = x[i][j][1];
                dxdt[i][j][1] = (-eta) * x[i][j][1] - dVdx(x[i][j][0], t);     // i could basically add interaction here
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
        // set initial time to be negative
        double t = -5 ;

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


            //mu += x[0];
            t += dt;
            // print position
            /*
            cout << "pos: ";
            cout << x[0] << endl;
            cout << "vel: ";
            cout << x[1] << endl;
             */
        }
        myfile.close();
        mu /= steps;
        cout << endl;
        cout << mu << endl;
    }
    lattice_model(double eta_val, double dt_val, double T_val, double w_val,  int n_val,
                  vector<vector<double>> x0_val,
                  vector<vector<double>> v0_val) {
        eta = eta_val;
        dt = dt_val;
        T = T_val;
        w = w_val;
        x0 = x0_val;
        v0 = v0_val;
        n = n_val;
    }
};

int main() {
    double eta = 1;
    double T = 0.05; // 10
    double dt = 0.008;
    double w = 2;   // 1
    int steps = 1500;
    int n = 20;

    vector<vector<double>> x0 = vector<vector<double>>(n, vector<double>(n, 0));
    vector<vector<double>> v0 = vector<vector<double>>(n, vector<double>(n, 0));

    lattice_model suite = lattice_model(eta, dt, T, w,n , x0, v0);
    // number of runs
    int runs = 1;
    for(int i = 0; i < runs; i++) {
        string name = "Brownian Motion Stochastic Lattice/" + to_string(i);
        suite.run(steps, name);
    }


    return 0;
}

