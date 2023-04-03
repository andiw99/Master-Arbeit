//
// Created by andi on 25.03.23.
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
typedef  vector<double> state_type;

class brownian_motion {
// unser state_type ist einfach 2d Vektor, lass uns ublas matrix nehmen? Oder normalen vector?
    double eta;
    //static const auto dice = bind(normal_distribution<double>(0, 1), default_random_engine());
    // _Bind<normal_distribution(linear_congruential_engine<unsigned long, 16807, 0, 2147483647>)> dice;
    double T;
    double dt;
    static inline auto normal_dice = bind(normal_distribution<double>(0, 1), default_random_engine(12));

    /*
    normal_distribution<double> normal_dist;
    default_random_engine generator;
    */
    void stochastic_euler_method(state_type &x, const state_type &dxdt, const state_type &theta) {
        // huh we dont have vector operations, maybe implement them
        // TODO for now we just do the step entry wise
        x[0] = x[0] + dxdt[0] * dt + theta[0] * sqrt(dt);
        /*
        cout << "normal: ";
        cout << theta[1] << endl;
        cout << "dvdt: ";
        cout << dxdt[1] << endl;
        */
         x[1] = x[1] + dxdt[1] * dt + theta[1] * sqrt(dt);
    }

    void brownian_step(const state_type &x, state_type &dxdt, state_type &theta, double t){
        dxdt[0] = x[1];
        // dxdt[1] = 0;
        dxdt[1] = (-eta) * x[1];
        theta[0] = 0;
        theta[1] = sqrt(2 * eta * T) * normal_dice();
        /*
        cout << "normal: ";
        cout << theta[1] << endl;
        */
         }

public:
    void run(int steps, string name="Brownian Motion Data") {
        // runs brownian motion with stepsize dt for steps steps
        // initialize
        state_type x = vector<double>(2);
        state_type dxdt = vector<double>(2);
        state_type theta = vector<double>(2);
        // set initial time to 0
        double t = 0;

        std::ofstream myfile;
        myfile.open (name + ".csv");

        // average value
        double mu = 0;
        for(int i = 0; i < steps; i++) {
            brownian_step(x, dxdt, theta, t);
            stochastic_euler_method(x, dxdt, theta);

            // write to file
            myfile << t << ",";
            myfile << x[0] << ",";
            myfile << x[1] << ", \n";

            mu += x[0];
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
    brownian_motion(double eta_val, double dt_val, double T_val) {
        eta = eta_val;
        dt = dt_val;
        T = T_val;
    }
};

int main() {
    std::cout << "Test" << endl;
    double eta = 10;
    double T = 1;
    double dt = 0.001;
    int steps = 100000;

    brownian_motion suite = brownian_motion(eta, dt, T);
    // number of runs
    int runs = 10;
    for(int i = 0; i < runs; i++) {
        string name = "Brownian Motion Data/" + to_string(i);
        suite.run(steps, name);
    }


    return 0;
}

