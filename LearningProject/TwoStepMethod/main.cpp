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
typedef  vector<vector<double>> state_type;

class two_step {
// unser state_type ist einfach 2d Vektor, lass uns ublas matrix nehmen? Oder normalen vector?
    double eta;
    //static const auto dice = bind(normal_distribution<double>(0, 1), default_random_engine());
    // _Bind<normal_distribution(linear_congruential_engine<unsigned long, 16807, 0, 2147483647>)> dice;
    double T;
    double dt;
    static inline auto normal_dice = bind(normal_distribution<double>(0, 1), default_random_engine(12));
    double w;
    double x0;
    double v0;


    void stochastic_two_step(state_type &x, const state_type &W, const state_type &theta, const state_type &A) {
        // q at time t
        double qt = x[0][1];
        double vt = x[1][1];
        // BGL for postition and velocity
        x[0][1] = x[0][1] + vt * dt;
        x[1][1] = vt + A[0][1] * dt - A[1][0] * theta[1][0] * dt * sqrt(dt) + W[1][0] + W[1][1];

        // set past values
        x[0][0] = qt;
        x[1][0] = vt;

    }

    double dVdq(double q) {
        return w * w * q;
    }

    void brownian_step(const state_type &x, state_type &A, state_type &theta, state_type &W, double t){
        // x is not altered here

        // A
        double At = A[0][1];
        double dAdqt = A[1][1];
        // calc new Values
        A[0][1] = A_func(x[1][1], x[0][1]);
        A[1][1] = A_func(x[1][1], x[0][1]);
        // replace past values
        A[0][0] = At;
        A[1][0] = dAdqt;

        //Theta
        // TODO das muss nicht jeden step gesetzt werden
        theta[0][0] = 0;
        theta[0][1] = 0;
        // theta at time t
        double tt = theta[1][1];
        // roll new theta
        theta[1][1] = sqrt(2 * eta * T) * normal_dice();
        // replace old theta
        theta[1][0] = tt;


        // W
        double Wt = W[1][1];
        W[0][0] = 0;
        W[0][1] = 0;
        W[1][1] = theta[1][1] * sqrt(dt) + dAdv(x[1][1]) * sqrt(2 * eta * T) * z() * dt * sqrt(dt);
        //replace past values
        W[1][0] = Wt;

        /*
        cout << "normal: ";
        cout << theta[1] << endl;
        */
    }

    // Hilfsfunktion
    double A_func(double v, double q) {
        return - (eta) * v - dVdq(q);
    }

    double dAdv(double v) const {
        return (-eta);
    }

    static double z() {
        return 1/2 * (normal_dice() + 1/sqrt(3) * normal_dice());
    }

public:
    void run(int steps, const string name="Brownian Motion in Potential Data") {
        // runs brownian motion with stepsize dt for steps steps
        // initialize
        state_type x = vector<vector<double>>(2, vector<double>(2));
        // set initial values
        x[0][1] = x0;
        x[1][1] = v0;
        cout << x[0][0] << ", " << x[0][1] << endl;
        cout << x[1][0] << ", " << x[1][1] << endl;

        state_type theta = vector<vector<double>>(2, vector<double>(2));
        state_type A = vector<vector<double>>(2, vector<double>(2));
        state_type W = vector<vector<double>>(2, vector<double>(2));
        // set initial time to 0
        double t = 0;

        std::ofstream myfile;
        myfile.open (name + ".csv");

        // average value
        double mu = 0;
        for(int i = 0; i < steps; i++) {
            brownian_step(x, A, theta, W, t);
            stochastic_two_step(x, W, theta, A);

            // write to file
            myfile << t << ",";
            myfile << x[0][1] << ",";
            myfile << x[1][1] << ", \n";

            mu += x[0][1];
            t += dt;
            // print position
        }
        myfile.close();
        mu /= steps;
        cout << endl;
        cout << mu << endl;

    }
    two_step(double eta_val, double dt_val, double T_val, double w_val, double x0_val, double v0_val) {
        eta = eta_val;
        dt = dt_val;
        T = T_val;
        w = w_val;
        x0 = x0_val;
        v0 = v0_val;
    }
};

int main() {
    double eta = 1;
    double T = 1; // 10
    double dt = 0.00001;
    double w = 5;   // 1
    int steps = 1000000;
    cout << "Test";

    two_step suite = two_step(eta, dt, T, w, 2, 0);
    // number of runs

    int runs = 2;
    for(int i = 0; i < runs; i++) {
        string name = "Brownian Motion Two Step/" + to_string(i);
        suite.run(steps, name);
    }


    return 0;
}

