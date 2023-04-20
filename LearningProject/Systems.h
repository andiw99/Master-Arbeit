//
// Created by andi on 14.04.23.
//

#ifndef LEARNINGPROJECT_SYSTEMS_H
#define LEARNINGPROJECT_SYSTEMS_H

#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <random>
#include <fstream>
#include <filesystem>
#include <utility>
#include "Helpfunctions and Classes.h"


// using namespaces!
using namespace std;
using namespace boost::numeric::odeint;
// new state is composed of every
typedef vector<double> entry_type;
typedef  vector<vector<entry_type>> state_type;


/**
 * Base Class for all Systems on a lattice
 */
class lattice_system {
protected:
    int n;      // lattice_system size, this is all taylored to lattices!
    double T;   // Temperature of the System
public:
    vector<vector<double>> x0;
    vector<vector<double>> v0;
    virtual void rhs(const state_type &x, state_type &dxdt, state_type &theta, double t) {
    };

    int get_n() const {
        return n;
    }

    virtual double get_temp() {
        return T;
    }

    lattice_system(int n_val, double T_val, vector<vector<double>> x0_val, vector<vector<double>> v0_val) {
        n = n_val;
        T = T_val;
        x0 = std::move(x0_val);
        v0 = std::move(v0_val);
    }
};

class bath: public lattice_system {
protected:
    // we have roughly the same parameters as the full interaction system
    double eta;
    //static const auto dice = bind(normal_distribution<double>(0, 1), default_random_engine());
    // _Bind<normal_distribution(linear_congruential_engine<unsigned long, 16807, 0, 2147483647>)> dice;
    // Temperature of the bath, will be cooled down to almost zero over the course of the phase Transition
    // The phase Transition happens roughly when the kick strength of the thermal fluctuations are of the order
    // of the height of the dip
    double current_T;              // current T
    // I would like to set the seed somehow but that doesnt seem to be possible so easily.
    // probably with a template but that seems overly complicated
    double tau;

    // TODO worry about that later...
    static inline auto normal_dice =
            bind(normal_distribution<double>(0, 1), default_random_engine(1234));

    double linear_T(double t) {
        // parametrisierung für die Temperatur
        // linearer Abfall
        current_T = max(T - t/tau, 1.0);
        return current_T;
    }

    virtual double dVdq(double qij, double qijp1, double qip1j, double qijm1, double qim1j) {
    }

public:
    void rhs(const state_type &x, state_type &dxdt, state_type &theta, double t) override{
        for(int i = 0; i < n; i++) {
            for(int j=0; j < n; j++) {
                dxdt[i][j][0] = x[i][j][1];
                // TODO Problem when calculating this are the border terms but here we can just use a modulo for now?
                // TODO Border terms with i - 1 have the problem to be (-1) sometimes, using modulo like in python
                // might be inefficient?
                dxdt[i][j][1] = (-eta) * x[i][j][1] - dVdq(x[i][j][0],
                                                           x[i][(j + 1) % n][0],
                                                           x[(i+1) % n][j][0],
                                                           x[i][pymod((j-1), n)][0],
                                                           x[pymod(i-1, n)][j][0]);
                theta[i][j][0] = 0;
                theta[i][j][1] = sqrt(2 * eta * linear_T(t)) * normal_dice();
            }
        }
    }

    double get_temp() override {
        return current_T;
    }

    bath(double eta_val, double T_val,  int n_val,
         vector<vector<double>> x0_val,
         vector<vector<double>> v0_val,
         double tau_val) : lattice_system(n_val, T_val, x0_val, v0_val) {
        eta = eta_val;
        current_T = T_val;
        tau = tau_val;
    }
};

class cooling_bath: public bath {

    // parameters of the potential
    double alpha;
    double beta;
    // strenght of the interaction
    double J;
    // timescale of the quench

    // TODO phenomenological derivation of the behaviour of the system can't have less parameters

    double dVsdq(double q) {
        return alpha * (2 * q * q * q - beta * q);
    }


    double dVidq(double qij, double qijp1, double qip1j, double qijm1, double qim1j) const {
        return J * (sin(qij - qijp1) + sin(qij - qip1j) + sin(qij - qijm1) + sin(qij - qim1j));
    }

    double dVdq(double qij, double qijp1, double qip1j, double qijm1, double qim1j) override{
        return dVsdq(qij) + dVidq(qij, qijp1, qip1j, qijm1, qim1j);
    }

public:
    cooling_bath( double eta_val, double T_val,  int n_val,
                  vector<vector<double>> x0_val,
                  vector<vector<double>> v0_val,
                  double alpha_val, double beta_val, double J_val, double tau_val)
            : bath(eta_val, T_val, n_val, x0_val, v0_val, tau_val) {
        alpha = alpha_val;
        beta = beta_val;
        J = J_val;
    }
};

class silicon_bath : public bath {
    double J_para;
    double J_ortho;
    double alpha;
    double beta;

    double dVsdq(double q) {
        return alpha * (2 * q * q * q - beta * q);
    }

    double dVidq(double qij, double qijp1, double qip1j, double qijm1, double qim1j) const {
        return - J_ortho * (sin(qij - qijp1) + sin(qij - qijm1)) -
        J_para * (sin(qij - qip1j)  + sin(qij - qim1j));
    }

    double dVdq(double qij, double qijp1, double qip1j, double qijm1, double qim1j) {
        return dVsdq(qij) + dVidq(qij, qijp1, qip1j, qijm1, qim1j);
    }

public:
    silicon_bath( double eta_val, double T_val,  int n_val,
                  vector<vector<double>> x0_val,
                  vector<vector<double>> v0_val,
                  double alpha_val, double beta_val, double J_para_val, double J_ortho_val, double tau_val)
            : bath(eta_val, T_val, n_val, x0_val, v0_val, tau_val) {
        alpha = alpha_val;
        beta = beta_val;
        J_para = J_para_val;
        J_ortho = J_ortho_val;
    }
};


class full_interaction_system: public lattice_system {
    double eta;
    //static const auto dice = bind(normal_distribution<double>(0, 1), default_random_engine());
    // _Bind<normal_distribution(linear_congruential_engine<unsigned long, 16807, 0, 2147483647>)> dice;
    double T;
    static inline auto normal_dice = bind(normal_distribution<double>(0, 1), default_random_engine(1234));
    double alpha;
    double beta;
    double J;
    double tau;
    int n;      // size of the lattice

    double dVsdq(double q, double t) {
        return alpha * (2 * q * q * q - beta * eps(t) * q);
    }

    double eps(double t) const {
        // TODO normally with time scale, but is this important here?
        return min(t/tau, 1.0);
    }

    double dVidq(double qij, double qijp1, double qip1j, double qijm1, double qim1j) const {
        return J * (sin(qij - qijp1) + sin(qij - qip1j) + sin(qij - qijm1) + sin(qij - qim1j));
    }

    double dVdq(double qij, double qijp1, double qip1j, double qijm1, double qim1j, double t) {
        return dVsdq(qij, t) + dVidq(qij, qijp1, qip1j, qijm1, qim1j);
    }

public:
    void rhs(const state_type &x, state_type &dxdt, state_type &theta, double t) override{
        for(int i = 0; i < n; i++) {
            for(int j=0; j < n; j++) {
                dxdt[i][j][0] = x[i][j][1];
                // TODO Problem when calculating this are the border terms but here we can just use a modulo for now?
                // TODO Border terms with i - 1 have the problem to be (-1) sometimes, using modulo like in python
                // might be inefficient?
                dxdt[i][j][1] = (-eta) * x[i][j][1] - dVdq(x[i][j][0],
                                                           x[i][(j + 1) % n][0],
                                                           x[(i+1) % n][j][0],
                                                           x[i][pymod((j-1), n)][0],
                                                           x[pymod(i-1, n)][j][0],
                                                           t);
                theta[i][j][0] = 0;
                theta[i][j][1] = sqrt(2 * eta * T) * normal_dice();
            }
        }
    }
    full_interaction_system( double eta_val, double T_val,  int n_val,
                             vector<vector<double>> x0_val,
                             vector<vector<double>> v0_val,
                             double alpha_val, double beta_val, double J_val, double tau_val)
            : lattice_system(n_val, T_val, x0_val, v0_val) {
        eta = eta_val;
        T = T_val;
        n = n_val;
        alpha = alpha_val;
        beta = beta_val;
        J = J_val;
        tau = tau_val;
    }
};

class harmonic_system: public lattice_system {
    double eta;
    //static const auto dice = bind(normal_distribution<double>(0, 1), default_random_engine());
    // _Bind<normal_distribution(linear_congruential_engine<unsigned long, 16807, 0, 2147483647>)> dice;
    double T;
    static inline auto normal_dice = bind(normal_distribution<double>(0, 1), default_random_engine(1234));
    double alpha;
    int n;      // size of the lattice

    double dVdq(double q, double t) {
        return alpha * q;
    }


public:
    void rhs(const state_type &x, state_type &dxdt, state_type &theta, double t) override{
        for(int i = 0; i < n; i++) {
            for(int j=0; j < n; j++) {
                dxdt[i][j][0] = x[i][j][1];
                dxdt[i][j][1] = (-eta) * x[i][j][1] - dVdq(x[i][j][0], t);
                theta[i][j][0] = 0;
                theta[i][j][1] = sqrt(2 * eta * T) * normal_dice();
            }
        }
    }
    harmonic_system( double eta_val, double T_val,  int n_val,
                     vector<vector<double>> x0_val,
                     vector<vector<double>> v0_val,
                     double alpha_val)
            : lattice_system(n_val, T_val, x0_val, v0_val) {
        eta = eta_val;
        T = T_val;
        n = n_val;
        alpha = alpha_val;
    }
};

class brown_system: public lattice_system {
    double eta;
    static inline auto normal_dice = bind(normal_distribution<double>(0, 1), default_random_engine(1234));

public:
    void rhs(const state_type &x, state_type &dxdt, state_type &theta, double t) override {
        for(int i = 0; i < n; i++) {
            for(int j=0; j < n; j++) {
                dxdt[i][j][0] = x[i][j][1];
                dxdt[i][j][1] = (-eta) * x[i][j][1];
                theta[i][j][0] = 0;
                theta[i][j][1] = sqrt(2 * eta * T) * normal_dice();
            }
        }
    }

    brown_system(double T_val,  int n_val,
                 vector<vector<double>> x0_val,
                 vector<vector<double>> v0_val,
                 double eta_val) : lattice_system(n_val, T_val, x0_val, v0_val) {
        eta = eta_val;
    }
};

#endif //LEARNINGPROJECT_SYSTEMS_H
