//
// Created by andi on 14.04.23.
//

#include "../Odeint/OdeintTest.h"
//
// Created by weitze73 on 13.03.23.
//

#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <cmath>

// using namespaces!
using namespace std;
using namespace boost::numeric::odeint;


// two-dimensional lattice, use ublas matrix with odeint
// double precision for now
// dont we need both the position and the momentum? i would think that the matrix should be filled with vectors
// or arrays? whatever has better performance

typedef vector<double> entry_type;
typedef boost::numeric::ublas::matrix<entry_type> state_type;

int py_mod(int a, int b) {
    // Implementation einer Funktion die sich verhält wie die python modulo funktion
    return (b + (a % b)) % b;
}

class kzm_silicon_2d {
    // I guess I set the "mass" in my ODE to be zero
    // so I have:
    // the damping eta
    // parameter alpha of Potential (steepness)
    // parameter beta of Potential
    // effective temperature? or maybe critical temp
    // Strength of Interaction in parallel direction J_para
    // Strength of Interaction in orthogonal direction J_ortho
    // Timescale of the Quench, tau_q
    double eta;
    double alpha;
    double beta;
    double eps;
    double J_para;
    double J_ortho;
    double tau_q;

    // Constructor, I use my classic syntax that i know from Java
    kzm_silicon_2d(double eta_val, double alpha_val, double beta_val, double J_para_val,
                   double J_ortho_val, double tau_q_val) {
        eta = eta_val;
        alpha = alpha_val;
        beta = beta_val;
        eps = -10.0;
        J_para = J_para_val;
        J_ortho = J_ortho_val;
        tau_q = tau_q_val;
    }

    // Implementation of custom potential

    double dV_dq(double q_ij, double q_ijm1, double q_im1j, double q_ijp1, double q_ip1j) {
        // "sombrero part"
        // splitting better for reading, but maybe significantly slower?
        // using std::pow faster or slower?
        double dV_dq_sombrero = alpha * (4 * pow(q_ij, 3) - 2 * beta * eps * q_ij);
        double dV_dq_interaction =  - J_para * (sin(q_ij - q_ijm1) + sin(q_ij - q_ijp1))
                                    - J_ortho * (sin(q_ij - q_im1j) + sin(q_ij - q_ip1j));
        return dV_dq_sombrero + dV_dq_interaction;
    }

    // Implementation of thermal fluctuations
    // TODO: Placeholder Implementation, real Implementation should draw the values randomly (from which dist?)
    /*
    double f(size_t i, size_t j, double t) {
        return 0.1;
    }
     */

    // Overload function f for integers
    double f(int i, int j, double t) {
        return 0.1;
    }

    // Implementation of ODE rhs
    // add const before {} ? What does it do?
    void operator() (const state_type &x, state_type &dxdt, double t) {
        // we only work with quadratic lattices, so only one size. Isnt that performance consuming?
        // we probably could make the lattice size a parameter of the class
        // for PBC we make a type conversion here
        int n = (int)x.size1();
        // parametrize eps? Linear dependency on t like near the phase transition
        // then we actually dont need to initialize eps at all. Should we maybe count the time in in multiples of
        // tau_q to save one calc everytime? -> probably not worth it...
        eps = t / tau_q;
        // think of periodic boundary conditions
        // nested loop i guess at least for the inner lattice sites
        for (int i=1; i < n - 1; i++) {
            for(int j=1; j < n - 1; j++) {
                step(i, j, x, dxdt, t);
            }
        }
        // PBC...
        for(int i = 0; i < n; i++) {
            boundary_step(i, x, dxdt, t, n);
        }
    }
    // one step of the iteration
    void step(int i, int j, const state_type &x, state_type &dxdt, double t) {
        // my ODE is of second order, so i need both equations, dont i? Why doesnt the Example on boost.org
        // for 2D lattices use it since it wants to deal with coupled oscillators?
        // ODE for angle
        dxdt(i, j)[0] = x(i, j)[1];
        // ODE for angle velocity
        dxdt(i, j)[1] = f(i, j, t) - eta * x(i, j)[0]
                        - dV_dq(x(i, j)[0], x(i, j - 1)[0],
                                x(i - 1, j)[0], x(i, j + 1)[0], x(i+1, j)[0]);
    }

    // special function for the boundaries
    void boundary_step(int i, const state_type &x, state_type &dxdt, double t, int n) {
        // first col
        dxdt(i, 0)[0] = x(i, 0)[1];
        dxdt(i, 0)[1] = f(i, 0, t) - eta * x(i, 0)[0]
                        - dV_dq(x(i, 0)[0], x(i, n - 1)[0],
                                x(py_mod(i, n), 0)[0], x(i, 0 + 1)[0], x((i+1) % n, 0)[0]);
        // last col
        dxdt(i, n-1)[0] = x(i, n-1)[1];
        dxdt(i, n-1)[1] = f(i, n-1, t) - eta * x(i, n-1)[0]
                          - dV_dq(x(i, n-1)[0], x(i, n - 2)[0],
                                  x(py_mod(i-1, n), n-1)[0], x(i, 0)[0], x((i+1) % n, n-1)[0]);
        // first row meaning i=0
        dxdt(0, i)[0] = x(0, i)[1];
        dxdt(0, i)[1] = f(0, i, t) - eta * x(0, i)[0]
                        - dV_dq(x(0, i)[0], x(0, py_mod(i-1, n))[0],
                                x(n - 1, i)[0], x(0, (i + 1) % n)[0], x(1, i)[0]);
        // last row (i = n-1)
        dxdt(n-1, i)[0] = x(n-1, i)[1];
        // ODE for angle velocn-1ty
        dxdt(n-1, i)[1] = f(n-1, i, t) - eta * x(n-1, i)[0]
                          - dV_dq(x(n-1, i)[0], x(n-1, py_mod(i-1, n))[0],
                                  x(n-2, i)[0], x(n-1, (i + 1) % n)[0], x(0, i)[0]);
    }
};

/*
int main() {
    return 0;
}
 */