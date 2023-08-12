//
// Created by weitze73 on 13.03.23.
//

#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <boost/random.hpp>
#include <fstream>
#include "../Header/Helpfunctions and Classes.h"
#include "../Header/Systems.h"
#include "../Header/Solvers.h"

// using namespaces!
using namespace std;
using namespace boost::numeric::odeint;


// two-dimensional lattice, use ublas matrix with odeint
// double precision for now
// dont we need both the position and the momentum? i would think that the matrix should be filled with vectors
// or arrays? whatever has better performance

// different typedefs
const static int N = 50;
typedef boost::array<double, 2> entry_type_odeint;
typedef  boost::array<boost::array<entry_type_odeint, N>, N> state_type_odeint;

string print_statetype(const state_type_odeint x) {
    string str = "";
    int n = x.size();
    str += "x: \n";
    str += "[";
    for (int i = 0; i < n; i++) {
        str += "[";
        for (int j = 0; j < n; j++) {
            if(j == n-1) {
                str += to_string(x[i][j][0]);
            } else {
                str += to_string(x[i][j][0]) + ", ";
            }
        }
        if(i < n-1) {
            str += "]\n ";
        }
    }
    str += "]]\n";
    str += "p: \n";
    str += "[";
    for (int i = 0; i < n; i++) {
        str += "[";
        for (int j = 0; j < n; j++) {
            if(j == n-1) {
                str += to_string(x[i][j][1]);
            } else {
                str += to_string(x[i][j][1]) + ", ";
            }
        }
        if(i < n-1) {
            str += "]\n ";
        }
    }
    str += "]]";
    return str;
}

class harmonic_system_odeint_det {
    double eta;
    double alpha;
    int n;      // size of the lattice
public:
    // Constructor, I use my classic syntax that i know from Java
    harmonic_system_odeint_det(double eta_val, double alpha_val) {
        eta = eta_val;
        alpha = alpha_val;

    }

    // Implementation of custom potential

    double dVdq(double q) {
        return alpha * q;
    }
    // Implementation of ODE rhs
    // add const before {} ? What does it do?
    void operator() (const state_type_odeint &x, state_type_odeint &dxdt) {

        n = x.size();
        for(int i = 0; i < n; i++) {
            for(int j = 0; j < n; j++) {
                dxdt[i][j][0] = x[i][j][1];
                dxdt[i][j][1] = (-eta) * x[i][j][1] - dVdq(x[i][j][0]);
            }
        }
    }
};

class harmonic_system_odeint_stoch {
    double T;
    double eta;
    boost::mt19937 &m_rng;
    boost::normal_distribution<> m_dist;
public:
    harmonic_system_odeint_stoch(double eta_val, double T_val, boost::mt19937 &rng , double sigma)
    : m_rng( rng ) , m_dist( 0.0 , sigma ) {
        eta = eta_val;
        T = T_val;
    };

    void operator() (const state_type_odeint &x, state_type_odeint &dxdt) {
        int n =x.size();

        for(int i = 0; i < n; i++) {
            for(int j=0; j < n; j++) {
                dxdt[i][j][0] = 0;
                dxdt[i][j][1] = sqrt(2 * eta * T) * m_dist(m_rng);

            }
        }

    }
};

template <int N>
class stochastic_euler
{
public:
    typedef boost::array<double, 2> entry_type;
    typedef boost::array<boost::array<entry_type, N>, N> state_type_odeint;
    typedef double time_type;
    typedef boost::numeric::odeint::stepper_tag stepper_category;


    template< class System >
    void do_step( System system , state_type_odeint &x , time_type t , time_type dt ) const
    {
        state_type_odeint det , stoch ;
        system.first( x , det );
        system.second( x , stoch );

        int n = x.size();

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                x[i][j][0] = det[i][j][0] * dt + stoch[i][j][0] * sqrt(dt);
                x[i][j][1] = det[i][j][1] * dt + stoch[i][j][1] * sqrt(dt);
            }
        }
    }
};

class streaming_observer
{
public:
    ofstream& file;
    size_t m_write_every;
    size_t m_count;
    int steps;
    state_type_odeint state;
    int test = 0;

    // now withougt fancy Temperature writing cause i am tired and
    streaming_observer(ofstream& out, int steps, int write_every = 100)
            : file(out), m_write_every( write_every ), m_count( 0 ), steps(steps) {}
    template< class State >
    void operator()( const State &x , double t )
    {
        if( ( m_count % m_write_every ) == 0 ) {
            // can we open and close the file again and again and not delete iti?
            int n = x.size();
            // write to file
            // everything?
            file << "t, " << t << ",";
            for(int i = 0; i < n; i++) {
                for(int j=0; j < n; j++) {
                    file << x[i][j][0] << ",";

                }
            }
            for(int i = 0; i < n; i++) {
                for(int j=0; j < n; j++) {
                    file << x[i][j][1] << ",";
                }
            }
            // Zeilenumbruch
            file << "\n";
        }
        ++m_count;
        // save last value
        if (m_count == steps) {
            test = 12;
            cout << m_count << endl;
            int n = x.size();
            for(int i = 0; i < n; i++) {
                for(int j=0; j < n; j++) {
                    state[i][j][0] = x[i][j][0];
                    state[i][j][1] = x[i][j][1];

                }
            }
            // Zeilenum;
        }
    }

    state_type_odeint get_state() {
        return state;
    }
};

int main() {
    boost::mt19937 rng;
    // N is now defined globally
    // int n = 20;
    double starting_t = 0.0;
    double end_t = 100.0;
    int steps = 10000;
    double dt = (end_t - starting_t) / steps;
    double eta = 5;
    double alpha = 1;
    double T = 100;
    string odeint_file = "../../Generated content/Odeint-Test/odeint.csv";
    string my_file = "../../Generated content/Odeint-Test/my.csv";

    std::ofstream file;
    file.open(odeint_file);
    state_type_odeint x0;
    state_type_odeint state;
    streaming_observer observer = streaming_observer(file, steps);
    integrate_const( stochastic_euler<N>() ,
            make_pair(harmonic_system_odeint_det(eta, alpha),
                      harmonic_system_odeint_stoch(eta, T, rng , 1.0 ) ),
                    x0 , starting_t , end_t , dt , observer);

    cout << print_statetype(observer.get_state()) << endl;
    cout << observer.test << endl;


    // How do we compare the two algorithms?
    vector<vector<double>> x02 =
            vector<vector<double>>(N, vector<double>(N, 0));
    vector<vector<double>> v0 =
            vector<vector<double>>(N, vector<double>(N, 0));
    harmonic_system sys = harmonic_system(eta, T, N, x02, v0, alpha);
}