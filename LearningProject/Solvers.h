//
// Created by andi on 14.04.23.
//

#ifndef LEARNINGPROJECT_SOLVERS_H
#define LEARNINGPROJECT_SOLVERS_H

#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <random>
#include <fstream>
#include <filesystem>
#include <utility>
#include "Systems.h"


// using namespaces!
using namespace std;
using namespace boost::numeric::odeint;
// new state is composed of every
typedef vector<double> entry_type;
typedef  vector<vector<entry_type>> state_type;


/**
 * Base class for my future solvers
 * Base class uses Euler-Mayurama method to update
 */
class solver {
protected:
// unser state_type ist einfach 2d Vektor, lass uns ublas matrix nehmen? Oder normalen vector?

    double dt;

    double tmin;       // starting time
    // TODO weil ich nicht weiß wie man das hier anständig macht muss das System beim solver init mit & übergeben werden
    lattice_system* LatticeSystem;
    int n;
    /*
    normal_distribution<double> normal_dist;
    default_random_engine generator;
    */
    void update(state_type &x, const state_type &dxdt, const state_type &theta) const {
        // huh we dont have vector operations, maybe implement them
        // We use the state_type here, so we need to cycle through every lattice
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                x[i][j][0] = x[i][j][0] + dxdt[i][j][0] * dt + theta[i][j][0] * sqrt(dt);
                x[i][j][1] = x[i][j][1] + dxdt[i][j][1] * dt + theta[i][j][1] * sqrt(dt);
            }
        }
    }

    void save_to_file(const state_type &x, double t, ofstream& file) {
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

        // Temperatur saven
        file << LatticeSystem->get_temp();

        // Zeilenumbruch
        file << "\n";
    }

    void set_IC(state_type& x) {
        for(int i = 0; i < n; i++) {
            for(int j=0; j < n; j++) {
                x[i][j][0] = LatticeSystem->x0[i][j];
                x[i][j][1] = LatticeSystem->v0[i][j];

            }
        }
    }

    void iterate_full(state_type& x, state_type& dxdt, state_type& theta, int steps, ofstream& file, bool save) {
        // set ininital time to starting time
        double t = tmin;

        if (save) {
            for(int k = 0; k < steps; k++) {
                LatticeSystem->rhs(x, dxdt, theta, t);
                update(x, dxdt, theta);

                save_to_file(x, t, file);

                t += dt;
            }
        } else {
            for(int k = 0; k < steps; k++) {
                LatticeSystem->rhs(x, dxdt, theta, t);
                update(x, dxdt, theta);
                t += dt;
            }
        }
    }

    void iterate_interval(state_type& x, state_type& dxdt, state_type& theta, int steps, ofstream& file, bool save,
                          int save_interval) {
        // set ininital time to starting time
        double t = tmin;

        if (save) {
            for(int k = 0; k < steps; k++) {
                LatticeSystem->rhs(x, dxdt, theta, t);
                update(x, dxdt, theta);
                // only save if k is a multiple of save_interval
                if(k % save_interval == 0 || k == steps - 1) {
                    save_to_file(x, t, file);
                }

                t += dt;
            }
        } else {
            for(int k = 0; k < steps; k++) {
                LatticeSystem->rhs(x, dxdt, theta, t);
                update(x, dxdt, theta);
                t += dt;
            }
        }
    }

public:
    /**
     *  runs the simulation
     * @param steps how many steps of the dgl are calculated
     * @param file the file
     * @return the first moments of the simulation, (averaged over every lattice site?)
     */
    state_type run_full_data(int steps, std::ofstream& file, bool save = true) {
        // runs brownian motion with stepsize dt for steps steps
        // initialize
        state_type x = state_type (n, vector<entry_type>(n, entry_type(2, 0)));
        // set initial values
        // naive setting because i don't know c++ to much atm and it only runs once when initializing
        set_IC(x);


        state_type dxdt = state_type (n, vector<entry_type>(n, entry_type(2, 0)));
        state_type theta = state_type (n, vector<entry_type>(n, entry_type(2, 0)));

        iterate_full(x, dxdt, theta, steps, file, save);

        return x;
    };
    /**
     *
     * @param steps
     * @param file
     * @param save
     * @param save_interval every "save_interval"-th value is saved in a file
     * @return
     */
    state_type run(int steps, std::ofstream& file, bool save = true, int save_interval=100) {
        // runs brownian motion with stepsize dt for steps steps
        // initialize
        state_type x = state_type (n, vector<entry_type>(n, entry_type(2, 0)));
        // set initial values
        // naive setting because i don't know c++ to much atm and it only runs once when initializing
        set_IC(x);


        state_type dxdt = state_type (n, vector<entry_type>(n, entry_type(2, 0)));
        state_type theta = state_type (n, vector<entry_type>(n, entry_type(2, 0)));

        iterate_interval(x, dxdt, theta, steps, file, save, save_interval);

        return x;
    };

    double calc_mu(const state_type& x) {
        double mu = 0;
        for(int i = 0; i < n; i++) {
            for(int j=0; j < n; j++) {
                mu += x[i][j][0];
            }
        }
        mu /= (n * n);
        return mu;
    };

    double calc_msd(const state_type& x, const double exp_value) {
        double msd = 0;
        for(int i = 0; i < n; i++) {
            for(int j=0; j < n; j++) {
                msd += pow((x[i][j][0] - exp_value), 2);
                // cout << pow((x[i][j][0] - exp_value), 2) << endl;
            }
        }
        msd = 1.0/(n * n) * msd;
        return msd;
    };

    double calc_msd(const state_type& x) {
        double mu = calc_mu(x);
        return calc_msd(x, mu);
    };

    vector<vector<double>> run_and_return {
            // TODO returns all values?
    };

    solver(double dt_val,
           class lattice_system *specific_system, double t_val = -5) {
        LatticeSystem = specific_system;
        dt = dt_val;
        tmin = t_val;
        n = LatticeSystem->get_n();
    }
};

/**
 * lm solver uses the Leimkuhler mathhews method to solve the SDE
 * LM-Method averages the generated random numbers over one step and achieves better convergence through that
 */
class lm_solver: public solver {
private:
    void update(state_type &x, const state_type &dxdt, const state_type &theta) {
        // huh we dont have vector operations, maybe implement them
        // We use the state_type here, so we need to cycle through every lattice
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                x[i][j][0] = x[i][j][0] + dxdt[i][j][0] * dt + 0.5 * (theta[i][j][0] + pre_theta[i][j][0]) * sqrt(dt);
                x[i][j][1] = x[i][j][1] + dxdt[i][j][1] * dt + 0.5 * (theta[i][j][1] + pre_theta[i][j][1]) * sqrt(dt);
                // set the previous values
                pre_theta[i][j][0] = theta[i][j][0];
                pre_theta[i][j][1] = theta[i][j][1];
            }
        }
    }


public:
    state_type pre_theta;
    lm_solver(double dt_val,
              class lattice_system *specific_system, double t_val = -5) :
            solver(dt_val, specific_system, t_val) {
        // IV of pre_theta zero
        pre_theta = vector<vector<entry_type>>(n, vector<vector<double>>(n, vector<double>(2, 0)));
    };
};

/**
 * Uses the fact that the DGL for x has no stochastic term
 */
class min_lm_solver: lm_solver {
    void update(state_type &x, const state_type &dxdt, const state_type &theta) {
        // huh we dont have vector operations, maybe implement them
        // We use the state_type here, so we need to cycle through every lattice
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                x[i][j][0] = x[i][j][0] + dxdt[i][j][0] ;
                x[i][j][1] = x[i][j][1] + dxdt[i][j][1] * dt + 0.5 * (theta[i][j][1] + pre_theta[i][j][1]) * sqrt(dt);
                // set the previous values
                pre_theta[i][j][1] = theta[i][j][1];
            }
        }
    }
public:
    min_lm_solver(double dt_val,
                  class lattice_system *specific_system, double t_val = -5) :
            lm_solver(dt_val, specific_system, t_val) {};
};


template <typename sys>
void init_and_run(double eta, double T, double dt, int steps, int n, double alpha, double beta, double J, double tau,
                  string storage_root, double starting_t) {
    vector<vector<double>> x0 =
            vector<vector<double>>(n, vector<double>(n, 0));
    vector<vector<double>> v0 =
            vector<vector<double>>(n, vector<double>(n, 0));

    sys specific_system =
            sys(eta, T, n, x0, v0, alpha, beta, J, tau);

    solver suite =
            solver(dt, &specific_system, starting_t);
    // number of runs
    int runs = 1;
    for (int i = 0; i < runs; i++) {
        // ok we construct directory tree and if
        // it is not empty, we count the number of files inside and then just
        // number our runs
        string dir_name = create_directory(eta, T, dt, n, alpha, beta, J, tau,
                                           storage_root);
        // since dir_name is per construction
        // a directory we dont have to check if its a directory?
        string name = dir_name + "/" + to_string(count_files(dir_name));
        std::ofstream file;
        std::ofstream para_file;
        file.open(name + ".csv");
        para_file.open(name + ".txt");
        suite.run(steps, file);
        write_parameters(para_file, eta, T, dt, n, alpha, beta, J, tau);
        file.close();
        para_file.close();

    }
}

template<typename sys>
void search_grid(vector<double> eta_values, vector<double> T_values, vector<double> dt_values, vector<int> steps_values, vector<int> n_values,
                 vector<double> alpha_values, vector<double> beta_values, vector<double> J_values, vector<double> tau_values,
                 string storage_root, double starting_t = 0) {
    int k = 0;
    int nr_configs = eta_values.size() * T_values.size() * dt_values.size() * steps_values.size() * n_values.size() *
                     alpha_values.size() * beta_values.size() * J_values.size() * tau_values.size();
    for(double eta : eta_values) {
        for(double T : T_values) {
            for(double dt : dt_values) {
                for(int steps : steps_values) {
                    for(int n : n_values) {
                        for(double alpha : alpha_values) {
                            for(double beta : beta_values) {
                                for(double J : J_values) {
                                    for(double tau : tau_values) {
                                        init_and_run<sys>(eta, T, dt, steps, n, alpha, beta, J, tau, storage_root,
                                                          starting_t);
                                        k++;
                                        cout << "run " << k << "/" << nr_configs << endl;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

#endif //LEARNINGPROJECT_SOLVERS_H
