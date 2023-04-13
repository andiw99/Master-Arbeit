//
// Created by andi on 06.04.23.
//

#ifndef LEARNINGPROJECT_FUNCTIONS_AND_CLASSES_H
#define LEARNINGPROJECT_FUNCTIONS_AND_CLASSES_H

//
// Created by andi on 28.03.23.
//


#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <random>
#include <fstream>
#include <filesystem>
#include <utility>


// using namespaces!
using namespace std;
using namespace boost::numeric::odeint;
// new state is composed of every
typedef vector<double> entry_type;
typedef  vector<vector<entry_type>> state_type;

int pymod(int a, int b) {
    return ((b + (a % b)) % b);
}

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
    lattice_system() {
        n = 0;
    }
};


class cooling_bath: public lattice_system {
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
    // TODO worry about that later...
    static inline auto normal_dice =
            bind(normal_distribution<double>(0, 1), default_random_engine(1234));
    // parameters of the potential
    double alpha;
    double beta;
    // strenght of the interaction
    double J;
    // timescale of the quench
    double tau;
    // TODO phenomenological derivation of the behaviour of the system can't have less parameters

    double dVsdq(double q) {
        return alpha * (2 * q * q * q - beta * q);
    }

    double linear_T(double t) {
        // parametrisierung für die Temperatur
        // linearer Abfall
        current_T = max(T - t/tau, 1.0);
        return current_T;
    }

    double dVidq(double qij, double qijp1, double qip1j, double qijm1, double qim1j) const {
        return J * (sin(qij - qijp1) + sin(qij - qip1j) + sin(qij - qijm1) + sin(qij - qim1j));
    }

    double dVdq(double qij, double qijp1, double qip1j, double qijm1, double qim1j) {
        return dVsdq(qij) + dVidq(qij, qijp1, qip1j, qijm1, qim1j);
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

    cooling_bath( double eta_val, double T_val,  int n_val,
                             vector<vector<double>> x0_val,
                             vector<vector<double>> v0_val,
                             double alpha_val, double beta_val, double J_val, double tau_val)
            : lattice_system(n_val, T_val, x0_val, v0_val) {
        eta = eta_val;
        current_T = T_val;
        alpha = alpha_val;
        beta = beta_val;
        J = J_val;
        tau = tau_val;
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


class solver {
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
    void stochastic_euler_method(state_type &x, const state_type &dxdt, const state_type &theta) const {
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
                stochastic_euler_method(x, dxdt, theta);

                save_to_file(x, t, file);

                t += dt;
            }
        } else {
            for(int k = 0; k < steps; k++) {
                LatticeSystem->rhs(x, dxdt, theta, t);
                stochastic_euler_method(x, dxdt, theta);
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
                stochastic_euler_method(x, dxdt, theta);
                // only save if k is a multiple of save_interval
                if(k % save_interval == 0 || k == steps - 1) {
                    save_to_file(x, t, file);
                }

                t += dt;
            }
        } else {
            for(int k = 0; k < steps; k++) {
                LatticeSystem->rhs(x, dxdt, theta, t);
                stochastic_euler_method(x, dxdt, theta);
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

string trunc_double(double a, int precision=2) {
    stringstream stream;
    stream << std::fixed << std::setprecision(precision) << a;
    return stream.str();
}


string create_directory(double eta, double T, double dt, int n, double alpha, double beta, double J, double tau,
                        const string root) {
    string dir_name = root + "eta=" + trunc_double(eta)
                      + "/T=" + trunc_double(T) + "/dt=" +
                      trunc_double(dt, 4) + "/n="+ to_string(n) + "/alpha=" + trunc_double(alpha) + "/beta=" +
                      trunc_double(beta) + "/J=" + trunc_double(J) + "/tau=" + trunc_double(tau);
    // check wheter the directory already exists, if not create it
    if(!filesystem::is_directory(dir_name) || !filesystem::exists(dir_name)) {
        filesystem::create_directories(dir_name);
    }
    return dir_name;
}


void write_parameters(ofstream& file, double eta, double T, double dt, int n, double alpha, double beta, double J,
                      double tau) {
    // insert the parameters
    file << "eta," << eta << ", \n";
    file << "T," << T << ", \n";
    file << "dt," << dt << ", \n";
    file << "n," << n << ", \n";
    file << "alpha," << alpha << ", \n";
    file << "beta," << beta << ", \n";
    file << "J," << J << ", \n";
    file << "tau," << tau << ", \n";
}


int count_files(string dir_name) {
    int i = 0;
    for(const auto & entry : filesystem::directory_iterator(dir_name)) {
        ++i;
    }
    return i;
}

/**
 * lattice_system is the CLASS (not an object) to create instances of
 * Okay i think what i wanted to do is impossible
 * Better Verision of the search grid function
 * @param paras vector of vectors of parameters
 */
template<typename T>
void search_grid_v2(vector<vector<double>> paras) {

    for(vector<double> para_set : paras) {

    }
}

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
#endif //LEARNINGPROJECT_FUNCTIONS_AND_CLASSES_H
