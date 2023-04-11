//
// Created by andi on 28.03.23.
//


#include <iostream>
#include <boost/numeric/odeint.hpp>
#include <cmath>
#include <random>
#include <fstream>
#include <filesystem>


// using namespaces!
using namespace std;
using namespace boost::numeric::odeint;
// new state is composed of every
typedef vector<double> entry_type;
typedef  vector<vector<entry_type>> state_type;

int pymod(int a, int b) {
    return ((b + (a % b)) % b);
}

class lattice_model {
// unser state_type ist einfach 2d Vektor, lass uns ublas matrix nehmen? Oder normalen vector?
    double eta;
    //static const auto dice = bind(normal_distribution<double>(0, 1), default_random_engine());
    // _Bind<normal_distribution(linear_congruential_engine<unsigned long, 16807, 0, 2147483647>)> dice;
    double T;
    double dt;
    static inline auto normal_dice = bind(normal_distribution<double>(0, 1), default_random_engine(1234));
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

    double dVidq(double qij, double qijp1, double qip1j, double qijm1, double qim1j) {
        return J * (sin(qij - qijp1) + sin(qij - qip1j) + sin(qij - qijm1) + sin(qij - qim1j));
    }

    double dVdq(double qij, double qijp1, double qip1j, double qijm1, double qim1j, double t) {
        return dVsdq(qij, t) + dVidq(qij, qijp1, qip1j, qijm1, qim1j);
    }


    void brownian_step(const state_type &x, state_type &dxdt, state_type &theta, double t){
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

public:
    void run(int steps, std::ofstream& file) {
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



        // average value
        double mu = 0;
        for(int k = 0; k < steps; k++) {
            brownian_step(x, dxdt, theta, t);
            stochastic_euler_method(x, dxdt, theta);

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

            t += dt;
        }
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


void search_grid(vector<double> eta_values, vector<double> T_values, vector<double> dt_values, vector<int> steps_values, vector<int> n_values,
                 vector<double> alpha_values, vector<double> beta_values, vector<double> J_values, vector<double> tau_values,
                 string storage_root) {
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
                                        double starting_t = -10;
                                        vector<vector<double>> x0 =
                                                vector<vector<double>>(n, vector<double>(n, 0));
                                        vector<vector<double>> v0 =
                                                vector<vector<double>>(n, vector<double>(n, 0));

                                        lattice_model suite =
                                                lattice_model(eta, dt, T,n , x0, v0,
                                                                            alpha, beta, J,
                                                                            tau, starting_t);
                                        // number of runs
                                        int runs = 1;
                                        for(int i = 0; i < runs; i++) {
                                            // ok we construct directory tree and if
                                            // it is not empty, we count the number of files inside and then just
                                            // number our runs
                                            string dir_name = create_directory(eta, T, dt, n, alpha, beta, J, tau,
                                                                               storage_root);
                                            // since dir_name is per construction
                                            // a directory we dont have to check if its a directory?
                                            string name = dir_name + "/" + to_string(count_files(dir_name));
                                            std::ofstream file;
                                            file.open (name + ".csv");
                                            suite.run(steps, file);
                                            write_parameters(file, eta, T, dt, n, alpha, beta, J, tau);
                                            file.close();
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
}



int main() {
    double eta = 5;
    double T = 0.01; // 10
    double dt = 0.002;
    int steps = 10001;
    int n = 20;
    double alpha = 1;   // 1
    double beta  = 5;
    double J = 10;
    double tau = 2;
    double starting_t = -10;



    /*
    // single search config
    vector<vector<double>> x0 = vector<vector<double>>(n, vector<double>(n, 0));
    vector<vector<double>> v0 = vector<vector<double>>(n, vector<double>(n, 0));

    lattice_model suite = lattice_model(eta, dt, T,n , x0, v0,
                                        alpha, beta, J, tau, starting_t);
    // number of runs
    int runs = 1;
    for(int i = 0; i < runs; i++) {
        // ok we construct directory tree and if it is not empty, we count the number of files inside and then just
        // number our runs
        string dir_name = create_directory(eta, T, dt, n, alpha, beta, J, tau);
        // since dir_name is per construction a directory we dont have to check if its a directory?
        string name = dir_name + "/" + to_string(count_files(dir_name));
        std::ofstream file;
        suite.run(steps, file, name);
        write_parameters(file, eta, T, dt, n, alpha, beta, J, tau, name);

    }
    */
    // important: choose the storage root outside of your git repo
    string storage_root = "../../../Generated content/";

    // TODO maybe replace dt with tmax
    vector<double> eta_values{5};
    vector<double> T_values{0.1, 1};
    vector<double> dt_values{0.00025};
    vector<int> steps_values{80001};
    vector<int> n_values{20};
    vector<double> alpha_values{1};
    vector<double> beta_values{5};
    vector<double> J_values{1, 5};
    vector<double> tau_values{1, 0.1};

    search_grid(eta_values, T_values, dt_values, steps_values, n_values, alpha_values, beta_values, J_values, tau_values,
                storage_root);

    return 0;
}

