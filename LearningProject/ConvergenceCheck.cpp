//
// Created by andi on 06.04.23.
//
#include "Helpfunctions and Classes.h"
#include "Systems.h"
#include "Solvers.h"

string create_name(solver* suite, double dt, int t_end) {
    string name = suite->get_name() + " dt=" +to_string(dt) + " t=" + to_string(t_end);
    return name;
};


void run_check(ofstream& file, int runs, double msd_harmonic_osc, double mu_harmonic_osc, int steps, solver* &suite) {
    // we still got to save this shit
    double error_mu = 0;
    double error_msd = 0;
    for(int i = 0; i < runs; i++) {
        state_type x = suite->run(steps);
        double mu = suite->calc_mu(x);
        double msd = suite->calc_msd(x, 0);
        file << "supposed value: mu = " << mu_harmonic_osc << "   actual value: mu = " << mu << endl;
        file << "supposed value: msd = " << msd_harmonic_osc << "   actual value: msd = " << msd << endl;
        error_mu += abs(mu_harmonic_osc - mu);
        error_msd += abs(msd_harmonic_osc - msd);
    }
    error_mu /= runs;
    error_msd /= runs;
    file << "Average deviation dmu = " << error_mu << endl;
    file << "Average deviation dmsd = " << error_msd << endl;

}


int main() {
    // we wont to check whether the stochastic euler method converges.
    // for brownian motion in harmonic oscillators, the
    // We want to search a grid again for different ending times, different dt and different solvers
    // We then want to save the results in probably different txt files
    double eta = 5;
    double T = 100; // 10
    double passed_time = 100;
    int n = 100;
    double alpha = 1;   // 1
    int runs = 1;
    // lets only make the here relevant parameters vectors
    // different ending times to see if we are already in the long time limit
    vector<int> ending_times = {10, 1, 100};
    // different number of steps to vary dt
    // search large phase space
    vector<int> steps_values = {10000, 100000, 1000000};
    // different methods of solving (here LM vs EM)
    // but for that we have to construct the suits first

    // first we need again a root directory where we want to save our results
    string root = "../../Generated content/Convergence Check Grid 2/";
    // create the directory if it doesnt exist alreade
    create_dir(root);


    // thoretical msd for brownian motion in a harmonic oscillator
    double msd_harmonic_osc = T / alpha;
    double mu_harmonic_osc = 0;

    // IC (don't really matter i think  but we have to have em)
    vector<vector<double>> x0 =
            vector<vector<double>>(n, vector<double>(n, 0));
    vector<vector<double>> v0 =
            vector<vector<double>>(n, vector<double>(n, 0));
    // harmonic system get  initialized once so we can keep that
    harmonic_system HarmonicSystem =
            harmonic_system(eta, T, n, x0, v0, alpha);

    // okay now we got to iterate
    int k = 1;
    int nr_configs = ending_times.size() * steps_values.size();
    cout << "Starting Check for " << nr_configs << " Configurations" << endl;
    for (int t_end : ending_times) {
        for(int steps : steps_values) {
            cout << "config " << k << "/" << nr_configs << endl;
            k++;
            double dt = static_cast<double>(t_end) / steps;
            cout << dt << endl;
            // weil ich nicht gut genug c++ schreiben kann mach ich die initalisierung für jeden solver einzeln
            // kurz die solver in nem container sammeln
            vector<solver*> suits;
            // i hope it gets initialized new as soon as the loop runs a second time?
            solver euler_suite = solver(dt, &HarmonicSystem);
            suits.push_back(&euler_suite);
            min_lm_solver lm_suite = min_lm_solver(dt, &HarmonicSystem);
            suits.push_back(&lm_suite);
            // We got to run here because otherways i don't know the steps value for the created suite
            // TODO: bzw. why is dt a class variable? I know it is comfortable but it would be more flexilbe to give it
            // as a parameter to the run method
            for(solver* suite: suits) {
                // create name
                cout << suite->get_name() << endl;
                string file_name = create_name(suite, dt, t_end);
                // whole name with root
                string name = root + file_name;
                // create ofstream
                ofstream file;
                // save everything in txt
                file.open(name + ".txt");
                // function to create a name for the txt file
                file << "Evaluation time t = " << dt * steps << endl;
                file << "Step size dt = " << dt << endl;
                file << "n = " << n << endl;
                file << "nr of runs: " << runs << endl;
                run_check(file, runs, msd_harmonic_osc, mu_harmonic_osc, steps, suite);
                file.close();
            }
        }
    }
    return 0;
}
