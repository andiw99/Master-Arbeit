//
// Created by andi on 18.04.23.
//

#include "ConvergenceCheckDiffusion.h"
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


void run_check(ofstream& file, ofstream& para_file, int runs, double msd_harmonic_osc, double mu_harmonic_osc, int steps, solver* &suite) {
    // we still got to save this shit
    double error_mu = 0;
    double error_msd = 0;
    double mu = 0;
    double msd = 0;
    for(int i = 0; i < runs; i++) {
        state_type x = suite->run(steps);
        // calculate the moments for this run
        double curr_mu = suite->calc_mu(x);
        double curr_msd = suite->calc_msd(x, 0);
        // add to the average value
        mu += curr_mu;
        msd += curr_msd;
        para_file << "supposed value: mu = " << mu_harmonic_osc << "   actual value: mu = " << curr_mu << endl;
        para_file << "supposed value: msd = " << msd_harmonic_osc << "   actual value: msd = " << curr_msd << endl;
        error_mu += abs(mu_harmonic_osc - mu);
        error_msd += abs(msd_harmonic_osc - msd);
    }
    // average errors and moments over the runs
    error_mu /= runs;
    error_msd /= runs;
    // use relative error
    error_mu /= mu_harmonic_osc;
    error_msd /= msd_harmonic_osc;

    mu /= runs;
    msd /= runs;
    // write to the files
    para_file << "Average deviation dmu = " << error_mu << endl;
    para_file << "Average deviation dmsd = " << error_msd << endl;

    file << mu << ", " << msd << ", " << error_mu << ", " << error_msd << ", ";

}


int main() {
    // we wont to check whether the stochastic euler method converges.
    // for brownian motion in harmonic oscillators, the
    // We want to search a grid again for different ending times, different dt and different solvers
    // We then want to save the results in probably different txt para_files
    double eta = 5;
    double T = 100; // 10
    double passed_time = 1;
    int n = 200;
    int runs = 1;
    // lets only make the here relevant parameters vectors
    // different ending times to see if we are already in the long time limit
    double t = 20;
    // different number of steps to vary dt
    // search large phase space
    vector<int> steps_values = {200, 400, 1000, 2000, 5000, 20000, 100000};
    // different methods of solving (here LM vs EM)
    // but for that we have to construct the suits first

    // first we need again a root directory where we want to save our results
    string root = "../../Generated content/Convergence Check Diffusion more time/";
    // create the directory if it doesnt exist alreade
    create_dir(root);


    // thoretical msd for brownian motion without potential 2 * D * t
    // D = k_B * T / eta
    double D = T / eta;
    double msd = 2 * D * t;
    double mu = 0;

    // IC (don't really matter i think  but we have to have em)
    vector<vector<double>> x0 =
            vector<vector<double>>(n, vector<double>(n, 0));
    vector<vector<double>> v0 =
            vector<vector<double>>(n, vector<double>(n, 0));
    // harmonic system get  initialized once so we can keep that
    brown_system brownian = brown_system(T, n, x0, v0, eta);

    // okay now we got to iterate
    int k = 1;
    int nr_configs = steps_values.size();
    cout << "Starting Check for " << nr_configs << " Configurations" << endl;
    // create ofstream for the results
    ofstream file;
    file.open(root + "results.csv");
    // write header into file
    // doch net
    for(int steps : steps_values) {
        cout << "config " << k << "/" << nr_configs << endl;
        k++;
        double dt = static_cast<double>(t) / steps;
        // dt into the file
        file << dt << ", ";
        cout << dt << endl;
        // weil ich nicht gut genug c++ schreiben kann mach ich die initalisierung für jeden solver einzeln
        // kurz die solver in nem container sammeln
        vector<solver*> suits;
        // i hope it gets initialized new as soon as the loop runs a second time?
        solver euler_suite = solver(dt, &brownian);
        suits.push_back(&euler_suite);
        min_lm_solver lm_suite = min_lm_solver(dt, &brownian);
        suits.push_back(&lm_suite);
        // We got to run here because otherways i don't know the steps value for the created suite
        // TODO: bzw. why is dt a class variable? I know it is comfortable but it would be more flexilbe to give it
        // as a parameter to the run method
        for(solver* suite: suits) {
            // create name
            cout << suite->get_name() << endl;
            string para_file_name = create_name(suite, dt, t);
            // whole name with root
            string name = root + para_file_name;
            // create ofstream for the parameters and stuff
            ofstream para_file;
            // save everything in txt
            para_file.open(name + ".txt");
            // function to create a name for the txt para_file
            para_file << "Evaluation time t = " << dt * steps << endl;
            para_file << "Step size dt = " << dt << endl;
            para_file << "n = " << n << endl;
            para_file << "nr of runs: " << runs << endl;
            run_check(file, para_file, runs, msd, mu, steps, suite);
            para_file.close();

        }
        // Zeilenumbruch
        file << "\n";
    }
    return 0;
}
