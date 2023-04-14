//
// Created by andi on 06.04.23.
//

#ifndef LEARNINGPROJECT_HELPFUNCTIONS_AND_CLASSES_H
#define LEARNINGPROJECT_HELPFUNCTIONS_AND_CLASSES_H

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


#endif //LEARNINGPROJECT_HELPFUNCTIONS_AND_CLASSES_H
