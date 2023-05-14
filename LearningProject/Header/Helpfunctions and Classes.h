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


string print_statetype(const state_type x) {
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

template <class T>
void print_vector(vector<T> vec) {
    cout << vec.size() << endl;
    for(int i = 0; i < vec.size(); i++) {
        if(i % 20 == 0) {
            cout << "\n";
        }
        cout << vec[i] << ", ";
    }
}


ifstream safe_read(string readpath) {
    cout << "reading: " << readpath << endl;
    ifstream file(readpath);
    // check if file is opened
    if(file.is_open()) {
        cout << "File successfully opened" << endl;
    } else {
        cout << "Failed to open file" << endl;
        // abort if file is not opened
        exit(0);
    }
    return file;
}


string trunc_double(double a, int precision=2) {
    stringstream stream;
    stream << std::fixed << std::setprecision(precision) << a;
    return stream.str();
}


void create_dir(const string dir_name) {
    // check wheter the directory already exists, if not create it
    if(!filesystem::is_directory(dir_name) || !filesystem::exists(dir_name)) {
        filesystem::create_directories(dir_name);
    }
}

string create_tree_name(double eta, double T, double dt, int n, double alpha, double beta, double J, double tau,
                        const string root) {
    string dir_name = root + "eta=" + trunc_double(eta)
                      + "/T=" + trunc_double(T) + "/dt=" +
                      trunc_double(dt, 4) + "/n="+ to_string(n) + "/alpha=" + trunc_double(alpha) + "/beta=" +
                      trunc_double(beta) + "/J=" + trunc_double(J) + "/tau=" + trunc_double(tau);
    create_dir(dir_name);
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


double ind_value(vector<double> paras, int ind) {
    return paras[ind % paras.size()];
}


int ind_value(vector<int> paras, int ind) {
    return paras[ind % paras.size()];
}

template <class T>
vector<T>& operator+ (vector<T> &a, vector<T> &b) {
    if (a.size() != b.size()) {
        throw std::invalid_argument("Vectors must be of equal length.");
    }
    for(int i = 0; i < a.size(); i++) {
        a[i] += b[i];
    }
    return a;
}

#endif //LEARNINGPROJECT_HELPFUNCTIONS_AND_CLASSES_H
