//
// Created by andi on 30.05.23.
//

#ifndef CUDAPROJECT_PARAMETERS_CUH
#define CUDAPROJECT_PARAMETERS_CUH

#include <map>

map<string, double> single_calc_standard = {
        {"steps", 100000.0},
        {"dt", 0.001},
        {"T", 30},
        {"J", 10},
        {"alpha", 5},
        {"beta", 10},
        {"tau", 10},
        {"eta", 1.2},
        {"N", 2},
        {"nr_save_values", 250},
        {"x0", 0},
        {"p0", 0}
};

string tempscan_root = "../../Generated content/ROCm/";

map<string, double> temp_scan_standard = {
        {"steps", 100000.0},
        {"end_t", 25.0},
        {"dt_max", 0.001},
        {"dt_start", 0.0001},
        {"t_relax", 30.0},
        {"J", 10.0},
        {"alpha", 5.0},
        {"beta", 20.0},
        {"tau", 10.0},
        {"eta", 0.8},
        {"N", 2.0},
        {"nr_save_values", 32.0},
        {"x0", 4.0},
        {"p0", 4.0},
        {"repeat_nr", 1.0},
        {"min_temp",100.0},
        {"max_temp", 100.0},
        {"nr_temps", 0.0},
};

string adaptive_tempscan_root = "../../Generated content/Coulomb Constant/real coulomb/J5a10b1/stepsize05";

map<string, double> adaptive_temp_scan_standard = {
        {"end_t", 50.0},
        {"dt_max", 0.0005},
        {"J", 5},
        {"alpha", 10},
        {"beta", 1.0},
        {"tau", 1},
        {"eta", 0.8},
        {"N", 2.0},
        {"nr_save_values", 16.0},
        {"x0", 1.4},
        {"p0", 1.4},
        {"repeat_nr", 1.0},
        {"min_temp", 0.01},
        {"max_temp", 2},
        {"nr_temps", 10.0},
        {"K", 10.0},
        {"tol", 0.5},
        {"logspace", 0.0}
};


string quench_root = "";

map<string, double> quench_paras = {
        {"starting_T", 120.0},
        {"end_T", 40.0},
        {"t_eq", 50},
        {"dt", 0.002},
        {"J", 10.0},
        {"alpha", 5.0},
        {"beta", 20.0},
        {"eta", 0.8},
        {"N", 2.0},
        {"nr_save_values", 16.0},
        {"nr_xis", 32.0},
        {"repeat", 10.0},
        {"min_tau_factor", 1.0},
        {"max_tau_factor", 3.0},
        {"nr_taus", 1.0},
        {"K", 3.0},
        {"tol", 0.01}
};

template <typename Key, typename Value>
void printMap(const std::map<Key, Value>& myMap) {
    for (const auto& pair : myMap) {
        std::cout << "Key: " << pair.first << ", Value: " << pair.second << std::endl;
    }
}



/*
const int       nr_temps = 40;
const size_t    lattice_dim = 100;
const size_t    n = lattice_dim * lattice_dim;
const size_t    N = 2;
*/








#endif //CUDAPROJECT_PARAMETERS_CUH
