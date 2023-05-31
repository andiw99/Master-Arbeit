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

string tempscan_root = "../../Generated content/Higher Temps/";

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
        {"nr_save_values", 200.0},
        {"x0", 4.0},
        {"p0", 4.0},
        {"repeat_nr", 1.0},
        {"min_temp",10.0},
        {"max_temp", 60.0},
        {"nr_temps", 10.0},
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
