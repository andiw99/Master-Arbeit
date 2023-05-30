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

map<string, double> temp_scan_standard = {
        {"steps", 100000.0},
        {"end_t", 200.0},
        {"dt_max", 0.001},
        {"dt_start", 0.0001},
        {"t_relax", 50},
        {"nr_temps", 40},
        {"J", 2},
        {"alpha", 1},
        {"beta", 10},
        {"tau", 10},
        {"eta", 0.2},
        {"N", 2},
        {"nr_save_values", 32},
        {"x0", 8.0},
        {"p0", 8.0},
        {"repeat_nr", 5},
        {"lattice_dim", 100}
};



/*
const int       nr_temps = 40;
const size_t    lattice_dim = 100;
const size_t    n = lattice_dim * lattice_dim;
const size_t    N = 2;
*/








#endif //CUDAPROJECT_PARAMETERS_CUH
