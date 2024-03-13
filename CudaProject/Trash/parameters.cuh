//
// Created by andi on 30.05.23.
//

#ifndef CUDAPROJECT_PARAMETERS_CUH
#define CUDAPROJECT_PARAMETERS_CUH

#include <map>

// we now create an enumeration of parameters to avoid mistakes like paras["Jx"]


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

const static map<string, double> temp_scan_standard = {
        {"lat_dim", 100},
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

string adaptive_tempscan_root = "../../Generated content/BBK Test";
typedef quadratic_trapped_lattice relax_system;
map<string, double> adaptive_temp_scan_standard = {
        {"lat_dim", 32.0},
        {"end_t", 100.0},
        {"dt", 0.001},
        {"J", 3},
        {"Jy", -1},
        {"alpha", 10},
        {"beta", 1.0},
        {"tau", 1},
        {"eta", 0.5},
        {"N", 2.0},
        {"nr_save_values", 16.0},
        {"x0", 1.4},
        {"p0", 1.4},
        {"repeat", 0.0},
        {"min_temp", 9.0},
        {"max_temp", 10.0},
        {"nr_temps", 2.0},
        {"K", 10.0},
        {"tol", 0.5},
        {"logspace", 0.0},
        {"random", 1.0}
};

const map<Parameter, double> temp_scan_paras = {
        {dim_size_x, 4.0},
        {dim_size_y, 4.0},
        {end_time, 200.0},
        {dt, 0.001},
        {J, 5},
        {Jy, 1.0},
        {alpha, 20},
        {Parameter::beta, 1.0},
        {tau, 1},
        {eta, 0.5},
        {nr_saves, 16.0},
        {x0, 1.4},
        {p0, 1.4},
        {nr_repeat, 0.0},
        {min_temp, 9.0},
        {max_temp, 10.0},
        {nr_runs, 2.0},
        {K, 10.0},
        {tol, 0.5},
        {logspaced, 0.0},
        {random_init, 1.0},
        {nr_ner_values, 100},
        {m, 1},
        {p, 2},
};


string quench_root = "../../Generated content/Trash/Testing Rectangular/Quench";
// TODO Good idea would have been to make the parameters an enumeration, refactoring idea for the far future
typedef anisotropic_coulomb_quench quench_system;
map<string, double> quench_string_paras = {
        {"lat_dim", 31.0},
        {"starting_T", 1},
        {"end_T", 0.05},
        {"t_eq", 100},
        {"dt", 0.001},
        {"J", -1.0},
        {"Jy", -3.0},
        {"alpha", 10.0},
        {"beta", 1.0},
        {"x0", 1},
        {"p0", 2},
        {"eta", 1.5},
        {"N", 2.0},
        {"nr_save_values", 128.0},
        {"nr_xis", 32.0},
        {"repeat", 0.0},
        {"min_tau_factor", 7.0},
        {"max_tau_factor", 11.0},
        {"nr_taus", 1.0},
        {"K", 3.0},
        {"tol", 0.01},
        {"video", 1.0},
        {"random", 1.0}
};

map<Parameter, double> quench_paras = {
        {dim_size_x, 256.0},
        {dim_size_y, 128.0},
        {starting_temp, 1},
        {end_temp, 0.05},
        {equil_time, 100},
        {dt, 0.001},
        {J, 4.0},
        {Jy, 1.0},
        {alpha, 10.0},
        {Parameter::beta, 1.0},
        {x0, 1},
        {p0, 2},
        {eta, 1.5},
        {nr_saves, 128.0},
        {nr_repeat, 0.0},
        {min_tau_factor, 6.0},
        {max_tau_factor, 10.0},
        {nr_runs, 4.0},
        {K, 3.0},
        {tol, 0.01},
        {random_init, 1.0}
};

string performance_root  = "../../Generated content/Perf Test";
// TODO Good idea would have been to make the parameters an enumeration, refactoring idea for the far future
map<string, double> performance_paras = {
        {"end_t", 1000.0},
        {"dt", 0.001},
        {"J", -1},
        {"Jy", -3},
        {"alpha", 10.0},
        {"beta", 1.0},
        {"eta", 1.5},
        {"N", 2.0},
        {"nr_save_values", 16.0},
        {"repeat", 1.0},
        {"min_lat_factor", 11.0},
        {"max_lat_factor", 11.0},
        {"nr_lat_dims", 0},
        {"K", 5.0},
        {"tol", 0.01},
};



template<typename key, typename value_type>
void printMapOfVectors(const std::map<key, std::vector<value_type>>& myMap) {
    for (const auto& entry : myMap) {
        std::cout << "Key: " << entry.first << ", Value: [";
        for (size_t i = 0; i < entry.second.size(); ++i) {
            std::cout << entry.second[i];
            if (i < entry.second.size() - 1) {
                std::cout << ", ";
            }
        }
        std::cout << "]" << std::endl;
    }
}



/*
const int       nr_temps = 40;
const size_t    lattice_dim = 100;
const size_t    n = lattice_dim * lattice_dim;
const size_t    N = 2;
*/








#endif //CUDAPROJECT_PARAMETERS_CUH
