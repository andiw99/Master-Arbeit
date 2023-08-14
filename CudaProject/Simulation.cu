//
// Created by andi on 13.08.23.
//

#include "Simulation.cuh"


int main() {
    // reading the parameters from the parameter file
    map<string, double> paras = quench_paras;
    int nr_save_values = (int)paras["nr_save_values"];
    fs::path simulation_path = quench_root;
    const size_t lat_dim = lattice_dim;
    // typedefs
    typedef thrust::device_vector<double> state_type;
    typedef thrust_algebra algebra;
    typedef thrust_operations operations;

    // Okay so we initialize the observer first haha
    observer* quench_obs = new quench_observer(nr_save_values);
    // templating..
    Simulation simulation = QuenchSimulation<   lat_dim, euler_combined,
                                                state_type,
                                                algebra, operations,
                                                anisotropic_coulomb_quench>(paras, simulation_path);
    simulation.register_observer(quench_obs);
    simulation.simulate();
    return 0;
}