//
// Created by andi on 17.08.23.
//
//
// Created by andi on 15.08.23.
//
//
// Created by andi on 13.08.23.
//

#include "../Simulation.cuh"


int main() {
    // reading the parameters from the parameter file
    map<string, double> paras = performance_paras;
    int nr_save_values = (int)paras["nr_save_values"];
    fs::path simulation_path = performance_root;
    const size_t lat_dim = lattice_dim;
    // typedefs
    typedef thrust::device_vector<double> state_type;
    typedef thrust_algebra algebra;
    typedef thrust_operations operations;

    // Okay so we initialize the observer first haha
    auto* relax_obs =
            new relax_observer<anisotropic_coulomb_constant, state_type>(nr_save_values);

    auto* runtime_obs =
            new runtime_observer<anisotropic_coulomb_constant, state_type>();

    // templating..
    PerformanceSimulation simulation = PerformanceSimulation<euler_combined,
            state_type,
            algebra, operations,
            anisotropic_coulomb_constant>(paras, simulation_path);
    simulation.register_observer(relax_obs);
    simulation.register_observer(runtime_obs);
    simulation.simulate();
    return 0;
}