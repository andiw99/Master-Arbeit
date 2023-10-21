//
// Created by andi on 13.08.23.
//

#include "Simulation.cuh"
#include "parameters.cuh"

int main() {
    // reading the parameters from the parameter file
    map<Parameter, double> paras = quench_paras;
    int nr_save_values = (int)paras[nr_saves];
    fs::path simulation_path = quench_root;
    // typedefs
    typedef thrust::device_vector<double> state_type;
    typedef thrust_algebra algebra;
    typedef thrust_operations operations;

    // We need new observers for the standard relaxation, i mean it basically does the same but the
    // end_T reading is different
    auto* quench_obs =
            new quench_observer<quench_system, state_type>(nr_save_values);
    auto* runtime_obs =
            new runtime_observer<quench_system, state_type>();
/*    quench_observer* quench_obs =
            new quench_observer(nr_save_values);*/
    // templating..
    QuenchSimulation simulation = QuenchSimulation<euler_combined,
                                                state_type,
                                                algebra, operations,
                                                quench_system>(paras, simulation_path);
    simulation.register_observer(quench_obs);
    simulation.register_observer(runtime_obs);
    simulation.simulate();
    return 0;
}