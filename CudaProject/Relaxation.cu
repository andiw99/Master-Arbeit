//
// Created by andi on 15.08.23.
//
//
// Created by andi on 13.08.23.
//

#include "Simulation.cuh"
#include "parameters.cuh"


int main() {
    // reading the parameters from the parameter file
    map<Parameter, double> paras = temp_scan_paras;
    int nr_save_values = (int)paras[Parameter::nr_saves];
    fs::path simulation_path = adaptive_tempscan_root;
    // typedefs
    typedef thrust::device_vector<double> state_type;
    typedef thrust_algebra algebra;
    typedef thrust_operations operations;

    // Okay so we initialize the observer first haha
    auto* relax_obs =
            new relax_observer<relax_system, state_type>(nr_save_values);

    auto* runtime_obs =
            new runtime_observer<relax_system, state_type>();

    auto* ner_obs = new NER_observer<relax_system, state_type>(paras[nr_ner_values]);

    // templating..
    RelaxationSimulation simulation = RelaxationSimulation<euler_combined,
            state_type,
            algebra, operations,
            relax_system>(paras, simulation_path);
    simulation.register_observer(relax_obs);
    simulation.register_observer(runtime_obs);
    simulation.register_observer(ner_obs);
    simulation.simulate();
    return 0;
}