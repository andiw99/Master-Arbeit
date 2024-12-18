//
// Created by andi on 15.08.23.
//
//
// Created by andi on 13.08.23.
//
#include "../Simulation.cuh"



int main(int argc, char* argv[]) {
    path filepath;
    typedef XY_Silicon relax_system;
    if (argc == 2) {
        filepath = "parameters/para_set_" + (string)argv[1] + ".txt" ;
    } else {
        filepath = "parameters/para_set_0.txt";
    }

    // reading the parameters from the parameter file
    map<Parameter, double> paras = readTxtFileToParameterMap(filepath, 1);  // start at the second line
    fs::path simulation_path = readTxtFileToString(filepath);
    int nr_save_values = (int)paras[Parameter::nr_saves];
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
    RelaxationSimulation simulation = RelaxationSimulation<bbk_stepper,
            state_type,
            algebra, operations,
            relax_system>(paras, simulation_path);
    simulation.register_observer(relax_obs);
    simulation.register_observer(runtime_obs);
    simulation.register_observer(ner_obs);
    simulation.simulate();
    return 0;
}