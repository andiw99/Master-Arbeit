//
// Created by andi on 15.08.23.
//
//
// Created by andi on 13.08.23.
//
#include "Simulation-cuda.cuh"



int main(int argc, char* argv[]) {
    path filepath;
    typedef XY_silicon_anisotrop_subsystems relax_system;
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
            new equilibration_observer<relax_system, state_type>();
    auto* corr_obs = new corr_equilibration_observer_adaptive<relax_system, state_type>(paras[nr_corr_values], paras[min_corr_nr],
                                                                               paras[corr_write_density], paras[equil_cutoff]);
    // As long as we dont have a density ft_observer that does not need the quench methods we dont need an ft observer at all
    //auto* ft_obs = new ft_observer<relax_system, state_type>(paras[nr_ft_values]);

    // templating..
    SubsystemRelaxationSimulation simulation = SubsystemRelaxationSimulation<bbk_stepper,
            state_type,
            algebra, operations,
            relax_system>(paras, simulation_path);
    simulation.register_observer(corr_obs);
    simulation.register_observer(relax_obs);
    simulation.simulate();
    return 0;
}