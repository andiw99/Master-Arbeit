//
// Created by andi on 13.08.23.
//

#ifndef CUDAPROJECT_SIMULATION_CUH
#define CUDAPROJECT_SIMULATION_CUH
#include "main.cuh"
#include "systems.cuh"
#include "parameters.cuh"


template <size_t lat_dim, template<class, class, class, class, class> class stepper>     // I think I don't do anything with the stepper so i should not template the base class?
class Simulation {
    observer obs;
    map<string, double> paras;
    fs::path simulation_path;
    int seed = 0;   //should the simulation be seedable?
public:
    Simulation(map<string, double> &paras, fs::path& simulation_path): paras(paras),
    simulation_path(simulation_path) {}
    Simulation(map<string, double> &paras, fs::path& simulation_path, int seed) :
    Simulation(paras, simulation_path) {
        this->seed = seed;
    }

    virtual void simulate() {}
    virtual void run() {}
    virtual void repeat() {}
    virtual void initialize() {}
};


template <size_t lat_dim, template<class, class, class, class, class> class stepper>
class QuenchSimulation : public Simulation<lat_dim, stepper> {
    using Simulation<lat_dim, stepper>::paras;
    quench<lat_dim> sys;
    vector<double> taus;

    void initialize() override {
        taus = logspace(paras["min_tau_factor"],
                                             paras["max_tau_factor"],
                                             (int)paras["nr_taus"] + 1);
    }

    void repeat(int runs) override{
        // repeat gets called just with the runs, we can look at the parameters?
        // but we need specific parameters for this repeat? I mean we can just alter the overall parameter map
        // I am thinking whether the repeat function is so general that i can put this into the base class
        // but I need to change the run for sure...
    }

    void run() {

    }

public:
    QuenchSimulation(map<string, double> &paras, fs::path& simulation_path) :
    Simulation<lat_dim, stepper>(paras, simulation_path) {
        // probably calling initialize here to set everything up
        initialize();
    }



};

#endif //CUDAPROJECT_SIMULATION_CUH
