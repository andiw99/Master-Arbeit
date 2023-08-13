//
// Created by andi on 13.08.23.
//

#ifndef CUDAPROJECT_SIMULATION_CUH
#define CUDAPROJECT_SIMULATION_CUH
#include "main.cuh"
#include "systems.cuh"
#include "parameters.cuh"


template <size_t lat_dim>     // I think I don't do anything with the stepper so i should not template the base class?
class Simulation {
    observer obs;
    map<string, double> paras;
    state_initializer* State_initializer;
    fs::path simulation_path;
    int seed = 0;   //should the simulation be seedable?
    int N;
    int n = lat_dim * lat_dim;
public:
    Simulation(map<string, double> &paras, fs::path& simulation_path): paras(paras),
    simulation_path(simulation_path), N((int)paras["N"]) {}
    Simulation(map<string, double> &paras, fs::path& simulation_path, int seed) :
    Simulation(paras, simulation_path) {
        this->seed = seed;
    }

    virtual void simulate() {
        // I probably will have some repeated code for simulate, since actually i only iterate over a tau or
        // T vector and call repeat, rest is handled by specific run implementation
    }
    virtual void run(int nr) {
        // We call the run with its number inside repeat to just know how to call our file
    }
    void repeat(int runs) {
        // repeat gets called just with the runs, we can look at the parameters?
        // but we need specific parameters for this repeat? I mean we can just alter the overall parameter map
        // I am thinking whether the repeat function is so general that i can put this into the base class
        // but I need to change the run for sure...
        // If the number of runs equals the repeat nr, which means we called repeat the first time with the current
        // parameters, we initialize a new state initializer
        if (runs == (int) paras["repeat_nr"]) {
            // i am not really happy with this if statement but it safes some redundant code which in don't know
            // where to put otherwise
            State_initializer = create_state_initializer((int)paras["random"], paras, simulation_path);
        }
        this->run();
        if(runs > 0) {
            repeat(runs-1);
        } else {
            return;
        }
    }
    virtual void initialize() {}
};


template <  size_t lat_dim,
            template<class, class, class, class, class> class stepper_type,
            class state_type, class alg, class oper,
            template<size_t> class quench_sys >
class QuenchSimulation : public Simulation<lat_dim> {
    using Simulation<lat_dim>::paras;
    using Simulation<lat_dim>::N;
    using Simulation<lat_dim>::n;
    using Simulation<lat_dim>::State_initializer;
    quench_sys<lat_dim> sys;
    vector<double> taus;
    fs::path folder_path;   // gets updated everytime repeat is called for the first time, depends on whether quench or not because of folder names
    stepper_type<state_type, alg, oper, double, double> stepper;
    size_t step_nr = 0;
    void initialize() override {
        // init the taus i want the simultion to run over
        taus = logspace(paras["min_tau_factor"],
                                             paras["max_tau_factor"],
                                             (int)paras["nr_taus"] + 1);
        // initialize the stepper, again the problem that depending on the type, I need to initialize differently
        // We do the heretical is_same initialization this time because partial template specialization is not a thing in c++
        stepper = create_stepper<state_type, alg, oper, double, double, stepper_type>(paras, n);
    }

    void run(int nr) {
        // okay what do we need to do for ever run?

        // 1. Initialize the state, either random, based on last run or in lowest energy state
        // How do we determine which? Is coded in random param, ... if else function is pretty lame again
        // The state initialization should also be exchangeable, or rather i mean it is code that is used by every
        // simulation, so we might have a 'state initializer' object as member of the base class
        state_type x(N * n);
        State_initializer->init_state(x);           // init state via correct state initializer that was created in repeat
        // 2. Initialize the stepper and the system
        // the stepper is already initialized, but i think it should have a reset button?
        stepper->reset();
        sys = create<lat_dim, quench_sys>(paras, step_nr);   // seed?
        // system? System should probably be created new since we need another tau for the system
        // reset the current T and stuff
        // 3. run the simulation haha, observing logic is handled by the observer
        // now we still need all the observing logic, damn that was more work than anticipated
        double end_t = sys.get_end_t();
        double t = 0;
        stepper.step_until(end_t, sys, x, paras["dt"], t);
        step_nr += sys.get_step_nr();
    }

public:
    QuenchSimulation(map<string, double> &paras, fs::path& simulation_path) :
    Simulation<lat_dim>(paras, simulation_path) {
        // probably calling initialize here to set everything up
        initialize();
    }
};

#endif //CUDAPROJECT_SIMULATION_CUH
