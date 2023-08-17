//
// Created by andi on 13.08.23.
//

#ifndef CUDAPROJECT_SIMULATION_CUH
#define CUDAPROJECT_SIMULATION_CUH
#include "main.cuh"
#include "systems.cuh"
#include "steppers.cuh"
#include "observers.cuh"
#include "parameters.cuh"


template <size_t lat_dim,
        template<size_t, class, class, class, template<size_t> class, class, class> class stepper_type,
        class state_type, class alg, class oper,
        template<size_t> class sys >     // I think I don't do anything with the stepper so i should not template the base class?
class Simulation {
public:
    // i think we should work with a vector of pointers to observers, this should work right? I dont really trust me to safely
    // handle observer** pointer
    typedef obsver<lat_dim, sys, state_type> observer_type;
    vector<observer_type*> obsvers;          // Observer might take many parameters to be flexibel in its logic, i think it is easiest to initialize it outside of the simulation class
    map<string, double> paras;
    state_initializer* State_initializer;
    fs::path simulation_path;
    int seed = 0;   //should the simulation be seedable?
    int N;
    int n = lat_dim * lat_dim;
    stepper_type<lat_dim, state_type, alg, oper, sys, double, double>* stepper;
    size_t step_nr = 0;
    fs::path folder_path;   // gets updated everytime repeat is called for the first time, depends on whether quench or not because of folder names
    int repeat_nr;
    Simulation(map<string, double> &paras, fs::path& simulation_path): paras(paras),
    simulation_path(simulation_path), N((int)paras["N"]) {
        stepper = create_stepper<lat_dim, state_type, alg, oper, sys, double, double, stepper_type>(paras, n);
    }
    Simulation(map<string, double> &paras, fs::path& simulation_path, int seed) :
    Simulation(paras, simulation_path) {
        this->seed = seed;
    }

    virtual void simulate() {
        // I probably will have some repeated code for simulate, since actually i only iterate over a tau or
        // T vector and call repeat, rest is handled by specific run implementation
        cout << "calling virtual function simulate" << endl;
    }

    virtual double get_end_t() {
        // this depends on the kind of simualtion, here we return dummy value? or nothing
        // it is either decoded in the system or in the parameters (actually both times we can get it out of the parameters...
        // but we would have to write another function for that...
        cout << "get end t of base Simualtion is called (WRONG!)" << endl;

    }

    void run(int nr) {
        // okay what do we need to do for ever run?
        // TODO okay actually all the stuff i have here is also done for every simulation.
        // What is actually different from run to run? In the Binder stuff I change the constant temperature, in
        // the quench stuff I change the Quench timescale. Then I just init the system and run it.
        // The system knows for it self whether it has to change its temperature and stuff and the observer
        // deals with where to write and what. What should be different is the folderpath construction which can
        // be done in simulate() and the for loop in simulate (which requires proper initialization). Isn't all the rest
        // the same or am I overlooking something?
        // 1. Initialize the state, either random, based on last run or in lowest energy state
        // How do we determine which? Is coded in random param, ... if else function is pretty lame again
        // The state initialization should also be exchangeable, or rather i mean it is code that is used by every
        // simulation, so we might have a 'state initializer' object as member of the base class
        state_type x(N * n);
        State_initializer->init_state(x);           // init state via correct state initializer that was created in repeat
        // 2. Initialize the stepper and the system (and open new file for the observer)
        // the stepper is already initialized, but i think it should have a reset button?
        stepper->reset();
        // sys<lat_dim> Sys = create<lat_dim, sys>(paras, step_nr);   // seed?
        // we try to create only with the map...
        paras["step_nr"] = step_nr;
        paras["run_nr"] = nr;
        sys<lat_dim> Sys = sys<lat_dim>(paras);
        for(auto obs : obsvers) {
            // calling init rather than open stream.
            obs->init(folder_path, paras, Sys);         // nah this is bad, i might have multiple observers that want to write to different files
        }
        // system? System should probably be created new since we need another tau for the system
        // reset the current T and stuff
        // 3. run the simulation haha, observing logic is handled by the observer
        // now we still need all the observing logic, damn that was more work than anticipated
        double end_t = this->get_end_t();
        double t = 0;
        cout << "I am at least trying to call the stepper?" << endl;
        stepper->step_until(end_t, Sys, x, paras["dt"], t, obsvers);
        step_nr += Sys.get_step_nr();
    }
    void repeat(int runs) {
        cout << "run " << repeat_nr - runs << "/" << repeat_nr << endl;
        // repeat gets called just with the runs, we can look at the parameters?
        // but we need specific parameters for this repeat? I mean we can just alter the overall parameter map
        // I am thinking whether the repeat function is so general that i can put this into the base class
        // but I need to change the run for sure...
        // If the number of runs equals the repeat nr, which means we called repeat the first time with the current
        // parameters, we initialize a new state initializer
        if (runs == repeat_nr) {
            // i am not really happy with this if statement but it safes some redundant code which in don't know
            // where to put otherwise
            State_initializer = create_state_initializer((int)paras["random"], paras, simulation_path);
        }
        run(repeat_nr - runs);
        if(runs > 0) {
            repeat(runs-1);
        } else {
            return;
        }
    }
    void initialize() {
        // initialize the stepper, again the problem that depending on the type, I need to initialize differently
        // We do the heretical is_same initialization this time because partial template specialization is not a thing in c++
        // We could theoretically also create the stepper in the base class, but then we would have
        // the huge templates everywhere, but on the other hand registering observer etc is stuff that has to be
        // done for every simulation and it would be nice to not have to write that again
        repeat_nr = (int) paras["repeat_nr"];
        paras["n"] = lat_dim * lat_dim;
        // the observer is supposed to be already alive here, so we just register it to the stepper
/*        for(auto obs : obsvers) {
            cout << "Registering observer to stepper: obs name:" << endl;
            cout << obs->get_name() << endl;
            stepper->register_observer(obs);
        }*/
    }

    template <class Obs>
    void register_observer(Obs* obs) {
        // this function needs to be called from outside and I have to register the custom initialized observers
        cout << obs->get_name() << endl;
        obsvers.push_back(obs);
        cout << "Type of first obs in obsvers:" << endl;
        cout << obsvers[0]->get_name() << endl;
    }

    ~Simulation() {
        // call the destroyers of the observers
        cout << "Deleting simulation" << endl;
        for(auto obs : obsvers) {
            delete obs;
        }
    }
};


template <  size_t lat_dim,
            template<size_t, class, class, class,template<size_t> class, class, class> class stepper_type,
            class state_type, class alg, class oper,
            template<size_t> class quench_sys >
class QuenchSimulation : public Simulation<lat_dim, stepper_type, state_type, alg, oper, quench_sys> {
    typedef Simulation<lat_dim, stepper_type, state_type, alg, oper, quench_sys> Simulation;
    using Simulation::paras;
    using Simulation::simulation_path;
    using Simulation::folder_path;
    vector<double> taus;

    void initialize() {
        // init the taus i want the simultion to run over
        taus = logspace(paras["min_tau_factor"],
                                             paras["max_tau_factor"],
                                             (int)paras["nr_taus"] + 1);
        // call the general initialization
        Simulation::initialize();
    }

public:

    void simulate() {
        // we need to initialize here or think of a different architecture since
        // if we initialize in the constructor we do not register the observers to the stepper
        this->initialize();
        // runs the whole simulation, the sequence depends on the type of simulation we do
        for(auto tau : taus) {
            // update the parameters so we create the correct systems and stuff
            paras["tau"] = tau;
            folder_path = simulation_path / to_string(tau);
            this->repeat((int)paras["repeat_nr"]);      // and actually this should be it
        }
    }

    double get_end_t() override {
        double eq_t = paras["t_eq"];
        double t_quench = (paras["starting_T"] - paras["end_T"]) * paras["tau"];
        double t_end = 2 * eq_t + t_quench;
        cout << "Simulating Quench until " <<  t_end << endl;
        return t_end;
    }

    QuenchSimulation(map<string, double> &paras, fs::path& simulation_path) :
            Simulation(paras, simulation_path) {
    }
};

template <  size_t lat_dim,
        template<size_t, class, class, class,template<size_t> class, class, class> class stepper_type,
        class state_type, class alg, class oper,
        template<size_t> class relax_sys >
class RelaxationSimulation : public Simulation<lat_dim, stepper_type, state_type, alg, oper, relax_sys> {
    typedef Simulation<lat_dim, stepper_type, state_type, alg, oper, relax_sys> Simulation;
    using Simulation::paras;
    using Simulation::simulation_path;
    using Simulation::folder_path;
    vector<double> temps;

    void initialize() {
        // init the taus i want the simultion to run over
        temps = linspace(paras["min_temp"],
                             paras["max_temp"], (int)paras["nr_temps"] + 1);
        // call the general initialization
        Simulation::initialize();
    }

public:

    void simulate() {
        // we need to initialize here or think of a different architecture since
        // if we initialize in the constructor we do not register the observers to the stepper
        this->initialize();
        // runs the whole simulation, the sequence depends on the type of simulation we do
        for(auto T : temps) {
            // update the parameters so we create the correct systems and stuff
            paras["T"] = T;
            folder_path = simulation_path / to_string(T);
            this->repeat((int)paras["repeat_nr"]);      // and actually this should be it
        }
    }

    double get_end_t() override {
        double t_end = paras["end_t"];
        cout << "Simulating until " << t_end << endl;
        return t_end;
    }
    RelaxationSimulation(map<string, double> &paras, fs::path& simulation_path) :
            Simulation(paras, simulation_path) {
    }
};


#endif //CUDAPROJECT_SIMULATION_CUH
