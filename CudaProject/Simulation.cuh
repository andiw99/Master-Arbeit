//
// Created by andi on 13.08.23.
//

#ifndef CUDAPROJECT_SIMULATION_CUH
#define CUDAPROJECT_SIMULATION_CUH
#include "main.cuh"
#include "systems.cuh"
#include "deprecated-systems.cuh"
#include "steppers.cuh"
#include "observers.cuh"


template <
        template<class, class, class, class, class, class> class stepper_type,
        class state_type, class alg, class oper,
        class sys >     // I think I don't do anything with the stepper so i should not template the base class?
class Simulation {
public:
    // i think we should work with a vector of pointers to observers, this should work right? I dont really trust me to safely
    // handle observer** pointer
    typedef obsver<sys, state_type> observer_type;
    vector<observer_type*> obsvers;          // Observer might take many parameters to be flexibel in its logic, i think it is easiest to initialize it outside of the simulation class
    map<Parameter, double> paras;
    fs::path simulation_path;
    int seed = 0;   //should the simulation be seedable?
    int n;          // current system size, but can change
    stepper_type<state_type, alg, oper, sys, double, double>* stepper;
    size_t step_nr = 0;
    fs::path folder_path;   // gets updated everytime repeat is called for the first time, depends on whether quench or not because of folder names
    int repeat_nr;
    int run_count = 0;          // stupid variable for simulations that are already present in the folder
    Simulation(map<Parameter, double> &paras, fs::path& simulation_path): paras(paras),
    simulation_path(simulation_path) {
        // Stepper cannot be created here anymore since it initializes vectors of size N * n which now can change
        // per setting
        // we need a random standard seed?
        chrono::milliseconds ms = chrono::duration_cast<chrono::milliseconds >(
                chrono::system_clock::now().time_since_epoch()
        );
        this->seed = (ms.count() % 10000) * 10000000;
    }
    Simulation(map<Parameter, double> &paras, fs::path& simulation_path, int seed) :
    Simulation(paras, simulation_path) {
        this->seed = seed;
    }
public:
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
        paras[Parameter::step_nr] = step_nr;
        paras[Parameter::run_nr] = nr;
        // 2. Initialize the stepper and the system (and open new file for the observer)
        // the stepper is already initialized, but i think it should have a reset button?
        stepper->print_stepper_info();
        stepper->reset();
        stepper->print_stepper_info();

        // sys<lat_dim> Sys = create<lat_dim, sys>(paras, step_nr);   // seed?
        // we try to create only with the map...

        state_type x(2 * n);                            // N is not a parameter anymore
        sys Sys = sys(paras);
        Sys.print_info();
        Sys.init_state(paras, x);           // init state via correct state initializer that was created in repeat

        double end_t = Sys.get_end_t();
        cout << "Starting run with parameters:" << endl;
        printMap(paras);

        for(auto obs : obsvers) {
            // calling init rather than open stream.
            obs->init(folder_path, paras, Sys);
        }
        // system? System should probably be created new since we need another tau for the system
        // reset the current T and stuff
        // 3. run the simulation haha, observing logic is handled by the observer
        // now we still need all the observing logic, damn that was more work than anticipated
        cout << "For a " << (int)paras[Parameter::dim_size_x] << " x " << (int)paras[Parameter::dim_size_y] << " System" << endl;
        double t = 0;
        stepper->step_until(end_t, Sys, x, paras[Parameter::dt], t, obsvers);
        cout << "Sys.get_step_nr():" << Sys.get_step_nr() << endl;
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
        // TODO it is important that the subroutines set the lat_dim parameter in paras
        run_count = findHighestCSVNumber(folder_path) + 1;
        if (runs == repeat_nr) {
            // i am not really happy with this if statement but it safes some redundant code which in don't know
            // where to put otherwise
            n = (int)(paras[Parameter::dim_size_x] * paras[Parameter::dim_size_y]);
            paras[Parameter::total_size] = n;
            // stepper = create_stepper<state_type, alg, oper, sys, double, double, stepper_type>(paras);
            stepper = new stepper_type<state_type, alg, oper, sys, double, double>(paras);
        }
        run(run_count);
        if(runs > 0) {
            repeat(runs-1);
        } else {
            return;
        }
    }

    virtual void initialize() {
        // initialize the stepper, again the problem that depending on the type, I need to initialize differently
        // We do the heretical is_same initialization this time because partial template specialization is not a thing in c++
        // We could theoretically also create the stepper in the base class, but then we would have
        // the huge templates everywhere, but on the other hand registering observer etc is stuff that has to be
        // done for every simulation and it would be nice to not have to write that again
        repeat_nr = (int) paras[Parameter::nr_repeat];
        step_nr = seed;
        cout << "initial step_nr = " << step_nr << endl;
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
    }

    ~Simulation() {
        // call the destroyers of the observers
        cout << "Deleting simulation" << endl;
        for(auto obs : obsvers) {
            delete obs;
        }
    }
};


template <
            template<class, class, class,class, class, class> class stepper_type,
            class state_type, class alg, class oper,
                    class quench_sys >
class QuenchSimulation : virtual public Simulation<stepper_type, state_type, alg, oper, quench_sys> {
    typedef Simulation<stepper_type, state_type, alg, oper, quench_sys> Simulation;
    using Simulation::paras;
    using Simulation::simulation_path;
    using Simulation::folder_path;
    vector<double> taus;

    void initialize() {
        // init the taus i want the simultion to run over
        taus = logspace(paras[min_tau_factor],
                                             paras[max_tau_factor],
                                             (int)paras[nr_runs] + 1);
        paras[T] = paras[starting_temp];
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
            paras[Parameter::tau] = tau;
            folder_path = simulation_path / to_string(tau);
            this->repeat((int)paras[nr_repeat]);      // and actually this should be it
        }
    }

    double get_end_t() override {
        double eq_t = paras[equil_time];
        double t_quench = (paras[starting_temp] - paras[end_temp]) * paras[tau];
        double t_end = 2 * eq_t + t_quench;
        cout << "Simulating Quench until " <<  t_end << endl;
        // set end t in paras for the observer
        paras[end_time] = t_end;
        return t_end;
    }

    QuenchSimulation(map<Parameter, double> &paras, fs::path& simulation_path) :
            Simulation(paras, simulation_path) {
    }
};

template <
        template<class, class, class, class, class, class> class stepper_type,
        class state_type, class alg, class oper,
        class relax_sys >
class RelaxationSimulation : virtual public Simulation<stepper_type, state_type, alg, oper, relax_sys> {
    typedef Simulation<stepper_type, state_type, alg, oper, relax_sys> Simulation;
    using Simulation::paras;
    using Simulation::simulation_path;
    using Simulation::folder_path;
    vector<double> temps;

    void initialize() override {
        // init the taus i want the simultion to run over
        temps = linspace(paras[min_temp],
                             paras[max_temp], (int)paras[nr_runs] + 1);
        // call the general initialization
        cout << "Temperatures:" << endl;
        print_vector(temps);
        Simulation::initialize();
    }

public:

    void simulate() {
        // we need to initialize here or think of a different architecture since
        // if we initialize in the constructor we do not register the observers to the stepper
        cout << "Calling RelaxationSimulation simulate" << endl;
        RelaxationSimulation::initialize();
        // runs the whole simulation, the sequence depends on the type of simulation we do
        for(auto T : temps) {
            // update the parameters so we create the correct systems and stuff
            paras[Parameter::T] = T;
            folder_path = simulation_path / to_string(T);
            this->repeat((int)paras[nr_repeat]);      // and actually this should be it
        }
    }

    double get_end_t() override {
        double t_end = paras[end_time];
        cout << "Simulating Relaxation until " << t_end << endl;
        return t_end;
    }
    RelaxationSimulation(map<Parameter, double> &paras, fs::path& simulation_path) :
            Simulation(paras, simulation_path) {
    }
};

template <
        template<class, class, class, class, class, class> class stepper_type,
        class state_type, class alg, class oper,
        class relax_sys >
class SubsystemRelaxationSimulation : public RelaxationSimulation<stepper_type, state_type, alg, oper, relax_sys> {
    typedef Simulation<stepper_type, state_type, alg, oper, relax_sys> Simulation;
    typedef RelaxationSimulation<stepper_type, state_type, alg, oper, relax_sys> RelaxationSimulation;
    using Simulation::paras;
    using Simulation::simulation_path;

    vector<pair<size_t, size_t>> subsystem_sizes{};
    double x_y_factor;
    path rootpath;

    void initialize() override {
        vector<size_t> Lx_s = linspace((size_t)paras[subsystem_min_Lx], (size_t)paras[subsystem_max_Lx], (int)paras[nr_subsystem_sizes] + 1);
        cout << "Sizes Lx:" << endl;
        print_vector(Lx_s);
        x_y_factor = paras[Parameter::x_y_factor];
        for(size_t Lx : Lx_s) {
            auto Ly = (size_t) ((double) Lx * x_y_factor);
            subsystem_sizes.push_back(pair<size_t, size_t>(Lx, Ly));
        }
        // RelaxationSimulation::initialize(); we dont even need to call that since the Relaxation Simulation calls it anyway?
    }
public:
    void simulate() {
        // we need to initialize here or think of a different architecture since
        // if we initialize in the constructor we do not register the observers to the stepper
        cout << "Calling SubsystemRelaxationSimulation simulate" << endl;
        initialize();
        // runs the whole simulation, the sequence depends on the type of simulation we do
        int nr_subsystems = (int)paras[Parameter::nr_subsystems];
        for(auto L_pair : subsystem_sizes) {
            // update the parameters so we create the correct systems and stuff
            paras[Parameter::dim_size_y] = L_pair.second;
            paras[Parameter::dim_size_x] = L_pair.first * nr_subsystems;
            paras[Parameter::subsystem_Ly] = L_pair.second;
            paras[Parameter::subsystem_Lx] = L_pair.first;
            simulation_path = rootpath / to_string(L_pair.first);
            RelaxationSimulation::simulate();
        }
    }
    SubsystemRelaxationSimulation(map<Parameter, double> &paras, fs::path& simulation_path) :
            RelaxationSimulation(paras, simulation_path), Simulation(paras, simulation_path), rootpath(simulation_path) {
    }

};

template <
        template<class, class, class, class, class, class> class stepper_type,
        class state_type, class alg, class oper,
        class relax_sys >
class SubsystemSimulation : virtual public Simulation<stepper_type, state_type, alg, oper, relax_sys> {
    typedef Simulation<stepper_type, state_type, alg, oper, relax_sys> Simulation;
    using Simulation::paras;
    using Simulation::simulation_path;

public:
    vector<pair<size_t, size_t>> subsystem_sizes{};
    double x_y_factor;

    path rootpath;
    void initialize_subsystems() {
        vector<size_t> Lx_s = linspace((size_t)paras[subsystem_min_Lx], (size_t)paras[subsystem_max_Lx], (int)paras[nr_subsystem_sizes] + 1);
        cout << "Sizes Lx:" << endl;
        print_vector(Lx_s);
        x_y_factor = paras[Parameter::x_y_factor];
        for(size_t Lx : Lx_s) {
            auto Ly = (size_t) ((double) Lx * x_y_factor);
            subsystem_sizes.push_back(pair<size_t, size_t>(Lx, Ly));
        }
        // RelaxationSimulation::initialize(); we dont even need to call that since the Relaxation Simulation calls it anyway?
    }
    void simulate_subsystems() {
        // we need to initialize here or think of a different architecture since
        // if we initialize in the constructor we do not register the observers to the stepper
        cout << "Calling SubsystemRelaxationSimulation simulate" << endl;
        initialize_subsystems();
        // runs the whole simulation, the sequence depends on the type of simulation we do
        int nr_subsystems = (int)paras[Parameter::nr_subsystems];
        for(auto L_pair : subsystem_sizes) {
            // update the parameters so we create the correct systems and stuff
            paras[Parameter::dim_size_y] = L_pair.second;
            paras[Parameter::dim_size_x] = L_pair.first * nr_subsystems;
            paras[Parameter::subsystem_Ly] = L_pair.second;
            paras[Parameter::subsystem_Lx] = L_pair.first;
            simulation_path = rootpath / to_string(L_pair.first);
            this->simulate();
        }
    }
    SubsystemSimulation(map<Parameter, double> &paras, fs::path& simulation_path) :
            rootpath(simulation_path), Simulation(paras, simulation_path) {
    }

};

/*
template <
        template<class, class, class, class, class, class> class stepper_type,
        class state_type, class alg, class oper,
        class relax_sys >
class SubsystemQuenchSimulation : public QuenchSimulation<stepper_type, state_type, alg, oper, relax_sys> {
    typedef Simulation<stepper_type, state_type, alg, oper, relax_sys> Simulation;
    typedef QuenchSimulation<stepper_type, state_type, alg, oper, relax_sys> QuenchSimulation;
    using Simulation::paras;
    using Simulation::simulation_path;

    vector<pair<size_t, size_t>> subsystem_sizes{};
    double x_y_factor;
    path rootpath;

    void initialize() override {
        vector<size_t> Lx_s = linspace((size_t)paras[subsystem_min_Lx], (size_t)paras[subsystem_max_Lx], (int)paras[nr_subsystem_sizes] + 1);
        cout << "Sizes Lx:" << endl;
        print_vector(Lx_s);
        x_y_factor = paras[Parameter::x_y_factor];
        for(size_t Lx : Lx_s) {
            auto Ly = (size_t) ((double) Lx * x_y_factor);
            subsystem_sizes.push_back(pair<size_t, size_t>(Lx, Ly));
        }
        // RelaxationSimulation::initialize(); we dont even need to call that since the Relaxation Simulation calls it anyway?
    }
public:
    void simulate() {
        // we need to initialize here or think of a different architecture since
        // if we initialize in the constructor we do not register the observers to the stepper
        cout << "Calling SubsystemRelaxationSimulation simulate" << endl;
        initialize();
        // runs the whole simulation, the sequence depends on the type of simulation we do
        int nr_subsystems = (int)paras[Parameter::nr_subsystems];
        for(auto L_pair : subsystem_sizes) {
            // update the parameters so we create the correct systems and stuff
            paras[Parameter::dim_size_y] = L_pair.second;
            paras[Parameter::dim_size_x] = L_pair.first * nr_subsystems;
            paras[Parameter::subsystem_Ly] = L_pair.second;
            paras[Parameter::subsystem_Lx] = L_pair.first;
            simulation_path = rootpath / to_string(L_pair.first);
            QuenchSimulation::simulate();
        }
    }
    SubsystemQuenchSimulation(map<Parameter, double> &paras, fs::path& simulation_path) :
            QuenchSimulation(paras, simulation_path), rootpath(simulation_path) {
    }

};
*/

template <
        template<class, class, class, class, class, class> class stepper_type,
        class state_type, class alg, class oper,
        class relax_sys >
class SubsystemQuenchSimulation : public QuenchSimulation<stepper_type, state_type, alg, oper, relax_sys>,
                                    public SubsystemSimulation<stepper_type, state_type, alg, oper, relax_sys> {
    typedef Simulation<stepper_type, state_type, alg, oper, relax_sys> Simulation;
    typedef QuenchSimulation<stepper_type, state_type, alg, oper, relax_sys> QuenchSimulation;
    typedef SubsystemSimulation<stepper_type, state_type, alg, oper, relax_sys> SubsystemSimulation;
public:
    using SubsystemSimulation::initialize;
    SubsystemQuenchSimulation(map<Parameter, double>& paras, fs::path& simulation_path):    QuenchSimulation(paras, simulation_path),
                                                                                            SubsystemSimulation(paras, simulation_path),
                                                                                            Simulation(paras, simulation_path) {}
};


template <
        template<class, class, class, class, class, class> class stepper_type,
        class state_type, class alg, class oper,
        class relax_sys >
class PerformanceSimulation : public Simulation<stepper_type, state_type, alg, oper, relax_sys> {
    typedef Simulation<stepper_type, state_type, alg, oper, relax_sys> Simulation;
    using Simulation::paras;
    using Simulation::simulation_path;
    using Simulation::folder_path;
    using Simulation::Simulation;
    vector<size_t> lat_dims;
    void initialize() {
        // init the taus i want the simultion to run over
        if(paras[logspaced] == 1.0) {
            lat_dims = logspace<size_t>(paras[min_lat_factor],
                                        paras[max_lat_factor],
                                        (int)paras[nr_runs] + 1);

        } else {
            lat_dims = linspace<size_t>(paras[min_lat_factor],
                                        paras[max_lat_factor],
                                        (int)paras[nr_runs] + 1);
        }
        // call the general initialization
        Simulation::initialize();
    }

public:

    void simulate() {
        // we need to initialize here or think of a different architecture since
        // if we initialize in the constructor we do not register the observers to the stepper
        this->initialize();
        // runs the whole simulation, the sequence depends on the type of simulation we do
        for(auto lat_dim : lat_dims) {
            // update the parameters so we create the correct systems and stuff
            paras[dim_size_x] = lat_dim;
            paras[dim_size_y] = lat_dim;
            folder_path = simulation_path / to_string(lat_dim);
            this->repeat((int)paras[nr_repeat]);      // and actually this should be it
        }
    }
    double get_end_t() override {
        double t_end = paras[end_time];
        cout << "Simulating until " << t_end << endl;
        return t_end;
    }
};


#endif //CUDAPROJECT_SIMULATION_CUH
