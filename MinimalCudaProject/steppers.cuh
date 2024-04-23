//
// Created by andi on 15.08.23.
//

#ifndef CUDAPROJECT_STEPPERS_CUH
#define CUDAPROJECT_STEPPERS_CUH
#include "../LearningProject/Header/Helpfunctions and Classes.h"

using namespace std;


/*
 * We start with the stepper, this stepper has to be templated since we want to exchange the container for our state type
 * We also template the 'Algebra' that outsources the iteration through my lattice sites?
 * And the 'Operation' that performs the operation on the lattice sites
 */

template<
        class state_type,
        class algebra,
        class operations,
        class System,
        class value_type = double,
        class time_type = value_type
>
class stepper {
public:
    typedef System Sys;

    ofstream debugging_file;

    stepper() {}
    stepper(size_t N) : N(N), dxdt(N), theta(N) {
        debugging_file.open("../../Generated content/debugging4");
    }
    stepper(map<Parameter, double> paras) : stepper((2 * (int)paras[total_size])){
        debugging_file.open("../../Generated content/debugging4");
    }
    // the system size, but what is N if i am in 2 dimensions? probably n * n. But how can i initialize the inital
    // values for all sites? not just by "stat_type(N)"
    // i still don't fully understand what the state_type is, is it just (x, p) for one lattice site? Or is it the
    // the state of the whole system? In the former, N would correspond with the dimensionality of my DGL system
    // in the latter i would not really know how to initialize the initial states independent of the dimensionality
    // of my problem

    virtual void reset() {
        //supposed to reset some variables, nothing to do here
    }

    const size_t N;
    state_type dxdt, theta;
/*    vector<observer*> obsvers = {};               // new, we can register an observer to a stepper

    void register_observer(observer* obs) {
        obsvers.push_back(obs);
    }*/

    virtual void do_step(Sys& sys, state_type& x, time_type dt_max, time_type& t) {
        cout << "dummy do step is called" << endl;
    }

/*    template<class Sys>
    void step_until(time_type end_time, Sys& sys, state_type& x, time_type dt_max, time_type &t) {
        for(auto obs : obsvers) {
            obs->operator()(sys, x, t); // Observing before anything happens
        }
        while (t < end_time){
            this->do_step(sys, x, dt_max, t); // it is important that the steper custom do step is called here
            for(auto obs : obsvers) {
                obs->operator()(sys, x, t); // Observing before anything happens
            }
        }
    }*/

    virtual void print_stepper_info() {
        cout << "Stepper Info:" << endl;
        cout << "N = " << N << endl;
    }


    template<class Obs>
    void step_until(time_type end_time, Sys& sys, state_type& x, time_type dt_max, time_type &t, vector<Obs*> obsvers) {
        for(auto obs : obsvers) {
            Singleton_timer::set_startpoint(obs->get_name());       // TODO i could imagine that this doesnt work since i have to pass a variable and not a string?
            obs->operator()(sys, x, t); // Observing before anything happens
            Singleton_timer::set_endpoint(obs->get_name());
        }
        cout << "Starting to step?" << endl;
        while (t < end_time){
/*            if (t > 0.95 * end_time) {
                for(int i = 0; i <= 1000; i++) {
                    debugging_file << theta[x.size() / 2 + i] << ",";
                }
                debugging_file << endl;
            }*/
            this->do_step(sys, x, dt_max, t); // it is important that the steper custom do step is called here
            //cout << endl << endl;
            for(auto obs : obsvers) {
                Singleton_timer::set_startpoint(obs->get_name());
                // cout << obs->get_name() << endl;
                obs->operator()(sys, x, t);
                Singleton_timer::set_endpoint(obs->get_name());
            }
            //cout << endl << endl;
            if(sys.is_equilibrated()) {
                break;
            }
        }
    }

    static string get_name() {
        return "base stepper";
    }
};

template<
        class state_type,
        class algebra,
        class operations,
        class System,
        class value_type = double,
        class time_type = value_type
>
class bbk_stepper : public stepper<state_type, algebra, operations, System, value_type, time_type> {
public:
    typedef stepper<state_type, algebra, operations, System, value_type, time_type> stepper;
    typedef System Sys;         // Why this typedef? ah just to shorten things i guess
    using stepper::dxdt;
    state_type F = dxdt;
    using stepper::theta;
    using stepper::N;       // this is 2 * n so the size of a normal vector
    using stepper::stepper;
    typedef typename operations::template apply_drift<time_type> apply_drift;
    typedef typename operations::template apply_bbk_v1<time_type> apply_bbk_v1;
    typedef typename operations::template apply_bbk_v2<time_type> apply_bbk_v2;
    size_t n;
    bool firststep = true;
    const string calc_force_checkpoint_name = " Calc Force";      // I think this should be optimized by the compiler??
    const string calc_diff_checkpoint_name = " Calc Diffusion";
    const string bbk_checkpoint_v1 = "BBK V1";
    const string bbk_checkpoint_v2 = "BBK V2";
    const string bbk_checkpoint_drift = "BBK apply drift";
    void do_step(Sys& sys, state_type& x, time_type dt, time_type& t) override {
        // we only have to calc F in the first run, otherwise we can reuse the F of the last step
        // also we only have to generate random numbers at the beginning if its first run
        if(firststep) {
            Singleton_timer::set_startpoint(sys.get_name() + calc_force_checkpoint_name);      // Danger to have a different string here, thats why you were binding it with variables
            sys.calc_force(x, F, t);
            Singleton_timer::set_endpoint(sys.get_name() + calc_force_checkpoint_name);

            Singleton_timer::set_startpoint(sys.get_name() + calc_diff_checkpoint_name);
            sys.calc_diff(theta, t);
            Singleton_timer::set_endpoint(sys.get_name() + calc_diff_checkpoint_name);

            firststep = false;
        }
        // TODO we should probalby also have a stepper that is not timed as this slowly seems to consume some time??
        // half a kick
        Singleton_timer::set_startpoint(bbk_checkpoint_v1);
        algebra::for_each(x.begin() + n, F.begin() + n, theta.begin() + n, n, apply_bbk_v1(dt, sys.get_eta()));
        Singleton_timer::set_endpoint(bbk_checkpoint_v1);

        // drift on q
        Singleton_timer::set_startpoint(bbk_checkpoint_drift);
        algebra::for_each(x.begin(), x.begin(), x.begin() + n, n, apply_drift(dt));
        Singleton_timer::set_endpoint(bbk_checkpoint_drift);

        sys.map_state(x);
        // now we can calculate the force again
        Singleton_timer::set_startpoint(sys.get_name() + calc_force_checkpoint_name);      // Danger to have a different string here, thats why you were binding it with variables
        sys.calc_force(x, F, t);
        Singleton_timer::set_endpoint(sys.get_name() + calc_force_checkpoint_name);

        Singleton_timer::set_startpoint(sys.get_name() + calc_diff_checkpoint_name);
        sys.calc_diff(theta, t);
        Singleton_timer::set_endpoint(sys.get_name() + calc_diff_checkpoint_name);
        // half a kick
        //                     v(n+1/2)           F(q_n+1)               R_n+1
        Singleton_timer::set_startpoint(bbk_checkpoint_v2);
        algebra::for_each(x.begin() + n, F.begin() + n, theta.begin() + n, n, apply_bbk_v2(dt, sys.get_eta()));
        Singleton_timer::set_endpoint(bbk_checkpoint_v2);

        // It is weird but it should be it i guess.
        // advance
        t += dt;
    }

    static string get_name() {
        return "bbk stepper";
    }
    void reset() override {
        //supposed to reset some variables, nothing to do here
    }
    bbk_stepper(size_t N) : stepper(N), n(N/2)
    {
    }
    bbk_stepper(map<Parameter, double> paras): bbk_stepper(2*(int)paras[total_size]){}


    void print_stepper_info() override {
        stepper::print_stepper_info();
    }
};


#endif //CUDAPROJECT_STEPPERS_CUH
