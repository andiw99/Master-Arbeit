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
class euler_mayurama_stepper : public stepper<state_type, algebra, operations, System, value_type, time_type> {
public:
    typedef System Sys;
    typedef stepper<state_type, algebra, operations, System, value_type, time_type> stepper;
    using stepper::dxdt;
    using stepper::theta;
    using stepper::N;
    using stepper::stepper;
    // observer* Observer;
    // the stepper needs a do_step method
    // I think our do step method needs an additional parameter for theta? Maybe not, we will see
    // We also template our do_step method to work with any system that we feed into it
    // i think later the system is called to perform the operation on the state types


    void do_step(Sys& sys, state_type& x, time_type dt, time_type& t) override {
        // okay we don't need this for_each3 stuff, we only need to apply x_(n+1) = x_n + k_1 dt + theta sqrt(dt)
        // first we need to calculate the derivative with the system and save it to temporary class things
        // i don't even think that we need dxdt as parameter

        // We call the system to calculate dxdt, should only calculate the deterministic part here?
        // No i would say we also calculate the stochastic part
        // so we give the system the state and istruct it later to save its calculations in dxdt and theta so that we
        // can later iterate over the lattice and apply the operation, meaning the update
        // timer.set_startpoint(system_name);
        // timer.set_startpoint(functor);
        sys.calc_drift(x, dxdt, t);
        //timer.set_endpoint(functor);
        //timer.set_startpoint(rng);
        sys.calc_diff(theta, t);
        // timer.set_endpoint(rng);
        // timer.set_endpoint(system_name);
        // this should set the correct values for dxdt and theta so that they can be applied in apply_em
        // can we print them here?

//        cout << "x[0] = " << x[0] << ", " << "x[1] = " << x[1] << endl;
//        cout << "dxdt[0] = " << dxdt[0] << ", " << "dxdt[1] = " << dxdt[1] << endl;
//        cout << "theta[0] = " << theta[0] << ", " << "theta[1] = " << theta[1] << endl;

        // for the update we need to define a type of a function, but i don't really understand the syntax
        // okay so the plan is to have a new name for the templated struct apply_em that is a member of operations
        // so that we can later more clearly instantiate the struct apply_em with the correct template
        // timer.set_startpoint(applying_name);
        typedef typename operations::template apply_em<time_type> apply_em;
        // this is the operation that i want to apply on every lattice site
        // cout << endl << endl;
        // cout << "THETA IN EULER MAYURAMA:" << endl;
        // cout << endl << endl;
        // cout << "Before application x = " << endl;
        // cout << endl;
        algebra::for_each(x, x, dxdt, theta, apply_em(dt));
        sys.map_state(x);
        // cout << "After application" << endl;
        // cout << endl;
        // timer.set_endpoint(applying_name);
        // and that should already be it?
        // Observe? How much time does this take?
        // so i also want to write down the temperature of the step, but i don't want to give it to the observer
        // since not every System has a temperature or at least one that is relevant
        // we could give the system to the observer, this would be possible, so that we can use sys.getT in a
        // special observer for the systems with temperature
        // Observer->operator()(sys, x, t);
        t += dt;
    }

    // i am currently not sure what parameters we additionally need, we don't have temporary x values like for the
    // runge kutta scheme, at least the system size should be a parameter i guess
    // I don't really get how this stuff is instantiated here
    euler_mayurama_stepper(size_t N) : stepper(N) //, Observer(Obs)
    {
        cout << "initializing euler mayurama stepper" << endl;
    }
    euler_mayurama_stepper(map<Parameter, double> paras): euler_mayurama_stepper(2*(int)paras[total_size]){}
/*    euler_mayurama_stepper(size_t N) : N(N) {
        Observer = new observer();
    }*/

    // now for the memory allocation. I am not sure if this ever changes for me but the implementation should not harm
    // on the other hand, my class doesn't have any saved state_types that could potentially be resized
    // so we skip this for now

    void print_stepper_info() override {
        stepper::print_stepper_info();
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
class euler_simple_adaptive : public stepper<state_type, algebra, operations, System, value_type, time_type>{
    int k;
    value_type tol;
    value_type error = 0;
    time_type dt;
    typedef typename operations::template apply_drift<time_type> apply_drift;
    typedef typename operations::calc_error calc_error;
    typedef typename operations::template sum<value_type> sum;
    typedef typename operations::template apply_diff<time_type> apply_diff;
    string drift_calc = "Second drift calc";
    string error_calc = "Error Calculation";
    string repetitions = "Repetitions";
    checkpoint_timer timer{};
public:
    typedef stepper<state_type, algebra, operations, System, value_type, time_type> stepper;
    typedef System Sys;
    using stepper::stepper;
    using stepper::dxdt;
    using stepper::theta;
    using stepper::N;

    // we now pass dt by reference, so that we can modify it
    void do_step(Sys& sys, state_type& x, time_type dt_max, time_type &t) override {

        // calc the stepsize with the current k
        dt = 1.0 / pow(2, k) * dt_max;

        // good thing is that we already split the drift and the diffusion
        // we are at t and can easily apply the system operations without thinking of it for now
        sys.calc_drift(x, dxdt, t);
        // Here i think we have to apply only dxdt for now to calculate f(x*) - f(x)
        // which means we need a new operations structure?
        algebra::for_each(x_drift, x, dxdt, apply_drift(dt));

        // we have x_drift now, now we need to calculate f(x_drift)
        // we really should not just call the system since the system will generate random numbers
        // so we just call again calc drift, but we need to store the result somewhere else than dxdt since
        // dxdt= f(x)
        timer.set_startpoint(drift_calc);
        sys.calc_drift(x_drift, dx_drift_dt, t);
        timer.set_endpoint(drift_calc);
        // now we need to calculate the difference between dx_drift_dt = f(x_drift) and dxdt=f(x)
        // how do we do that? we actually for a simple case just need to calc the difference for every
        // entry of dxdt and dx_drift_dt and then sum it up / average it. this should actually be a very simple
        // thrust operation
        timer.set_startpoint(error_calc);
        algebra::for_each(dx_drift_dt, dxdt, calc_error());
        // now the error is in dx_drift_dt, now we got to sum and average it
        error = sum()(dx_drift_dt) / N;
        timer.set_endpoint(error_calc);
        // now we have to check whether the error is small enough
        if(error < tol) {
            // if error is smaller than the tolerance we apply everything
            sys.calc_diff(theta, t);
            algebra::for_each(x, x_drift, theta, apply_diff(dt));
            // we also increase the time
            t += dt;
            // and we reduce k for the next step
            // do we have to check that k does not get smaller than zero?
            // TODO is there a faster way to do this than with if?
            if (k > 0) {
                k--;
            }
        } else {
            // everytime we enter else we failed the error test, so if we time this, we now how long the repetitions
            // took
            timer.set_startpoint(repetitions);
            // if the error is to large, we have to do the step again with increased k
            k++;
            do_step(sys, x, dt_max, t);
            timer.set_endpoint(repetitions);
            // I thianak thats it?
        }

    }

    euler_simple_adaptive(size_t N, int K, double tol) : dx_drift_dt(N), x_drift(N), k(K), tol(tol),
                                                         stepper(N){}

    int get_k() {
        return k;
    }

    double get_error() {
        return error;
    }


    void print_stepper_info() override {
        stepper::print_stepper_info();
        cout << "k = " << k;
    }

private:
    state_type x_drift, dx_drift_dt;
};

template<
        class state_type,
        class algebra,
        class operations,
        class System,
        class value_type = double,
        class time_type = value_type
>
class euler_combined : public euler_mayurama_stepper<state_type, algebra, operations, System, value_type, time_type>{
    int k;
    int prev_accepted_k;
    const int first_k;
    value_type tol;
    value_type error = 0;
    time_type dt;
    int switch_counter = 0;
    int switch_count;
    double reduction_factor;
    bool switched = false;
    typedef typename operations::template apply_drift<time_type> apply_drift;
    typedef typename operations::calc_error calc_error;
    typedef typename operations::template sum<value_type> sum;
    typedef typename operations::template apply_diff<time_type> apply_diff;
    typedef euler_mayurama_stepper<state_type, algebra, operations, System, value_type, time_type> euler_mayurama_stepper;
    using euler_mayurama_stepper::theta;
    using euler_mayurama_stepper::N;
    using euler_mayurama_stepper::dxdt;
/*    using do_euler_step = euler_mayurama_stepper<state_type, algebra, operations,
            value_type, time_type>::do_step;*/
    string drift_calc = "Second drift calc";
    string error_calc = "Error Calculation";
    string repetitions = "Repetitions";
    string euler_steps = "Euler-Mayurama Steps";
    string adaptive_steps = "Adaptive Steps";
    string rng = "Random Number Generation";
    checkpoint_timer timer{};
public:
    typedef System Sys;

    void print_stepper_info() override {
        euler_mayurama_stepper::print_stepper_info();
        cout << "k = " << k << endl;
        cout << "prev_accepted_k = " << prev_accepted_k << endl;
        cout << "switched = " << switched << endl;
        cout << "switch_count = " << switch_count << endl;
    }

    void do_step(Sys& sys, state_type& x, time_type dt_max, time_type &t) override {
        if(switch_counter > switch_count * k) {
            timer.set_startpoint(euler_steps);
            if (!switched) {
                cout << "switched: k = " << k << " dt = " << dt << endl;
                cout << "or are we at prev_k? prev_k = " << prev_accepted_k << endl;
                switched = true;
            }
            euler_mayurama_stepper::do_step(sys, x, dt, t);
            t += dt;
            timer.set_endpoint(euler_steps);
        } else {
            timer.set_startpoint(adaptive_steps);
            do_adaptive_step(sys, x, dt_max, t);
            timer.set_endpoint(adaptive_steps);
        }
    }

/*
    template<class Sys>
    void step_until(time_type end_time, Sys& sys, state_type& x, time_type dt_max, time_type &t) {
        while (t < end_time){
            do_step(sys, x, dt_max, t);
        }
    }
*/

    // we now pass dt by reference, so that we can modify it
    void do_adaptive_step(Sys& sys, state_type& x, time_type dt_max, time_type &t) {

        // calc the stepsize with the current k
        dt = 1.0 / pow(reduction_factor, k) * dt_max;

        // good thing is that we already split the drift and the diffusion
        // we are at t and can easily apply the system operations without thinking of it for now
        sys.calc_drift(x, dxdt, t);
        // Here i think we have to apply only dxdt for now to calculate f(x*) - f(x)
        // which means we need a new operations structure?
        algebra::for_each(x_drift, x, dxdt, apply_drift(dt));
        sys.map_state(x_drift);
        // we have x_drift now, now we need to calculate f(x_drift)
        // we really should not just call the system since the system will generate random numbers
        // so we just call again calc drift, but we need to store the result somewhere else than dxdt since
        // dxdt= f(x)
        timer.set_startpoint(drift_calc);
        sys.calc_drift(x_drift, dx_drift_dt, t);
        timer.set_endpoint(drift_calc);
        // now we need to calculate the difference between dx_drift_dt = f(x_drift) and dxdt=f(x)
        // how do we do that? we actually for a simple case just need to calc the difference for every
        // entry of dxdt and dx_drift_dt and then sum it up / average it. this should actually be a very simple
        // thrust operation
        timer.set_startpoint(error_calc);
        algebra::for_each(dx_drift_dt, dxdt, calc_error());
        // now the error is in dx_drift_dt, now we got to sum and average it
        error = sum()(dx_drift_dt) / N;
//        cout << "Error: " << error << endl;
//        cout << "x[0]: " << x[0] << endl;
//        cout << "dxdt[0]: " << dxdt[0] << endl;
        timer.set_endpoint(error_calc);
        // now we have to check whether the error is small enough
        if(error < tol) {
            // if error is smaller than the tolerance we apply everything
            timer.set_startpoint(rng);
            // cout << "calcing diff"<< endl;
            sys.calc_diff(theta, t);
            // cout << "THETA IN ADAPTIVE STEPPER:" << endl;
            // print_container(theta);
            // cout << endl << endl;
            timer.set_endpoint(rng);
            algebra::for_each(x, x_drift, theta, apply_diff(dt));
            sys.map_state(x);
            // we also increase the time
            t += dt;
            // and we reduce k for the next step
            // do we have to check that k does not get smaller than zero?
            // wo know that it was accepted, so if our accepted k is equal to the last accepted k, we increase the counter
            if (k == prev_accepted_k) {
                ++switch_counter;
            } else {
                // if athe accepted k is not equal to the last accepted one, we set the new k and reset the counter
                // cout << "Resetting because prev_accepted k:" << prev_accepted_k << " vs k: " << k << endl;
                prev_accepted_k = k;
                switch_counter = 0;
            }

            // TODO is there a faster way to do this than with if?
            if (k > 0) {
                k--;
            }
        } else {
            // everytime we enter else we failed the error test, so if we time this, we now how long the repetitions
            // took
            timer.set_startpoint(repetitions);
            // if the error is to large, we have to do the step again with increased k
            k++;
            do_step(sys, x, dt_max, t);
            timer.set_endpoint(repetitions);
            // I thianak thats it?
        }

    }

    euler_combined() {} // TODO default constructor, could be inherited but...

    euler_combined(size_t N, int K, double tol, int switch_count = 1, double reduction_factor=1.5) : euler_mayurama_stepper(N),
                                                                                                         dx_drift_dt(N), x_drift(N), reduction_factor(reduction_factor),
                                                                                                         k(K), tol(tol), prev_accepted_k(K), switch_count(switch_count), first_k(K)
    {
        cout << "creating euler combined stepper" << endl;
    }

    euler_combined(map<Parameter, double> paras): euler_combined(2*(int)paras[total_size], paras[K], paras[Parameter::tol]){
        cout << "Initializing euler_combined stepper from paras" << endl;
        cout << "N =" << N << endl;
        switch_count = 0;
        switch_counter = 1;
    }

    int get_k() {
        return k;
    }

    double get_error() {
        return error;
    }

    void reset() {
        // resetting the stepper for another use, the switch counter at least should be set to zero again
        switch_counter = 0;
        // TODO i think dxdt, theta etc don't have to be reset as they get overwritten every step anyway
        // TODO i think we don't have to reset prev_accepted k, but the normal one
        k = first_k;
        // something else?
    }


private:
    state_type x_drift, dx_drift_dt;
};

/*template<
        class state_type,
        class algebra,
        class operations,
        class System,
        class value_type = double,
        class time_type = value_type
>
class euler_stepper : public euler_mayurama_stepper<state_type, algebra, operations, System, value_type, time_type>{
    typedef typename operations::template apply_drift<time_type> apply_drift;
    typedef typename operations::template sum<value_type> sum;
    typedef typename operations::template apply_diff<time_type> apply_diff;
    typedef euler_mayurama_stepper<state_type, algebra, operations, System, value_type, time_type> euler_mayurama_stepper;
    typedef stepper<state_type, algebra, operations, System, value_type, time_type> stepper;
    using stepper::dxdt;
    state_type F = dxdt;
    using stepper::theta;
    using stepper::N;       // this is 2 * n so the size of a normal vector
    using stepper::stepper;
    using euler_mayurama_stepper::theta;
    using euler_mayurama_stepper::N;
    using euler_mayurama_stepper::dxdt;
*//*    using do_euler_step = euler_mayurama_stepper<state_type, algebra, operations,
            value_type, time_type>::do_step;*//*
public:
    typedef System Sys;

    void print_stepper_info() override {
        euler_mayurama_stepper::print_stepper_info();
    }

    void do_step(Sys& sys, state_type& x, time_type dt, time_type &t) override {
        euler_mayurama_stepper::do_step(sys, x, dt, t);
        t += dt;
    }

*//*
    template<class Sys>
    void step_until(time_type end_time, Sys& sys, state_type& x, time_type dt_max, time_type &t) {
        while (t < end_time){
            do_step(sys, x, dt_max, t);
        }
    }
*//*


    euler_stepper(size_t N) : euler_mayurama_stepper(N)
    {
    }
    euler_stepper(map<Parameter, double> paras): euler_stepper(2*(int)paras[total_size]){
        cout << "initializing euler stepper from paras" << endl;
    }

private:
    state_type x_drift, dx_drift_dt;
};*/

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
        Singleton_timer::set_startpoint(bbk_checkpoint_v2);
        //                     v(n+1/2)           F(q_n+1)               R_n+1
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

template<
        class state_type,
        class algebra,
        class operations,
        class System,
        class value_type = double,
        class time_type = value_type
>
class bbk_stepper2 : public stepper<state_type, algebra, operations, System, value_type, time_type> {
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

    void do_step(Sys& sys, state_type& x, time_type dt, time_type& t) override {
        // We somehow need to wave in this p~ into our architecture
        // and what is fundamentally different is that we need two random variables per side per step
        // and we actually cannot really multiply it with the temperature beforehand which might be a problem...
        // seems like we sadly have to restructure a bit? New methods for bbk aswell as splitting old stuff to be more
        // modular

        // I am a bit unsure of how i am supposed to implement this shhhiii... I think dx_dt will now be used for F(x(t))
        // the p~ will be saved in the p spots of x?

        // I somehow need to calculate v~(t) then x(t +dt) then v(t +dt) sequientially? I cannot do it simoultanously anymore?
        // since i cannot do it simoultaneously, there is no reason to have q and p in the same vector anymore?
        // we should probably still leave it like that because we otherwise get compatibility problems everywhere?
        // weird thing is just that we now have the first half of dxdt always empty so pretty useless
        // I guess it was the same deal with theta before...
        // okay so we use our first half of x to calculate the force which we save in the second half of F
        sys.calc_force(x, F, t);
        // now we have to apply it, i think we can use apply_drift for this
        //                      v~          v                   F
        algebra::for_each(x.begin() + n, x.begin() + n, F.begin() + n, n, apply_drift(0.5 * dt));
        // now this should have done x = x and p = p + dt/ 2 * F
        // now we need to calculate q (t + dt), therefore we need the theta vector
        sys.calc_diff_bbk(theta, dt);
        // The theta vector should now be filled with the temp * zeta values
        // now to apply this we definetly need a new integrater, this one it even has to know the dampening
        // easiest would probably be to get the dampening out of the system.
        // Watch out! you need zeta_2 first but you wrot it into second half of theta.
        //                  q           v~              zeta_2
        algebra::for_each(x.begin(), x.begin() + n, theta.begin() + n, n, apply_bbk_q(dt, sys.get_eta()));
        sys.map_state(x);
        // now we can calculate the force again
        sys.calc_force(x, F, t);
        // now we need to integrate p(t + dt) which needs again its own integrator
        //                     v~           F               zeta_1
        algebra::for_each(x.begin() + n, F.begin() + n, theta.begin(), n, apply_bbk_p(dt, sys.get_eta()));
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
    bbk_stepper2(size_t N) : stepper(N), n(N/2)
    {
    }
    bbk_stepper2(map<Parameter, double> paras): bbk_stepper2(2*(int)paras[total_size]){}


    void print_stepper_info() override {
        stepper::print_stepper_info();
    }
};

template<
        class state_type,
        class algebra,
        class operations,
        class System,
        class value_type = double,
        class time_type = value_type,
        template<class, class, class, class, class, class> class stepper_type
>
stepper_type<state_type, algebra, operations, System, value_type, time_type>* create_stepper(map<Parameter, double>& paras, int n) {
    int N = 2;   // N for now no parameter
    cout << stepper_type<state_type, algebra, operations, System, value_type, time_type>::get_name();
    if (stepper_type<state_type, algebra, operations, System, value_type, time_type>::get_name() == "bbk stepper") {
        return new stepper_type<state_type, algebra, operations, System, value_type, time_type>(N * n);
    }
    else if(is_same<stepper_type<state_type, algebra, operations, System, value_type, time_type>,
            euler_combined<state_type, algebra, operations, System, value_type, time_type>>::value){
        // damn tahts ugly#
        // TODO This cannot! possibly work, i see it already. probably we have to exchange one 'class' with template<size_t> class
        double K = paras[Parameter::K];
        double tol = paras[Parameter::tol];
        // return new stepper_type<state_type, algebra, operations, System, value_type, time_type>(N * n, K ,tol);
        return new stepper_type<state_type, algebra, operations, System, value_type, time_type>(N * n);
    } else {
        return new stepper_type<state_type, algebra, operations, System, value_type, time_type>(N * n);

    }// TODO other steppers
}

template<
        class state_type,
        class algebra,
        class operations,
        class System,
        class value_type = double,
        class time_type = value_type,
        template<class, class, class, class, class, class> class stepper_type
>
stepper_type<state_type, algebra, operations, System, value_type, time_type>* create_stepper(map<Parameter, double>& paras) {
    int n =(int) (paras[dim_size_x] * paras[dim_size_y]);
    return create_stepper<state_type, algebra, operations, System, value_type, time_type, stepper_type>(paras, n);
}


#endif //CUDAPROJECT_STEPPERS_CUH
