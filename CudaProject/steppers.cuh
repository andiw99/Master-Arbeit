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

    stepper() {}
    stepper(size_t N) : N(N), dxdt(N), theta(N) {}
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
            obs->operator()(sys, x, t); // Observing before anything happens
        }
        while (t < end_time){
            this->do_step(sys, x, dt_max, t); // it is important that the steper custom do step is called here
            for(auto obs : obsvers) {
                obs->operator()(sys, x, t); // Observing before anything happens
            }
        }
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
    string system_name = "System";
    string applying_name = "Applying Euler";
    string rng = "RNG during Euler";
    string functor = "Functor during Euler";
    checkpoint_timer timer{{system_name, applying_name, rng, functor}};
public:
    typedef System Sys;
    using stepper<state_type, algebra, operations, System, value_type, time_type>::dxdt;
    using stepper<state_type, algebra, operations, System, value_type, time_type>::theta;
    using stepper<state_type, algebra, operations, System, value_type, time_type>::N;
    using stepper<state_type, algebra, operations, System, value_type, time_type>::stepper;
    // observer* Observer;
    // the stepper needs a do_step method
    // I think our do step method needs an additional parameter for theta? Maybe not, we will see
    // We also template our do_step method to work with any system that we feed into it
    // i think later the system is called to perform the operation on the state types


    void do_step(Sys& sys, state_type& x, time_type dt, time_type t) {
        // okay we don't need this for_each3 stuff, we only need to apply x_(n+1) = x_n + k_1 dt + theta sqrt(dt)
        // first we need to calculate the derivative with the system and save it to temporary class things
        // i don't even think that we need dxdt as parameter

        // We call the system to calculate dxdt, should only calculate the deterministic part here?
        // No i would say we also calculate the stochastic part
        // so we give the system the state and istruct it later to save its calculations in dxdt and theta so that we
        // can later iterate over the lattice and apply the operation, meaning the update
        timer.set_startpoint(system_name);
        timer.set_startpoint(functor);
        sys.calc_drift(x, dxdt, t);
        timer.set_endpoint(functor);
        timer.set_startpoint(rng);
        sys.calc_diff(theta, t);
        timer.set_endpoint(rng);
        timer.set_endpoint(system_name);
        // this should set the correct values for dxdt and theta so that they can be applied in apply_em
        // can we print them here?

//        cout << "x[0] = " << x[0] << ", " << "x[1] = " << x[1] << endl;
//        cout << "dxdt[0] = " << dxdt[0] << ", " << "dxdt[1] = " << dxdt[1] << endl;
//        cout << "theta[0] = " << theta[0] << ", " << "theta[1] = " << theta[1] << endl;

        // for the update we need to define a type of a function, but i don't really understand the syntax
        // okay so the plan is to have a new name for the templated struct apply_em that is a member of operations
        // so that we can later more clearly instantiate the struct apply_em with the correct template
        timer.set_startpoint(applying_name);
        typedef typename operations::template apply_em<time_type> apply_em;
        // this is the operation that i want to apply on every lattice site
        algebra::for_each(x, x, dxdt, theta, apply_em(dt));
        timer.set_endpoint(applying_name);
        // and that should already be it?
        // Observe? How much time does this take?
        // so i also want to write down the temperature of the step, but i don't want to give it to the observer
        // since not every System has a temperature or at least one that is relevant
        // we could give the system to the observer, this would be possible, so that we can use sys.getT in a
        // special observer for the systems with temperature
        // Observer->operator()(sys, x, t);
    }

    // i am currently not sure what parameters we additionally need, we don't have temporary x values like for the
    // runge kutta scheme, at least the system size should be a parameter i guess
    // I don't really get how this stuff is instantiated here
    euler_mayurama_stepper(size_t N) : stepper<state_type, algebra, operations, System, value_type, time_type>(N) //, Observer(Obs)
    {
    }

/*    euler_mayurama_stepper(size_t N) : N(N) {
        Observer = new observer();
    }*/

    // now for the memory allocation. I am not sure if this ever changes for me but the implementation should not harm
    // on the other hand, my class doesn't have any saved state_types that could potentially be resized
    // so we skip this for now

    void print_stepper_info() override {
        stepper<state_type, algebra, operations, System, value_type, time_type>::print_stepper_info();
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
    checkpoint_timer timer{{drift_calc, error_calc, repetitions}};
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
    checkpoint_timer timer{{drift_calc, error_calc, repetitions, euler_steps, adaptive_steps, rng}};
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
            cout << "THETA IN ADAPTIVE STEPPER:" << endl;
            print_container(theta);
            cout << endl << endl;
            timer.set_endpoint(rng);
            algebra::for_each(x, x_drift, theta, apply_diff(dt));
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

    euler_combined(size_t N, int K, double tol, int switch_count = 10000, double reduction_factor=1.5) : euler_mayurama_stepper(N),
                                                                                                         dx_drift_dt(N), x_drift(N), reduction_factor(reduction_factor),
                                                                                                         k(K), tol(tol), prev_accepted_k(K), switch_count(switch_count), first_k(K)
    {
        cout << "creating euler combined stepper" << endl;
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

template<
        class state_type,
        class algebra,
        class operations,
        class System,
        class value_type = double,
        class time_type = value_type,
        template<class, class, class, class, class, class> class stepper_type
>
stepper_type<state_type, algebra, operations, System, value_type, time_type>* create_stepper(map<string, double>& paras, int n) {
    if(is_same<stepper_type<state_type, algebra, operations, System, value_type, time_type>,
            euler_combined<state_type, algebra, operations, System, value_type, time_type>>::value){
        // damn tahts ugly#
        // TODO This cannot! possibly work, i see it already. probably we have to exchange one 'class' with template<size_t> class
        int N = (int) paras["N"];
        double K = paras["K"];
        double tol = paras["tol"];
        return new stepper_type<state_type, algebra, operations, System, value_type, time_type>(N * n, K ,tol);
    } // TODO other steppers
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
stepper_type<state_type, algebra, operations, System, value_type, time_type>* create_stepper(map<string, double>& paras) {
    int n =(int) (paras["lat_dim"] * paras["lat_dim"]);
    return create_stepper<state_type, algebra, operations, System, value_type, time_type, stepper_type>(paras, n);
}


#endif //CUDAPROJECT_STEPPERS_CUH
