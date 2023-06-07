//
// Created by andi on 30.04.23.
//

#ifndef CUDAPROJECT_MAIN_CUH
#define CUDAPROJECT_MAIN_CUH
//

#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <functional>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <boost/typeof/typeof.hpp>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/random.h>
#include <thrust/reduce.h>
#include <thrust/sequence.h>
#include <thrust/execution_policy.h>
#include <fstream>
#include <filesystem>
#include <map>


// timing stuff from stackoverflow
#include <thrust/host_vector.h>
#include <thrust/generate.h>
#include <thrust/for_each.h>
#include <thrust/execution_policy.h>
#include <thrust/random.h>


using namespace std;
/*
 * So instead of the complicated runge kutta scheme we start with the euler mayurama
 * i hopoe this gpu stuff works with SDEs, i am worried about the generation of the random numbers.
 * instead of x_(n+1) = x_n + 1/6 * dt * (k_1 + 2_k_2 + 2k_3 + k_4) we
 * only need  x_(n+1) = x_n + k_1 dt + theta sqrt(dt)
 * I think i will try to achieve this flexible implementation of the paper
 * what do we need?
 * - Representation: needs to be an template parameter since we want to swap the containers from cpu to gpu and vice versa
 * - We need memeory management, meaning the allocation of the state vectors etc, but i think cuda and thrust will handle
 *      this on its own
 * - We need the vector iteration which is not entirely trivial for execution on gpu's since we have to iterate
 *      simultaneously over all lattice sites
 * - We need a gpu implementation of the operation that takes place on the lattice sites.
 */


/*
 * Now lets think about the algebra for a 'normal' implementation first
 * Am I mistaken or do I need to make sure that I iterate over every lattice site in the algebra?
 * But i thought the algebra that was given in the paper would work for arbitrary state types, but i cannot imagine
 * how this would work here. maybe i would have to typedef or define a custom state_type that can be initialized
 * with just N and can be iterated through with only using container[i] instead of container[i, j]
 * Or i am overthinking it and i just specialize for 2D
 *
 */
struct container_algebra {
    template<class S1, class S2, class S3, class S4, class Op>
    static void for_each(S1 &s1, S2 &s2, S3 &s3, S4 &s4, Op op) {
        const size_t dim = s1.size();
        for(size_t n = 0; n < dim; ++n)
            op(s1[n], s2[n], s3[n], s4[n]);
    }
};

/*
 * Now we need a "Thrust Algebra", will this be doable? I still need to figure out what my statetype is and how I
 * iterate over the lattice. I hope in the coupled oscillator example this will get clear
 */
struct thrust_algebra {
    template<class S1, class S2, class S3, class S4, class Op>
    static void for_each(S1 &s1, S2 &s2, S3 &s3, S4 &s4, Op op) {
        thrust::for_each(
                thrust::make_zip_iterator( thrust::make_tuple(
                        s1.begin(), s2.begin(), s3.begin(), s4.begin()) ),
                thrust::make_zip_iterator( thrust::make_tuple(
                        s1.end(), s2.end(), s3.end(), s4.begin())),
                op);
    }
    template<class S1, class S2, class S3, class Op>
    static void for_each(S1 &s1, S2 &s2, S3 &s3, Op op) {
        thrust::for_each(
                thrust::make_zip_iterator( thrust::make_tuple(
                        s1.begin(), s2.begin(), s3.begin()) ),
                thrust::make_zip_iterator( thrust::make_tuple(
                        s1.end(), s2.end(), s3.end())),
                op);
    }
    template<class S1, class S2, class Op>
    static void for_each(S1 &s1, S2 &s2, Op op) {
        thrust::for_each(
                thrust::make_zip_iterator( thrust::make_tuple(
                        s1.begin(), s2.begin()) ),
                thrust::make_zip_iterator( thrust::make_tuple(
                        s1.end(), s2.end())),
                op);
    }
};

/*
 * Next up is the operations, but this should be simple, i shuld only need one operation that applies the em scheme
 */
struct default_operations {
    template<class time_type = double>
    struct apply_em {
        const time_type dt;
        apply_em(time_type dt)
                : dt(dt) { }
        template<class T0, class T1, class T2>
        void operator()(T0 &x, const T0 &x0, const T1 &dxdt, const T2 &theta) const {
            // we somehow need to make sure that those vector operations work?
            // i dont see how this would apply to my 2D lattice with (x, p) on every site atm
            x = x0 + dt * dxdt + sqrt(dt) * theta;
        }
    };
};

/*
 * Now thrust operations that have to apply em with only one factor but we have to make a sqrt transform somewhere
 */

struct thrust_operations {
    template<class time_type = double>
    struct apply_em {
        const time_type dt;
        apply_em(time_type dt)
                : dt(dt) { }
        template< class Tuple >
        __host__ __device__ void operator()(Tuple tup) const {
            // The most simple tuple would be for only one lattice site, then the tuple would be
            // tup = ( x = (q, p), dxdt = (dqdt, dpdt), theta = (theta_q = 0, theta_p = rand))
            // and since for_each does the following code for every entry of every vector in the tuple
            thrust::get<0>(tup) = thrust::get<1>(tup) +
                                  dt * thrust::get<2>(tup) +
                                  sqrt(dt) * thrust::get<3>(tup);
            // the above code translates into
            // q = q + dt * dqdt + theta_q * sqrt(dt)
            // p = p + dt * dpdt + theta_p * sqrt(dt)
            // if we would have more lattice sites that are ordered correctly, the next equations would be for
            // q2 = ...
            // p2 = ...
        }
    };

    template<class time_type = double>
    struct apply_drift {
        const time_type dt;
        apply_drift(time_type dt)
                : dt(dt) { }
        template< class Tuple >
        __host__ __device__ void operator()(Tuple tup) const {
            thrust::get<0>(tup) = thrust::get<1>(tup) +
                                  dt * thrust::get<2>(tup);

        }
    };

    struct calc_error {
        calc_error() {}
        template< class Tuple >
        __host__ __device__ void operator()(Tuple tup) const {
            // here we only have f(x_drift) and f(x)
            //TODO  we just save the difference in the dx_drift_dt vector to save time? is this a good idea?
            thrust::get<0>(tup) = abs(thrust::get<0>(tup) - thrust::get<1>(tup));
        }
    };

    // TODO somehow i have to specify thrust::host here since i am to dumb to get it working withoug
    template <class valuetype = double>
    struct sum {
        sum() {}
        template <class statetype>
        __host__ __device__ valuetype operator()(statetype &error)  {
            return thrust::reduce(error.begin(), error.end(),  (valuetype) 0.0, thrust::plus<valuetype>());
        }
    };

    template<class time_type = double>
    struct apply_diff {
        const time_type dt;
        apply_diff(time_type dt)
                : dt(dt) { }
        template< class Tuple >
        __host__ __device__ void operator()(Tuple tup) const {
            thrust::get<0>(tup) = thrust::get<1>(tup) +
                                  sqrt(dt) * thrust::get<2>(tup);

        }
    };
};





// We need those observers that do stuff like writing the values to a file
class observer {
public:
    template<class State, class System>
    void operator()(const System sys, const State &x , double t ) {
        cout << "Base Observer gets called" << endl;
    }

    observer() {
    }

};

// Observer for lattic on bath specific for gpu_bath system?
class bath_observer : public observer {
    // this one has a filestream
    ofstream& file;
    // and a write_every
    const size_t write_every;
    // and a count? Is that all efficient?
    size_t count;

public:
    bath_observer(ofstream& out, size_t write_every = 1000) : file(out), write_every(write_every) {
        count = 1;
    }

    template<class State, class System>
    void operator()(System &sys, const State &x , double t ) {
        // so i think if statements are not to costly on c++
        if (count % write_every == 0 || count == 1) {
            // now we write the state to the file
            // i call a function here that would not be available for all Systems, not to clean and a bit fishy
            // but i guess it works for now
            size_t lat_dim = sys.get_lattice_dim();
            double T = sys.get_cur_T();

            // writing, but only q for the moment
            file << "t, " << t << ",";
            for(int i = 0; i < lat_dim * lat_dim; i++) {
                file << x[i] << ",";
                // cout << x[i] << endl;

            }
            // for the last we write T? I mean thats how i did it, right?
            // Temperatur saven
            file << T;

            // Zeilenumbruch
            file << "\n";
        }
        count++;
    }
    template<class State, class System>
    void write(System &sys, const State &x , double t ) {
        size_t lat_dim = sys.get_lattice_dim();
        double T = sys.get_cur_T();

        file << "t, " << t << ",";
        for(int i = 0; i < lat_dim * lat_dim; i++) {
            file << x[i] << ",";
            // cout << x[i] << endl;

        }
        // for the last we write T? I mean thats how i did it, right?
        // Temperatur saven
        file << T;

        // Zeilenumbruch
        file << "\n";
    }
};

// Observer for lattic on bath specific for gpu_bath system?
class chain_observer : public observer {
    // this one has a filestream
    ofstream& file;
    // and a write_every
    const size_t write_every;
    // and a count? Is that all efficient?
    size_t count;

public:
    chain_observer(ofstream& out, size_t write_every = 1000) : file(out), write_every(write_every) {
        count = 1;
    }

    template<class State, class System>
    void operator()(const System sys, const State &x , double t ) {
        // so i think if statements are not to costly on c++
        if (count % write_every == 0 || count == 1) {
            // now we write the state to the file
            // i call a function here that would not be available for all Systems, not to clean and a bit fishy
            // but i guess it works for now
            size_t lat_dim = sys.get_lattice_dim();
            double T = sys.get_cur_T();

            // writing, but only q for the moment
            file << "t, " << t << ",";
            for(int i = 0; i < lat_dim; i++) {
                file << x[i] << ",";
                // cout << x[i] << endl;

            }
            // for the last we write T? I mean thats how i did it, right?
            // Temperatur saven
            file << T;

            // Zeilenumbruch
            file << "\n";
        }
        count++;
    }

    template<class State, class System>
    void write(const System sys, const State &x , double t ) {
        size_t lat_dim = sys.get_lattice_dim();
        double T = sys.get_cur_T();

        file << "t, " << t << ",";
        for(int i = 0; i < lat_dim; i++) {
            file << x[i] << ",";
            // cout << x[i] << endl;

        }
        // for the last we write T? I mean thats how i did it, right?
        // Temperatur saven
        file << T;

        // Zeilenumbruch
        file << "\n";
    }

};

struct timer {
    chrono::time_point<chrono::high_resolution_clock> starttime;
public:
    timer() {
        starttime = chrono::high_resolution_clock::now();
    }
    ~timer() {
        auto endtime = chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                endtime - starttime);
        auto duration_count = duration.count();
        cout << "total execution took " << duration_count << "ms" << endl;
    }
};
// TODO this is basically the same as the system timer, we could implement a "Checkpoint timer" that would be universal
// to use
struct checkpoint_timer : public timer {
    map<string, long> times;
    map<string, chrono::time_point<chrono::high_resolution_clock>> startpoints;
    map<string, chrono::time_point<chrono::high_resolution_clock>> endpoints;
    int nr_checkpoints;
public:
    // we now make a function, that sets the startpoint of checkpoint i
    void set_startpoint(const string& checkpoint) {
        startpoints[checkpoint] = chrono::high_resolution_clock::now();
    }

    void set_endpoint(const string& checkpoint) {
        endpoints[checkpoint] = chrono::high_resolution_clock::now();
        // add the duration
        times[checkpoint] += std::chrono::duration_cast<std::chrono::microseconds>(
                endpoints[checkpoint] - startpoints[checkpoint]).count();
    }

    checkpoint_timer(const vector<string>& checkpoint_names) : timer() {
        // constructor takes a list of names of the checkpoints
        // TODO we could use a map that maps the names to the start and endpoints? Or is that slower than working
        // with indices?
        // initialize the maps
        for(const string& name : checkpoint_names) {
            startpoints[name] = chrono::high_resolution_clock::now();
            endpoints[name] = chrono::high_resolution_clock::now();
            times[name] = 0;
        }
    }

    ~checkpoint_timer() {
        // print the durations
        for(const auto& timepair : times) {
            cout << timepair.first << " took " << (double)timepair.second * 0.001 << " ms" << endl;
        }
    }

};

/*
 * We start with the stepper, this stepper has to be templated since we want to exchange the container for our state type
 * We also template the 'Algebra' that outsources the iteration through my lattice sites?
 * And the 'Operation' that performs the operation on the lattice sites
 */


template<
        class state_type,
        class algebra,
        class operations,
        class value_type = double,
        class time_type = value_type
>
class euler_mayurama_stepper{
    string system_name = "System";
    string applying_name = "Applying Euler";
    checkpoint_timer timer{{system_name, applying_name}};
public:
    // observer* Observer;
    // the stepper needs a do_step method
    // I think our do step method needs an additional parameter for theta? Maybe not, we will see
    // We also template our do_step method to work with any system that we feed into it
    // i think later the system is called to perform the operation on the state types


    template<class Sys>
    void do_step(Sys& sys, state_type& x, time_type dt, time_type t) {
        // okay we don't need this for_each3 stuff, we only need to apply x_(n+1) = x_n + k_1 dt + theta sqrt(dt)
        // first we need to calculate the derivative with the system and save it to temporary class things
        // i don't even think that we need dxdt as parameter

        // We call the system to calculate dxdt, should only calculate the deterministic part here?
        // No i would say we also calculate the stochastic part
        // so we give the system the state and istruct it later to save its calculations in dxdt and theta so that we
        // can later iterate over the lattice and apply the operation, meaning the update
        timer.set_startpoint(system_name);
        sys.calc_drift(x, dxdt, t);
        sys.calc_diff(theta, t);
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
    euler_mayurama_stepper(size_t N) : N(N), dxdt(N), theta(N) //, Observer(Obs)
    {
    }

/*    euler_mayurama_stepper(size_t N) : N(N) {
        Observer = new observer();
    }*/

    // now for the memory allocation. I am not sure if this ever changes for me but the implementation should not harm
    // on the other hand, my class doesn't have any saved state_types that could potentially be resized
    // so we skip this for now
protected:
    // the system size, but what is N if i am in 2 dimensions? probably n * n. But how can i initialize the inital
    // values for all sites? not just by "stat_type(N)"
    // i still don't fully understand what the state_type is, is it just (x, p) for one lattice site? Or is it the
    // the state of the whole system? In the former, N would correspond with the dimensionality of my DGL system
    // in the latter i would not really know how to initialize the initial states independent of the dimensionality
    // of my problem
    const size_t N;
    state_type dxdt, theta;
};

template<
        class state_type,
        class algebra,
        class operations,
        class value_type = double,
        class time_type = value_type
>
class euler_simple_adaptive{
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

    // we now pass dt by reference, so that we can modify it
    template<class Sys>
    void do_step(Sys& sys, state_type& x, time_type dt_max, time_type &t) {

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

    euler_simple_adaptive(size_t N, int K, double tol) : N(N), dxdt(N), dx_drift_dt(N), x_drift(N), theta(N), k(K), tol(tol)
    {}

    int get_k() {
        return k;
    }

    double get_error() {
        return error;
    }

private:
    const size_t N;
    state_type x_drift, dxdt, dx_drift_dt, theta;
};

template<
        class state_type,
        class algebra,
        class operations,
        class value_type = double,
        class time_type = value_type
>
class euler_combined : public euler_mayurama_stepper<state_type, algebra, operations, value_type, time_type>{
    int k;
    int prev_accepted_k;
    value_type tol;
    value_type error = 0;
    time_type dt;
    int switch_counter;
    int switch_count;
    double reduction_factor;
    bool switched = false;
    typedef typename operations::template apply_drift<time_type> apply_drift;
    typedef typename operations::calc_error calc_error;
    typedef typename operations::template sum<value_type> sum;
    typedef typename operations::template apply_diff<time_type> apply_diff;
    using euler_mayurama_stepper<state_type, algebra, operations,
            value_type, time_type>::theta;
    using euler_mayurama_stepper<state_type, algebra, operations,
            value_type, time_type>::N;
    using euler_mayurama_stepper<state_type, algebra, operations,
            value_type, time_type>::dxdt;
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

    template<class Sys>
    void do_step(Sys& sys, state_type& x, time_type dt_max, time_type &t) {
        if(switch_counter > switch_count * k) {
            timer.set_startpoint(euler_steps);
            if (!switched) {
                cout << "switched: k = " << k << " dt = " << dt << endl;
                switched = true;
            }
            euler_mayurama_stepper<state_type, algebra, operations,
                    value_type, time_type>::do_step(sys, x, dt, t);
            t += dt;
            timer.set_endpoint(euler_steps);
        } else {
            timer.set_startpoint(adaptive_steps);
            do_adaptive_step(sys, x, dt_max, t);
            timer.set_endpoint(adaptive_steps);
        }
    }

    // we now pass dt by reference, so that we can modify it
    template<class Sys>
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
        timer.set_endpoint(error_calc);
        // now we have to check whether the error is small enough
        if(error < tol) {
            // if error is smaller than the tolerance we apply everything
            timer.set_startpoint(rng);
            sys.calc_diff(theta, t);
            timer.set_endpoint(rng);
            algebra::for_each(x, x_drift, theta, apply_diff(dt));
            // we also increase the time
            t += dt;
            // and we reduce k for the next step
            // do we have to check that k does not get smaller than zero?
            // wo know that it was accepted, so if our accepted k is equal to the last accepted k, we increase the counter
            if (k == prev_accepted_k) {
                ++switch_count;
            } else {
                // if athe accepted k is not equal to the last accepted one, we set the new k and reset the counter
                prev_accepted_k = k;
                switch_count = 0;
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

    euler_combined(size_t N, int K, double tol, int switch_count = 100, double reduction_factor=1.5) : euler_mayurama_stepper<state_type, algebra, operations, value_type, time_type>(N),
             dx_drift_dt(N), x_drift(N), reduction_factor(reduction_factor),
    k(K), tol(tol), prev_accepted_k(K), switch_count(switch_count)
    {}

    int get_k() {
        return k;
    }

    double get_error() {
        return error;
    }

private:
    state_type x_drift, dx_drift_dt;
};

// I don't know if this is so smart what i did here since i don't know whether this grid searching will work with
// this templated stuff. I mean we could probably write a shell script that compiles and executes the code for the
// different values but I doubt that it will be worth it
string trunc_double(double a, int precision=2) {
    stringstream stream;
    stream << std::fixed << std::setprecision(precision) << a;
    return stream.str();
}

void create_dir(const string dir_name) {
    // check wheter the directory already exists, if not create it
    if(!filesystem::is_directory(dir_name) || !filesystem::exists(dir_name)) {
        filesystem::create_directories(dir_name);
    }
}

string create_tree_name(double eta, double T, double dt, int n, double alpha, double beta, double J, double tau,
                        const string root) {
    string dir_name = root + "eta=" + trunc_double(eta)
                      + "/T=" + trunc_double(T) + "/dt=" +
                      trunc_double(dt, 4) + "/n="+ to_string(n) + "/alpha=" + trunc_double(alpha) + "/beta=" +
                      trunc_double(beta) + "/J=" + trunc_double(J) + "/tau=" + trunc_double(tau);
    create_dir(dir_name);
    return dir_name;
}

void write_parameters(ofstream& file, double eta, double T, double dt, int n, double alpha, double beta, double J,
                      double tau) {
    // insert the parameters
    file << "eta," << eta << ", \n";
    file << "T," << T << ", \n";
    file << "dt," << dt << ", \n";
    file << "n," << n << ", \n";
    file << "alpha," << alpha << ", \n";
    file << "beta," << beta << ", \n";
    file << "J," << J << ", \n";
    file << "tau," << tau << ", \n";
}

template<typename T>
std::vector<T> linspace(T start_in, T end_in, int num_in)
{

    std::vector<T> linspaced;

    double start = static_cast<double>(start_in);
    double end = static_cast<double>(end_in);
    double num = static_cast<double>(num_in);

    if (num == 0) { return linspaced; }
    if (num == 1)
    {
        linspaced.push_back(start);
        return linspaced;
    }

    double delta = (end - start) / (num - 1);

    for(int i=0; i < num-1; ++i)
    {
        linspaced.push_back((T)(start + delta * i));
    }
    linspaced.push_back(end); // I want to ensure that start and end
    // are exactly the same as the input
    return linspaced;
}
void print_vector(std::vector<double> vec)
{
    std::cout << "size: " << vec.size() << std::endl;
    for (double d : vec)
        std::cout << d << " ";
    std::cout << std::endl;
}

double testvariable = 2.0;

struct rand_init_values
{
    double mu, sigma, ampl;

    __host__ __device__
    rand_init_values(double ampl, double mu = 0.0, double sigma = 1.0) : ampl(ampl), mu(mu), sigma(sigma) {};

    __host__ __device__
    float operator()(const unsigned int ind) const
    {
        thrust::default_random_engine rng;
        thrust::normal_distribution<double> dist(mu, sigma);
        rng.discard(ind);

        // dist * ampl zu returnen ist wie... aus dist mit std ampl zu ziehen: b * N(m, o) = N(m, b*o)
        auto rnr = (double)dist(rng);
        // if the drawn number is greater than 3 * sigma we cut it
        if(abs(rnr) > 3.0 * sigma) {
            rnr = (rnr < 0) ? (-1.0) * 3.0 * sigma : 3.0 * sigma;
        }
        return (double)ampl * rnr;
    }
};

template <size_t n>
void fill_init_values(thrust::device_vector<double>& state, float x0, float p0, int run = 0, double mu=0, double sigma=1) {

    thrust::counting_iterator<size_t> index_sequence_begin(run * state.size());
    // thrust::fill(theta.begin(), theta.begin() + n, 0);
    // n is system size
    // fill the starting positions
    thrust::transform(index_sequence_begin,
                      index_sequence_begin + n,
                      state.begin(),
                      rand_init_values(x0, mu, sigma));
    // fill starting impulses
    thrust::transform(index_sequence_begin + n,
                      index_sequence_begin + 2*n,
                      state.begin() + n,
                      rand_init_values(p0, mu, sigma));
}


#endif //CUDAPROJECT_MAIN_CUH
