//
// Created by andi on 30.04.23.
//

#ifndef CUDAPROJECT_MAIN_CUH
#define CUDAPROJECT_MAIN_CUH
//

#include "../LearningProject/Header/Helpfunctions and Classes.h"
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
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <unsupported/Eigen/NonLinearOptimization>


// timing stuff from stackoverflow
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

template<typename scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
struct Functor
{
    typedef scalar Scalar;
    enum {
        InputsAtCompileTime = NX,
        ValuesAtCompileTime = NY
    };
    typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
    typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;

    int m_inputs, m_values;

    Functor() : m_inputs(InputsAtCompileTime), m_values(ValuesAtCompileTime) {}
    Functor(int inputs, int values) : m_inputs(inputs), m_values(values) {}

    int inputs() const { return m_inputs; }
    int values() const { return m_values; }

};

struct ExponentialDecayFunctor : Functor<double> {
private:
    // Observations for a sample.
    Eigen::MatrixXd C_x;
public:
    ExponentialDecayFunctor(const Eigen::MatrixXd &C_x): Functor(), C_x(C_x) {}

    int operator()(const Eigen::VectorXd &paras, Eigen::VectorXd &fvec) const {
        // fvec is the vector with residuals
        for(int i = 0; i < this->values(); i++) {
            fvec(i) = C_x(i, 1) - paras(0) * exp(- C_x(i, 0) / paras(1));
        }
        return 0;
    }
};

struct NumericalExpDecayDiffFunctor : Eigen::NumericalDiff<ExponentialDecayFunctor> {
public:
    NumericalExpDecayDiffFunctor(const Eigen::MatrixXd &C_x) : Eigen::NumericalDiff<ExponentialDecayFunctor>(C_x)
            {}
};

class Fitter {
private:
    int kNumObservations;
public:
    Fitter() {

    }

    template <class Residual>
    vector<double> fit() {

    }
};

template <template<class valuetype> class vec, class value_type>
vec<value_type> average_two(const vec<value_type> &a, const vec<value_type> &b) {
    // vec needs to be initializable like a vector
    vec<value_type> result(a.size());

    // averaging
    for(int i = 0; i < a.size(); i++) {
        result[i] = (a[i] + b[i]) / (value_type) 2.0;
    }

    return result;
}

template <class vec>
vec average_two(const vec &a, const vec &b) {
    // vec needs to be initializable like a vector
    vec result(a.size());

    // averaging
    for(int i = 0; i < a.size(); i++) {
        result[i] = (a[i] + b[i]) / 2.0;
    }

    return result;
}

template <template<typename valuetype> typename vec, typename value_type>
Eigen::MatrixXd construct_matrix(const vec<value_type> &a, const vec<value_type> &b) {
    Eigen::VectorXd a_eigen = Eigen::Map<Eigen::VectorXd>(a.data(), a.size());
    Eigen::VectorXd b_eigen = Eigen::Map<Eigen::VectorXd>(b.data(), b.size());

    Eigen::MatrixXd result_matrix(a.size(), 2);
    result_matrix << a_eigen, b_eigen;

    return result_matrix;
}

template <typename vec>
Eigen::MatrixXd construct_matrix(const vec &a, const vec &b) {
    Eigen::VectorXd a_eigen = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(a.data(), a.size());
    Eigen::VectorXd b_eigen = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(b.data(), b.size());

    Eigen::MatrixXd result_matrix(a.size(), 2);
    result_matrix << a_eigen, b_eigen;

    return result_matrix;
}

Eigen::MatrixXd construct_matrix(vector<double> &a, vector<double> &b) {
    Eigen::VectorXd a_eigen = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(a.data(), (long)a.size());
    Eigen::VectorXd b_eigen = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(b.data(), (long)b.size());

    Eigen::MatrixXd result_matrix(a.size(), 2);
    result_matrix << a_eigen, b_eigen;

    return result_matrix;
}

Eigen::MatrixXd construct_matrix(Eigen::VectorXd &a, Eigen::VectorXd &b) {

    Eigen::MatrixXd result_matrix(a.size(), 2);
    result_matrix << a, b;

    return result_matrix;
}


class quench_observer : public bath_observer {
    ofstream &xi_out;
public:
    quench_observer(ofstream& sys_out, ofstream& xi_out) : bath_observer(sys_out), xi_out(xi_out){
        // write the header
        xi_out << "t, xi" << endl;
    }

    template<class State, class System>
    void write_xi(System &sys, const State &x, double t) {
        // we need to calculate the correlation function and the correlation length
        // correlation function is easy but length is hard, need to do fit
        // the thing is, x is a gpu vector and to copy it to host to use my old corr_func might be slow
        // rewriting in thrust will take some time tho...
        // probably for now best to use slow variant since we probably switch to fftw sometime anyway?
        // so we use thrust copy to do this
        size_t lat_dim = sys.get_lattice_dim();
        size_t n = lat_dim * lat_dim;
        size_t nr_dists = lat_dim / 2 + 1;
        vector<double> vals(n);         // init standard vector with size of x_values which is x.size()/2

        // copy from x to f
        thrust::copy(x.begin(), x.begin() + n, vals.begin());      // okay, should be copied, we have one d f
        vector<vector<double>> f = oneD_to_twoD(vals);
        // ready to calc corr func, need vectors to save it
        Eigen::VectorXd C_x = Eigen::VectorXd::Zero((long)nr_dists);
        Eigen::VectorXd C_y = Eigen::VectorXd::Zero((long)nr_dists);

        calc_corr(f, C_x, C_y);
        // average them for now
        Eigen::VectorXd C = average_two(C_x, C_y);
        // vector with distances?
        Eigen::VectorXd dists = Eigen::VectorXd::LinSpaced((long)nr_dists, 0, (int)nr_dists - 1);
        // construct matrix out of the two
        Eigen::MatrixXd C_x_vals(nr_dists, 2);
        C_x_vals = construct_matrix(C, dists);
        // now we need to do the fit
        // init starting values
        Eigen::VectorXd params(2);
        params << 1.0, 5.0;
        NumericalExpDecayDiffFunctor functor(C_x_vals);
        Eigen::LevenbergMarquardt<NumericalExpDecayDiffFunctor> lm(functor);
        Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(params);
        std::cout << "status: " << status << std::endl;

        //std::cout << "info: " << lm.info() << std::endl;

        std::cout << "params that minimizes the function: " << std::endl << params << std::endl;

        // now printing to file
        // TODO i am now done and now i realize that it isnt even that smart to fit for every system...
        // TODO we should actually average the Correlation function and fit afterwards
        // but good that you now know how to fit in c++...
        // params(1) is the correlation length that we want, so we write this to the file
        xi_out << t << "," << params(1) << endl;
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
    string rng = "RNG during Euler";
    string functor = "Functor during Euler";
    checkpoint_timer timer{{system_name, applying_name, rng, functor}};
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
    int switch_counter = 0;
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
                cout << "or are we at prev_k? prev_k = " << prev_accepted_k << endl;
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
            // cout << "done calcing diff" << endl;
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

    euler_combined(size_t N, int K, double tol, int switch_count = 10000, double reduction_factor=1.5) : euler_mayurama_stepper<state_type, algebra, operations, value_type, time_type>(N),
             dx_drift_dt(N), x_drift(N), reduction_factor(reduction_factor),
    k(K), tol(tol), prev_accepted_k(K), switch_count(switch_count)
    {
        cout << "creating euler combined stepper" << endl;
    }

    int get_k() {
        return k;
    }

    double get_error() {
        return error;
    }

private:
    state_type x_drift, dx_drift_dt;
};


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

template<typename T>
std::vector<T> logspace(T start_in, T end_in, int num_in, T base_in = 2.0)
{
    std::vector<T> logspaced;

    double start = static_cast<double>(start_in);
    double end = static_cast<double>(end_in);
    double base = static_cast<double>(base_in);
    double num = static_cast<double>(num_in);

    if (num == 0) { return logspaced; }
    if (num == 1)
    {
        logspaced.push_back((T)pow(base, start));
        return logspaced;
    }

    double delta = (end - start) / (num - 1);

    for(int i=0; i < num-1; ++i)
    {
        logspaced.push_back((T)pow(base, start + delta * i));
    }
    logspaced.push_back(end); // I want to ensure that start and end
    // are exactly the same as the input
    return logspaced;
}

template <typename T>
std::vector<T> geomspace(T start, T stop, int num, bool endpoint = true) {
    std::vector<T> result;

    if (num < 1) {
        return result;
    }

    if (endpoint) {
        T factor = std::pow(stop / start, 1.0 / (num - 1));
        for (int i = 0; i < num; ++i) {
            T value = start * std::pow(factor, i);
            result.push_back(value);
        }
    } else {
        T factor = std::pow(stop / start, 1.0 / num);
        for (int i = 0; i < num; ++i) {
            T value = start * std::pow(factor, i);
            result.push_back(value);
        }
    }

    return result;
}

void print_vector(std::vector<double> vec)
{
    std::cout << "size: " << vec.size() << std::endl;
    for (double d : vec)
        std::cout << d << " ";
    std::cout << std::endl;
}

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
        double rnr = (double)dist(rng);
        // if the drawn number is greater than 3 * sigma we cut it
/*        printf(
        "%f  ", rnr
        );
        printf("\n");*/
        if(abs(rnr) > 2.0 * sigma) {
            rnr = (rnr < 0) ? (-1.0) * 2.0 * sigma : 2.0 * sigma;
/*            printf(
                    "rnr changed to %f", rnr
                    );*/
        }
        return (double)ampl * rnr;
    }
};

template <class State, size_t n>
void fill_init_values(State &state, float x0, float p0, int run = 0, double mu=0, double sigma=1) {

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
/*    for(int i = 0; i < state.size(); i++) {
        cout << "x[" << i << "]=" << state[i] << endl;
    }*/
}


#endif //CUDAPROJECT_MAIN_CUH
