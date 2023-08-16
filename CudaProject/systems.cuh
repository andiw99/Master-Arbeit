//
// Created by andi on 30.04.23.
//

#ifndef CUDAPROJECT_SYSTEMS_CUH
#define CUDAPROJECT_SYSTEMS_CUH

#include "main.cuh"
#include <cmath>
#include <random>
#include <functional>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <boost/typeof/typeof.hpp>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/random.h>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <stdio.h>

using namespace std;

template <size_t lat_dim>
struct System {
public:
    const size_t n;
    // needed to "advance the random engine" and generate different random numbers for the different steps
    size_t step_nr;
    double T;
    const double eta;
    double D;
    checkpoint_timer timer {{}};           // checkpoint timer with no names, for now only total time

    struct rand
    {
        double mu, sigma, D;

        __host__ __device__
        rand(double D, double mu = 0.0, double sigma = 1.0) : D(D), mu(mu), sigma(sigma) {};

        __host__ __device__
        float operator()(const unsigned int ind) const
        {
            thrust::default_random_engine rng;
            thrust::normal_distribution<double> dist(mu, sigma);
            rng.discard(ind);

            return D * dist(rng);
        }
    };
    struct rand_float
    {
        float mu, sigma, D;

        __host__ __device__
        rand_float(float D, float mu = 0.0f, float sigma = 1.0f) : D(D), mu(mu), sigma(sigma) {};

        __host__ __device__
        float operator()(const unsigned int ind) const
        {
            thrust::default_random_engine rng;
            thrust::normal_distribution<float> dist(mu, sigma);
            rng.discard(ind);

            return D * dist(rng);
        }
    };
    /*
    __global__ curandState global_state;
    struct curand {
        float D;
        curand(float D) : D(D) {}

        __host__ __device__
        float operator()() const
        {

            return D;
        }
    };
*/

    struct left : thrust::unary_function<size_t, size_t> {
        __host__ __device__ size_t operator()(size_t i) const {
            // Here we implement logic that return the index of the left neighbor
            // we have to think about that we are actually in 2D and i guess we want to use PBC?
            // so we have to know the system size
            // would we do that with a template oder with a attribute?
            // Another thing is how to implement the logic
            // with modulo we don't need any logic but will this be faster?
            // lat_dim is the sqrt of n
            // size_t j;
            // j is always just i-1, except when it is on the left side of the lattice, then it is i + lat_dim (-1?)
            // if i is on the left side of the lattice, i % lat_dim = 0
            // j = (i % lat_dim == 0) ? i + lat_dim - 1 : i - 1;

            return (i % lat_dim == 0) ? i + lat_dim - 1 : i - 1;
        }
    };

    struct right : thrust::unary_function<size_t, size_t> {
        __host__ __device__ size_t operator()(size_t i) const {
            // j is always i+1, expect when i is on the right side of the lattice
            // if i is on the right side of the lattice, j is i - (d - 1)
            // if i is one the right side of the lattice i % lat_dim = lat_dim - 1

            return (i % lat_dim == lat_dim - 1) ? i - (lat_dim - 1) : i + 1;
        }
    };

    struct up : thrust::unary_function<size_t, size_t> {
        __host__ __device__ size_t operator()(size_t i) const {
            // j is always i - d, except when i is on the upper bound of the lattice
            // if it is on the upper bound, j will be i + d(d-1)
            // if i is on the upper bound, i will be smaller than d
            return (i < lat_dim) ? i + lat_dim * (lat_dim - 1) : i - lat_dim;
        }
    };

    struct down : thrust::unary_function<size_t, size_t> {
        __host__ __device__ size_t operator()(size_t i) const {
            // j is always i + d, except when i is on the lower bound of the lattice
            // if it is on the lower bound, j will be i - d(d-1)
            // if i is on the lower bound, i will be larger than d * (d-1) - 1 = d*d - d - 1
            return (i >= lat_dim * (lat_dim -1 )) ? i - lat_dim * (lat_dim - 1) : i + lat_dim;
        }
    };


    template<class State, class Deriv, class Stoch>
    void operator()(const State &x, Deriv &dxdt, Stoch &theta, double t) {}

    struct Functor {
        Functor(){}

        virtual void operator()() {}

    };

    template<class State, class Deriv, class FunctorType>
    void universalStepOperations(const State &x, Deriv &dxdt, double t, FunctorType functor) {
        BOOST_AUTO(start, thrust::make_zip_iterator(thrust::make_tuple(
                x.begin(),
                x.begin() + n,
                dxdt.begin(),
                dxdt.begin() + n,
                thrust::make_permutation_iterator(
                        x.begin(),
                        thrust::make_transform_iterator(
                                thrust::counting_iterator<size_t>(0),
                                left()
                        )
                ),
                thrust::make_permutation_iterator(
                        x.begin(),
                        thrust::make_transform_iterator(
                                thrust::counting_iterator<size_t>(0),
                                right()
                        )
                ),
                thrust::make_permutation_iterator(
                        x.begin(),
                        thrust::make_transform_iterator(
                                thrust::counting_iterator<size_t>(0),
                                up()
                        )
                ),
                thrust::make_permutation_iterator(
                        x.begin(),
                        thrust::make_transform_iterator(
                                thrust::counting_iterator<size_t>(0),
                                down()
                        )
                )
        )));
        thrust::for_each(start, start + n, functor);
        step_nr++;
    }

    template<class Stoch>
    void calc_diff(Stoch &theta, double t) {
        // TODO actually t is never needed here with the current architecture, but i am too lazy to fix that
        // as it will probably improve nothing
        thrust::counting_iterator<size_t> index_sequence_begin(step_nr * n);
        if(step_nr == 0) {
            thrust::fill(theta.begin(), theta.begin() + n, 0);
        }
        thrust::transform(index_sequence_begin,
                          index_sequence_begin + n,
                          theta.begin() + n,
                          rand(D));
    }

    System(size_t step_nr, const double eta, const double T) : step_nr(step_nr), n(lat_dim * lat_dim), eta(eta), T(T), D(sqrt(2 * T * eta)) {
        cout << "System constructor is called with eta = " << eta << "  T = " << T << endl;
    }
    System() {
        cout << "Probably wrong system constructor called" << endl;
    }

    size_t get_lattice_dim() const{
        return lat_dim;
    }

    double get_cur_T() const{
        // If we Quench, we just have to keep this T up to date and we do that I think, then we can reuse this function
        return T;
    }

    double get_T(double t) {
        return T;
    }

    template <class T>
    struct square
    {
        __host__ __device__
        T operator()(const T& x) const {
            return x * x;
        }
    };

    template<class State>
    double calc_kinetic_energy(State &x) {
        double E_kin = 0.5 * thrust::transform_reduce(x.begin() + n, x.end(), square<double>(), 0.0, thrust::plus<double>());
        return E_kin;
    }

    size_t get_step_nr() {
        return step_nr;
    }
};




template <size_t lat_dim>
struct constant_bath : public System<lat_dim> {
    using left = typename System<lat_dim>::left;
    using right = typename System<lat_dim>::right;
    using up = typename System<lat_dim>::up;
    using down = typename System<lat_dim>::down;
public:
    const double alpha;
    const double beta;
    const double J;
    // systemsize should probably be a template argument?
    // In contrast to the brownian system we need one more parameter for the timescale of the cooling down
    // the current temperature


    string rng = "RNG";
    string theta_filling = "Filling of theta";
    string functor_point = "Functor Calc";
    checkpoint_timer timer {{rng, functor_point, theta_filling}};
    // parameters of the potential and of the Interaction
    struct bath_functor {
        // I think also the potential and interaction parameters have to be set in the functor
        // I mean i could template everything and this would probably also give a bit of potential but is it really
        // worth it?
        // why would you not use templates in c++ instead of parameters, only when the parameter is not clear at
        // runtime, am i right? since lattice size, potential parameters, etc. don't change during runtime we
        // could just template everything
        const double alpha, beta, J, eta;

        bath_functor(const double eta, const double alpha,
                     const double beta, const double J) : alpha(alpha), beta(beta), J(J), eta(eta) { }

        template<class Tup>
        __host__ __device__ void operator()(Tup tup) {
            double q = thrust::get<0>( tup );
            double p = thrust::get<1>( tup );
            thrust::get<2>( tup ) = p;
            double q_left = thrust::get<4>(tup);
            double q_right = thrust::get<5>(tup);
            double q_up = thrust::get<6>(tup);
            double q_down = thrust::get<7>(tup);
            double interaction = J * ((q - q_left) + (q - q_right) + (q - q_up) + (q - q_down));
            thrust::get<3>( tup ) = (-eta) * p                                                                                  // Friction
                                    - alpha * (2 * q * q * q - beta * q)                                                        // double well potential
                                    - interaction;       // Interaction
        }
    };


    template<class State, class Deriv>
    void calc_drift(const State &x, Deriv &dxdt, double t) {
        // cout << "calc_drift is called" << endl;
        timer.set_startpoint(functor_point);

        bath_functor functor = bath_functor(System<lat_dim>::eta, alpha, beta, J);

        this->universalStepOperations(x, dxdt, t, functor);
        timer.set_endpoint(functor_point);
    }

    template<class Stoch>
    void calc_diff(Stoch &theta, double t) {
        // cout << "calc_diff is called" << endl;
        timer.set_startpoint(rng);
        System<lat_dim>::calc_diff(theta, t);
        timer.set_endpoint(rng);
    }
public:
    constant_bath(const double T, const double eta, const double alpha, const double beta, const double J, const int init_step=0)
            : System<lat_dim>(init_step, eta, T), alpha(alpha), beta(beta), J(J) {
        cout << "constant Bath System is constructed" << endl;
    }

    double get_cur_T() const{
        return System<lat_dim>::T;
    }
    ~constant_bath() {
        cout << "Bath System is destroyed" << endl;
    }
};


template <size_t lat_dim>
struct coulomb_interaction : public System<lat_dim> {
    using left = typename System<lat_dim>::left;
    using right = typename System<lat_dim>::right;
    using up = typename System<lat_dim>::up;
    using down = typename System<lat_dim>::down;
public:
    const double alpha;
    const double beta;
    const double J;
    // systemsize should probably be a template argument?
    // needed to "advance the random engine" and generate different random numbers for the different steps
    // In contrast to the brownian system we need one more parameter for the timescale of the cooling down
    // the current temperature
    // parameters of the potential and of the Interaction
    struct coulomb_functor {
        // I think also the potential and interaction parameters have to be set in the functor
        // I mean i could template everything and this would probably also give a bit of potential but is it really
        // worth it?
        // why would you not use templates in c++ instead of parameters, only when the parameter is not clear at
        // runtime, am i right? since lattice size, potential parameters, etc. don't change during runtime we
        // could just template everything
        const double alpha, beta, J, eta;

        coulomb_functor(const double eta, const double alpha,
                        const double beta, const double J) : alpha(alpha), beta(beta), J(J), eta(eta) { }

        template<class Tup>
        __host__ __device__ double coulomb_interaction(Tup tup) {
            double q = thrust::get<0>( tup );
            double q_left = thrust::get<4>(tup);
            double q_right = thrust::get<5>(tup);
            double q_up = thrust::get<6>(tup);
            double q_down = thrust::get<7>(tup);

            return J * (
                    ((q - q_left) / pow(1.0 + (q - q_left) * (q - q_left), 1.5))
                    +   ((q - q_right) / pow(1.0 + (q - q_right) * (q - q_right), 1.5))
                    +   ((q - q_up) / pow(1.0 + (q - q_up) * (q - q_up), 1.5))
                    +   ((q - q_down) / pow(1.0 + (q - q_down) * (q - q_down), 1.5))
            );
        }
        template<class Tup>
        __host__ __device__ void operator()(Tup tup) {
            double q = thrust::get<0>( tup );
            double p = thrust::get<1>( tup );
            thrust::get<2>( tup ) = p;
            double q_left = thrust::get<4>(tup);
            double q_right = thrust::get<5>(tup);
            double q_up = thrust::get<6>(tup);
            double q_down = thrust::get<7>(tup);
            double interaction = J * (
                        ((q - q_left)   / pow(1.0 + (q - q_left)    * (q - q_left),  1.5))
                    +   ((q - q_right)  / pow(1.0 + (q - q_right)   * (q - q_right), 1.5))
                    +   ((q - q_up)     / pow(1.0 + (q - q_up)      * (q - q_up),    1.5))
                    +   ((q - q_down)   / pow(1.0 + (q - q_down)    * (q - q_down),  1.5))
            );
            thrust::get<3>( tup ) = (-eta) * p                                                                                  // Friction
                                    - alpha * (2 * q * q * q - beta * q)                                                        // double well potential
                                    - interaction;       // Interaction
        }
    };

public:
    template<class State, class Deriv>
    void calc_drift(const State &x, Deriv &dxdt, double t) {
        coulomb_functor functor = coulomb_functor(System<lat_dim>::eta, alpha, beta, J);
        this->universalStepOperations(x, dxdt, t, functor);
    }



    coulomb_interaction(const double T, const double eta, const double alpha, const double beta, const double J, const int init_step=0)
            : System<lat_dim>(init_step, eta, T), alpha(alpha),
            beta(beta), J(J) {
    }
};

template <size_t lat_dim>
class coulomb_constant : public coulomb_interaction<lat_dim> {
public:
    coulomb_constant(const double T, const double eta, const double alpha, const double beta, const double J, const int init_step=0)
    : coulomb_interaction<lat_dim>(T, eta, alpha, beta, J, init_step) {

    }
    template<class Stoch>
    void calc_diff(Stoch &theta, double t) {
        // cout << "calc_diff is called" << endl;
        System<lat_dim>::calc_diff(theta, t);
    }
};


template <size_t lat_dim>
struct gpu_bath : public coulomb_interaction<lat_dim> {
    // actually I think we can write the gpu_bath as extension of coulomb or coulomb constant?
public:
    const double T_start;       // Start-Temperture of the Quenching For example: 10
    const double T_end;         // End-Temperature of the Quencheing. For example: 1
    const double s_eq_t;        // start equilibration time: Time the system gets to equilibrate at T_start
    const double e_eq_t;        // end equilibration time: Time the system gets to equilibrate at T_end (after Quench)
    double t_quench;            // total time the quench takes, is determined by the quench timescale tau and the temp difference
    double end_quench_t;        // timepoint at which the Quenching ends = s_eq_t + t_quench
    const double tau;
    // systemsize should probably be a template argument?
    // In contrast to the brownian system we need one more parameter for the timescale of the cooling down
    // the current temperature
    // parameters of the potential and of the Interaction

    double simple_linear_T(double t) {
        // parametrisierung für die Temperatur
        // linearer Abfall
        System<lat_dim>::T = max(T_start - t/tau, T_end);
        return System<lat_dim>::T;
    }

    double linear_T(double t) {
        if(s_eq_t < t && t < end_quench_t) {
            // if we are in the quench phase, we reduce T
            System<lat_dim>::T = T_start - (t - s_eq_t)/tau;
        }
        return System<lat_dim>::T;
    }

    template<class Stoch>
    void calc_diff(Stoch &theta, double t) {
        // I think this should work? We change the D of the System and then just calc the random numbers
        System<lat_dim>::D = sqrt(2 * linear_T(t) * System<lat_dim>::eta);
        System<lat_dim>::calc_diff(theta, t);
    }
public:
    gpu_bath(const double T, const double T_end, const double eta, const double alpha, const double beta, const double J, const double tau, size_t init_step = 0, double eq_t = 30)
            : coulomb_interaction<lat_dim>(T, eta, alpha, beta, J, init_step), T_start(T),
              tau(tau), T_end(T_end), s_eq_t(eq_t), e_eq_t(eq_t) {
        t_quench = (get_quench_time());
        end_quench_t = t_quench + s_eq_t;       // end of quench is the quench time + equilibrate at beginning
    }

    double get_quench_time() {
        // returns the time it takes to do the quench
        // in this system, we use a linear quench
        return (T_start - T_end) * tau;
    }

    double get_end_t(){
        // the total time are the two equilibriate times + the quench time
        return s_eq_t + e_eq_t + t_quench;
    }

    double get_end_quench_time() {
        return end_quench_t;
    }
};

template <size_t lat_dim>
struct quench : virtual public System<lat_dim> {
    const double T_start;       // Start-Temperture of the Quenching For example: 10
    const double T_end;         // End-Temperature of the Quencheing. For example: 1
    const double s_eq_t;        // start equilibration time: Time the system gets to equilibrate at T_start
    const double e_eq_t;        // end equilibration time: Time the system gets to equilibrate at T_end (after Quench)
    double t_quench;            // total time the quench takes, is determined by the quench timescale tau and the temp difference
    double end_quench_t;        // timepoint at which the Quenching ends = s_eq_t + t_quench
    const double tau;

    double linear_T(double t) {
        if(s_eq_t < t && t < end_quench_t) {
            // if we are in the quench phase, we reduce T
            System<lat_dim>::T = T_start - (t - s_eq_t)/tau;
        }
        return System<lat_dim>::T;
    }

    template<class Stoch>
    void calc_diff(Stoch &theta, double t) {
        // I think this should work? We change the D of the System and then just calc the random numbers
        System<lat_dim>::D = sqrt(2 * linear_T(t) * System<lat_dim>::eta);
        System<lat_dim>::calc_diff(theta, t);
    }

public:
   quench(const double T, const double T_end, const double eta, const double tau, size_t init_step = 0, double eq_t = 30)
            : System<lat_dim>(init_step, eta, T), T_start(T), tau(tau), T_end(T_end), s_eq_t(eq_t), e_eq_t(eq_t) {
        t_quench = (get_quench_time());
        end_quench_t = t_quench + s_eq_t;       // end of quench is the quench time + equilibrate at beginning
    }

    double get_quench_time() {
        // returns the time it takes to do the quench
        // in this system, we use a linear quench
        return (T_start - T_end) * tau;
    }

    double const get_end_t() const{
        // the total time are the two equilibriate times + the quench time
        return s_eq_t + e_eq_t + t_quench;
    }

    double get_end_quench_time() {
        return end_quench_t;
    }
};

template <size_t lat_dim>
struct anisotropic_coulomb_interaction : virtual public System<lat_dim> {
    using left = typename System<lat_dim>::left;
    using right = typename System<lat_dim>::right;
    using up = typename System<lat_dim>::up;
    using down = typename System<lat_dim>::down;
public:
    const double alpha;
    const double beta;
    const double Jx, Jy;
    struct ani_coulomb_functor {
        const double alpha, beta, Jx, Jy, eta;

        ani_coulomb_functor(const double eta, const double alpha,
                        const double beta, const double Jx, const double Jy) : alpha(alpha), beta(beta), Jx(Jx), Jy(Jy), eta(eta) { }

        template<class Tup>
        __host__ __device__ void operator()(Tup tup) {
            double q = thrust::get<0>( tup );
            double p = thrust::get<1>( tup );
            thrust::get<2>( tup ) = p;
            double q_left = thrust::get<4>(tup);
            double q_right = thrust::get<5>(tup);
            double q_up = thrust::get<6>(tup);
            double q_down = thrust::get<7>(tup);
            double interaction_x = Jx * (
                    ((q - q_left)   / pow(1.0 + (q - q_left)    * (q - q_left),  1.5))
                    +   ((q - q_right)  / pow(1.0 + (q - q_right)   * (q - q_right), 1.5)));
            double interaction_y = Jy * (
                    +   ((q - q_up)     / pow(1.0 + (q - q_up)      * (q - q_up),    1.5))
                    +   ((q - q_down)   / pow(1.0 + (q - q_down)    * (q - q_down),  1.5))
            );
            thrust::get<3>( tup ) = (-eta) * p                                                                                  // Friction
                                    - alpha * (2 * q * q * q - beta * q)                                                        // double well potential
                                    - interaction_x - interaction_y;       // Interaction
        }
    };

public:
    template<class State, class Deriv>
    void calc_drift(const State &x, Deriv &dxdt, double t) {
        ani_coulomb_functor functor = ani_coulomb_functor(System<lat_dim>::eta, alpha, beta, Jx, Jy);
        this->universalStepOperations(x, dxdt, t, functor);
    }

    anisotropic_coulomb_interaction(const double T, const double eta, const double alpha, const double beta, const double Jx, const double Jy, const int init_step=0)
            : System<lat_dim>(init_step, eta, T), alpha(alpha),
              beta(beta), Jx(Jx), Jy(Jy) {
    }
};


template <size_t lat_dim>
class anisotropic_coulomb_constant : public anisotropic_coulomb_interaction<lat_dim> {
public:
    anisotropic_coulomb_constant(const double T, const double eta, const double alpha, const double beta, const double Jx, const double Jy, const int init_step=0)
            : System<lat_dim>(init_step, eta, T), anisotropic_coulomb_interaction<lat_dim>(T, eta, alpha, beta, Jx, Jy, init_step) {

    }
    template<class Stoch>
    void calc_diff(Stoch &theta, double t) {
        // cout << "calc_diff is called" << endl;
        System<lat_dim>::calc_diff(theta, t);
    }
};

template <size_t lat_dim>
class anisotropic_coulomb_quench: public anisotropic_coulomb_interaction<lat_dim>, public quench<lat_dim> {
public:
    anisotropic_coulomb_quench(const double T, const double T_end, const double eta, const double alpha,
                               const double beta, const double Jx,const double Jy, const double tau,
                               size_t init_step = 0, double eq_t = 30)
            : anisotropic_coulomb_interaction<lat_dim>(T, eta, alpha, beta, Jx, Jy, init_step),
                    quench<lat_dim>(T, T_end, eta, tau, init_step, eq_t),
                            System<lat_dim>(init_step, eta, T){

    }
};

/*template <size_t lat_dim>
struct sys_factory {
public:
    template<template<class> class sys>
    sys<size_t> init_sys(map<string, double> &paras) {
    }

    template<>
    anisotropic_coulomb_quench<lat_dim> init_sys<anisotropic_coulomb_quench<lat_dim >>(map<string, double> &paras) {
        return anisotropic_coulomb_quench<lat_dim>(paras["starting_T"],
                                                   paras["end_T"],
                                                   paras["eta"],
                                                   paras["alpha"],
                                                   paras["beta"],
                                                   paras["J"],
                                                   paras["Jy"],
                                                   paras["tau"],
                                                   paras["seed"],
                                                   paras["t_eq"]);
    }

    template<>
    gpu_bath<lat_dim> init_sys<gpu_bath<lat_dim>>(map<string, double> &paras) {
        return gpu_bath<lat_dim>(paras["starting_T"],
                                 paras["end_T"],
                                 paras["eta"],
                                 paras["alpha"],
                                 paras["beta"],
                                 paras["J"],
                                 paras["tau"],
                                 paras["seed"],
                                 paras["t_eq"]);
    }
};*/

template<size_t lat_dim, template<size_t> class sys>
sys<lat_dim> create(map<string, double> &paras, size_t seed){
    size_t seed_val = seed;
    if(std::is_same<sys<lat_dim>, anisotropic_coulomb_constant<lat_dim>>::value) {
            cout << std::is_same<sys<lat_dim>, gpu_bath<lat_dim>>::value << endl;
            cout << std::is_same<sys<lat_dim>, anisotropic_coulomb_quench<lat_dim>>::value << endl;
            cout << "Create function for anisotropic coulomb constant called" << endl;
            return sys<lat_dim>(paras["T"],
                              paras["eta"],
                              paras["alpha"],
                              paras["beta"],
                              paras["J"],
                              paras["Jy"],
                              seed_val);}
       if(std::is_same<sys<lat_dim>, gpu_bath<lat_dim>>::value) {
        cout << "Create function for gpu_bath called" << endl;
        return sys<lat_dim>(paras["starting_T"],
                            paras["end_T"],
                            paras["eta"],
                            paras["alpha"],
                            paras["beta"],
                            paras["J"],
                            paras["tau"],
                            seed_val,
                            paras["t_eq"]);
    } else if(std::is_same<sys<lat_dim>, anisotropic_coulomb_quench<lat_dim>>::value) {
        cout << "Create function for anisotropic_coulomb_quench called" << endl;
        return sys<lat_dim>(paras["starting_T"],
                            paras["end_T"],
                            paras["eta"],
                            paras["alpha"],
                            paras["beta"],
                            paras["J"],
                            paras["Jy"],
                            paras["tau"],
                            seed_val,
                            paras["t_eq"]);
    }
};

template<size_t lat_dim, template<size_t> class sys>
sys<lat_dim> create(map<string, double> &paras){
    return create<lat_dim, sys>(paras, paras["seed"]);
}



/*

template <template <class> class sys>
sys<size_t> create()
{
    sys<size_t>* object = new sys<size_t>();

    // Do some generic initialization after construction.

    return *object;
};
*/

/*template<class sys>
sys init_sys(map<string, double> &paras) {
    return NULL;
}

template<>
anisotropic_coulomb_quench init_sys<anisotropic_coulomb_interaction>(map<string, double> &paras) {
    return anisotropic_coulomb_quench<lat_dim>(paras["starting_T"],
                                      paras["end_T"],
                                      paras["eta"],
                                      paras["alpha"],
                                      paras["beta"],
                                      paras["J"],
                                      paras["Jy"],
                                      paras["tau"],
                                      paras["seed"],
                                      paras["t_eq"]);
}

template <size_t lat_dim>
gpu_bath<lat_dim> init_sys(map<string, double> &paras) {
    return gpu_bath<lat_dim>(paras["starting_T"],
                                               paras["end_T"],
                                               paras["eta"],
                                               paras["alpha"],
                                               paras["beta"],
                                               paras["J"],
                                               paras["tau"],
                                               paras["seed"],
                                               paras["t_eq"]);
}*/

template <size_t lat_dim>
struct quadratic_trapped_lattice : public System<lat_dim> {
    using rand = typename System<lat_dim>::rand;
    using left = typename System<lat_dim>::left;
    using right = typename System<lat_dim>::right;
    using up = typename System<lat_dim>::up;
    using down = typename System<lat_dim>::down;
public:
    const double alpha;
    const double J;
    // systemsize should probably be a template argument?
    // In contrast to the brownian system we need one more parameter for the timescale of the cooling down
    // the current temperature
    string rng = "RNG";
    string theta_filling = "Filling of theta";
    string functor_point = "Functor Calc";
    checkpoint_timer timer {{rng, functor_point, theta_filling}};
    // parameters of the potential and of the Interaction
    struct harmonic_trap_functor {
        // I think also the potential and interaction parameters have to be set in the functor
        // I mean i could template everything and this would probably also give a bit of potential but is it really
        // worth it?
        // why would you not use templates in c++ instead of parameters, only when the parameter is not clear at
        // runtime, am i right? since lattice size, potential parameters, etc. don't change during runtime we
        // could just template everything
        const double alpha, J, eta;

        harmonic_trap_functor(const double eta, const double alpha,
                      const double J) : alpha(alpha), J(J), eta(eta) { }

        template<class Tup>
        __host__ __device__ void operator()(Tup tup) {
            double q = thrust::get<0>( tup );
            double p = thrust::get<1>( tup );
            thrust::get<2>( tup ) = p;
            double q_left = thrust::get<4>(tup);
            double q_right = thrust::get<5>(tup);
            double q_up = thrust::get<6>(tup);
            double q_down = thrust::get<7>(tup);
            double interaction = J * ((q - q_left) + (q - q_right) + (q - q_up) + (q - q_down));

            thrust::get<3>( tup ) = (-eta) * p                                                                                  // Friction
                                    - alpha * (2 * q)                                                        // double well potential
                                    + interaction;       // Interaction
        }
    };


    template<class State, class Deriv>
    void calc_drift(const State &x, Deriv &dxdt, double t) {
        timer.set_startpoint(functor_point);

        harmonic_trap_functor functor = harmonic_trap_functor(System<lat_dim>::eta, alpha, J);

        this->universalStepOperations(x, dxdt, t, functor);
        timer.set_endpoint(functor_point);
    }

    template<class Stoch>
    void calc_diff(Stoch &theta, double t) {
        timer.set_startpoint(rng);
        System<lat_dim>::calc_diff(theta, t);
        timer.set_endpoint(rng);
    }

public:
    quadratic_trapped_lattice(const double T, const double eta, const double alpha, const double J, const int init_step=0)
            : System<lat_dim>(init_step, eta, T), alpha(alpha), J(J) {
        cout << "Lattice quadratic Trap System is constructed" << endl;
    }

    double get_cur_T() const{
        return System<lat_dim>::T;
    }

    size_t get_lattice_dim() const{
        return lat_dim;
    }

};


template <size_t lat_dim>
struct quadratic_chain {
    using rand = typename System<lat_dim>::rand;
public:
    const double eta;
    const double J;
    const double D;
    // systemsize should probably be a template argument?
    const size_t n;
    // needed to "advance the random engine" and generate different random numbers for the different steps
    size_t step_nr;
    // In contrast to the brownian system we need one more parameter for the timescale of the cooling down
    // the current temperature
    double T;
    // parameters of the potential and of the Interaction
    struct functor {
        // I think also the potential and interaction parameters have to be set in the functor
        // I mean i could template everything and this would probably also give a bit of potential but is it really
        // worth it?
        // why would you not use templates in c++ instead of parameters, only when the parameter is not clear at
        // runtime, am i right? since lattice size, potential parameters, etc. don't change during runtime we
        // could just template everything
        const double J, eta;

        functor(const double eta, const double J) : J(J), eta(eta) { }

        template<class Tup>
        __host__ __device__ void operator()(Tup tup) {
            double q = thrust::get<0>( tup );
            double p = thrust::get<1>( tup );
            thrust::get<2>( tup ) = p;
            double q_left = thrust::get<4>(tup);
            double q_right = thrust::get<5>(tup);


            thrust::get<3>( tup ) = (-eta) * p                                                                                  // Friction
                                    - J * ((q - q_left) + (q - q_right));       // Interaction
        }
    };

    struct left : thrust::unary_function<size_t, size_t> {
        __host__ __device__ size_t operator()(size_t i) const {
            // if we are at the left end of our chain, we need to use the right end as left neighbor
            return (i == 0) ? i + lat_dim - 1 : i - 1;
        }
    };

    struct right : thrust::unary_function<size_t, size_t> {
        __host__ __device__ size_t operator()(size_t i) const {
            return (i == lat_dim - 1) ? i - (lat_dim - 1) : i + 1;
        }
    };


    template<class State, class Deriv>
    void calc_drift(const State &x, Deriv &dxdt, double t) {


        BOOST_AUTO(start, thrust::make_zip_iterator(thrust::make_tuple(
                x.begin(),
                x.begin() + n,
                dxdt.begin(),
                dxdt.begin() + n,
                thrust::make_permutation_iterator(
                        x.begin(),
                        thrust::make_transform_iterator(
                                thrust::counting_iterator<size_t>(0),
                                left()
                        )
                ),
                thrust::make_permutation_iterator(
                        x.begin(),
                        thrust::make_transform_iterator(
                                thrust::counting_iterator<size_t>(0),
                                right()
                        )
                )
        )));
        thrust::for_each(start, start + n, functor(eta, J));
        step_nr++;
        // cout << "x[1] = " << x[1] << endl;
        // the problem is here, actually dxdt[0] = x[1] should be
        // somehow dxdt is not set but i dont really get why
        // cout << "dxdt[0] = " << dxdt[0] << endl;
    }

    template<class Stoch>
    void calc_diff(Stoch &theta, double t) {
        thrust::counting_iterator<size_t> index_sequence_begin(step_nr * n);
        if(step_nr == 0) {
            // TODO will this already be initialized to zero without this statement?
            thrust::fill(theta.begin(), theta.begin() + n, 0);
        }

        thrust::transform(index_sequence_begin,
                          index_sequence_begin + n,
                          theta.begin() + n,
                          rand(D));
    }

public:
    quadratic_chain(const double T, const double eta, const double J, const int init_step=0)
            : T(T), step_nr(init_step), n(lat_dim), eta(eta), J(J), D(sqrt(2 * T * eta)) {
    }

    double get_cur_T() const{
        return T;
    }

    size_t get_lattice_dim() const{
        return lat_dim;
    }

    template<class State>
    double calc_energy(const State &x) {
        double E = 0;
        for(int i = 0; i < lat_dim; i++) {
            // add potential energy
            E += J/2 * (x[i] - x[(i + 1) % lat_dim]) * (x[i] - x[(i + 1) % lat_dim]);
            // add kinetic energy
            E += 0.5 * x[i + n] * x[i + n];
        }
        return E;
    }

    template<class State>
    double calc_kinetic_energy(const State &x) {
        double Ekin = 0;
        for(int i = 0; i < lat_dim; i++) {
            // add kinetic energy
            Ekin += (0.5 * (x[i + n] * x[i + n]));
        }
        return Ekin;
    }

    template<class State>
    double calc_potential_energy(const State &x) {
        double Epot = 0;
        for(int i = 0; i < lat_dim; i++) {
            // add kinetic energy
            Epot += J/2 * (x[i] - x[(i + 1) % lat_dim]) * (x[i] - x[(i + 1) % lat_dim]);

        }
        return Epot;
    }
    template<class State>
    double calc_total_squared_dist(const State &x) {
        double d2 = 0;
        for(int i = 0; i < lat_dim - 1; i++) {
            // add kinetic energy
            d2 += (x[i] - x[i + 1]) * (x[i] - x[i + 1]);
        }
        return d2;
    }


};


template <size_t lat_dim>
struct gpu_oscillator_chain {
    using rand = typename System<lat_dim>::rand;
public:
    const double T;
    const double alpha;
    const double eta;
    // no tau, constant temp
    // no beta, just harmonic potential
    // no interaction
    // systemsize should probably be a template argument?
    const size_t n;
    size_t step_nr;

    struct oscillator_chain_functor {

        const double alpha, eta;

        oscillator_chain_functor(const double eta, const double alpha) : alpha(alpha), eta(eta) { }

        template<class Tup>
        __host__ __device__ void operator()(Tup tup) {
            double q = thrust::get<0>( tup );
            double p = thrust::get<1>( tup );
            thrust::get<2>( tup ) = p;
            thrust::get<3>( tup ) = (-eta) * p                                                                                  // Friction
                                    - alpha * q;                                                         // harmonic potential
        }
    };



    template<class State, class Deriv, class Stoch>
    void operator()(const State &x, Deriv &dxdt, Stoch &theta, double t) {
        thrust::counting_iterator<size_t> index_sequence_begin(step_nr * n);
        if(step_nr == 0) {
            thrust::fill(theta.begin(), theta.begin() + n, 0);
        }

        double D = sqrt(2 * T * eta);

        thrust::transform(index_sequence_begin,
                          index_sequence_begin + n,
                          theta.begin() + n,
                          rand(D));

        BOOST_AUTO(start, thrust::make_zip_iterator(thrust::make_tuple(
                x.begin(),
                x.begin() + n,
                dxdt.begin(),
                dxdt.begin() + n
        )));
        thrust::for_each(start, start + n, oscillator_chain_functor(eta, alpha));
        step_nr++;
    }
public:
    gpu_oscillator_chain(const double T, const double eta, const double alpha)
            : T(T), step_nr(0), n(lat_dim * lat_dim), eta(eta), alpha(alpha) {
    }

    double get_cur_T() const{
        return T;
    }

    size_t get_lattice_dim() const{
        return lat_dim;
    }
};





struct brownian_particel {
    const double eta, T;
    static inline auto normal_dice = bind(normal_distribution<double>(0, 1), default_random_engine(13));

    brownian_particel(const double eta, const double T)
            : eta(eta), T(T) { }

    /*
     * call of the system needs to take in the theta state type
     */
    template <class state_type>
    void operator()(const state_type &x, state_type &dxdt, state_type &theta, double t)
    {
        dxdt[0] = x[1];
        dxdt[1] = (-eta) * x[1];
        theta[0] = 0;
        theta[1] = sqrt(2 * eta * T) * normal_dice();
    }
};

/*
 * The most difficult part seems to be the implementation of the system, since we need to use thrust operations
 */
typedef thrust::device_vector<double> gpu_state_type;
struct gpu_brownian_system {
public:

    const double eta, T;
    // i don't fucking know, maybe also as template
    // systemsize
    const size_t n;
    // needed to "advance the random engine" and generate different random numbers for the different steps
    size_t step_nr;
    struct brownian_functor {
        const double eta, T;
        //const size_t step_nr;

        brownian_functor(const double eta, const double T)
                : eta(eta), T(T) { }
        template<class Tup>
        __host__ __device__ void operator()(Tup tup) {
            thrust::random::default_random_engine rng;
            thrust::random::normal_distribution<double> dist(0.0, 1.0);
            // the first n=2 (q, p) entries of the tuple correspond to x, the second pair to dxdt and the third to theta
            double q = thrust::get<0>( tup );
            double p = thrust::get<1>( tup );
            // this line has to be the dgl for x
            // dqdt =
            thrust::get<2>( tup ) = p;
            // this the one for p
            // dpdt=
            thrust::get<3>( tup ) = (-eta) * p;
            // and now the ones for theta
            // thrust::get<4>( tup ) = 0;
            // okay, here random number generation, this might be a bit problematic
            // this won't work on the device. This is good i think? Seems to me because multiple random numbers are
            // generated
            // I think I am going to continue tomorrow? I don't think that I am able to solve this now
            // cout << dist(rng) << endl;
            // if i write it like this, i think i set the value of every lattice to be dist(rng) so this is not what
            // I want. I think i need to generate somewhere a vector with 2n random numbers and feed it into the functor?
            // rng.discard(step_nr);
            // thrust::get<5>( tup ) = sqrt(2 * T * eta) * dist(rng);
            // and this is the third that i dont need
            // thrust::get<6>( t ) = -b * z + x * y;
        }
    };

    struct rand
    {
        double mu, sigma, D;

        __host__ __device__
        rand(double D, double mu = 0.0, double sigma = 1.0) : D(D), mu(mu), sigma(sigma) {};

        __host__ __device__
        float operator()(const unsigned int ind) const
        {
            thrust::default_random_engine rng;
            thrust::normal_distribution<double> dist(mu, sigma);
            rng.discard(ind);

            return D * dist(rng);
        }
    };

    template<class State, class Deriv, class Stoch>
    void operator()(const State &x, Deriv &dxdt, Stoch &theta, double t) {
        // I think we need to generate the random numbers here and put them into the start iterator so that we can set
        // the theta values. Another possibility would be to just not put theta into the start iterator at all?
        // I mean it IS completely independent from the rest
        // how do I know whether the random numbers that i generate here are generated in parallel and on the gpu?
        thrust::counting_iterator<size_t> index_sequence_begin(step_nr * n);
        // What i want to do now is fill theta half with zeros
        // I actually only have to do this once. maybe when initializing
        // i mean i could use if..
        if(step_nr == 0) {
            // otherwise it will be filled with zeros since it isn't changed
            // TODO will this already be initialized to zero without this statement?
            thrust::fill(theta.begin(), theta.begin() + n, 0);
        }
        // and the other half with random numbers with eta, T
        // i need to use transform then
        // i think the easiest would be if i give eta and T to rand, but I really cannot judge how this will be
        // performance wise
        // i will at least do the computation only once
        // could be const at the moment but in later Stuff T will change.
        double D = sqrt(2 * T * eta);
        thrust::transform(index_sequence_begin,
                          index_sequence_begin + n,
                          theta.begin() + n,
                          rand(D));
        // Now i hope theta will just be fine and I don't even have to put it into start?
        BOOST_AUTO(start, thrust::make_zip_iterator(thrust::make_tuple(
                                                            // x begin has all q values
                                                            x.begin(),
                                                            // x begin + n has all p values
                                                            x.begin() + n,
                                                            // dxdt begin has all dqdt values
                                                            dxdt.begin(),
                                                            // dxdt begin + n has all dpdt values
                                                            dxdt.begin() + n
                                                            // theta begin hass all theta[0] values
                                                            // theta.begin(),
                                                            //theta begin + n has all theta[1] values
                                                            // theta.begin() + n
                                                    )
        )
        );
        thrust::for_each(start, start + n, brownian_functor(eta, T));
        // this line also doesn't work
        // cout << theta[n] << ", "  << theta[n+1] << endl;
        // advance step, can't there be a better way? Is this inefficient? We lose the const because of that
        // We can also somehow make it a parameter of the stepper
        step_nr++;
    }
public:
    gpu_brownian_system(const double eta, const double T, const size_t n)
            : eta(eta), T(T), n(n), step_nr(0) { }
};

#endif //CUDAPROJECT_SYSTEMS_CUH