//
// Created by andi on 30.04.23.
//

#ifndef CUDAPROJECT_SYSTEMS_CUDA_CUH
#define CUDAPROJECT_SYSTEMS_CUDA_CUH

#include "main-cuda.cuh"

// #include "parameters.cuh"
#include <cmath>
#include <random>
#include <functional>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <boost/typeof/typeof.hpp>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/random.h>
#include <thrust/iterator/discard_iterator.h>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <stdio.h>
#include <curand_kernel.h>
#include <cuda_runtime.h>
#include <cufft.h>


// #include "cufft_utils.h"

using namespace std;

struct System {
public:
    const size_t n;

    // needed to "advance the random engine" and generate different random numbers for the different steps
    double T, D, end_t;
    const double eta;
    // possibility of different sizes in the different directions
    const size_t dim_size_x;
    const size_t dim_size_y;
    thrust::device_vector<curandState> curand_states;
    bool curand_random = false;
    bool equilibrated = false;

    template<class State>
    void init_state(map<Parameter, double>& paras, State &x) {
        cout << "standard init state called" << endl;
    }

    virtual void print_info() {
        cout << "Printing system info:" << endl;
        cout << "T = " << T << endl;
        cout << "eta = " << eta << endl;
        cout << "D = " << D << endl;
        cout << "dim_size_x = " << dim_size_x << endl;
        cout << "dim_size_y = " << dim_size_y << endl;
        cout << "n = " << n << endl;
    }

    struct curand
    {
        double mu, sigma, D;

        __host__ __device__
        curand(double D, double mu = 0.0, double sigma = 1.0) : D(D), mu(mu), sigma(sigma) {};

        template<class Tuple>
        __device__
        float operator()(Tuple tup) const
        {
            curandState local_state = thrust::get<1>(tup);
            thrust::get<0>(tup) = D * curand_normal(&local_state);
            thrust::get<1>(tup) = local_state;
        }
    };

    struct curand_setup
    {
        long seed = 0;

        curand_setup(long seed): seed(seed) {}

        using init_tuple = thrust::tuple<int, curandState &>;
        __device__
        void operator()(init_tuple t) const{
            curandState s;
            int id = thrust::get<0>(t);
            curand_init(seed, id, 0, &s);
            thrust::get<1>(t) = s;
        }
    };

    struct hor_neighbor : thrust::unary_function<size_t, size_t> {
        const size_t dim_size;
        hor_neighbor(size_t dimension_size): thrust::unary_function<size_t, size_t>(), dim_size(dimension_size){
        }
    };

    struct vert_neighbor : thrust::unary_function<size_t, size_t> {
        const size_t dim_size_x;
        const size_t dim_size_y;
        vert_neighbor(size_t dimension_size_x, size_t dimension_size_y):
        thrust::unary_function<size_t, size_t>(), dim_size_x(dimension_size_x), dim_size_y(dimension_size_y){
        }
    };

    struct left : public hor_neighbor {
        using hor_neighbor::hor_neighbor;
        virtual __host__ __device__ size_t operator()(size_t i) const {
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

            return (i % dim_size == 0) ? i + dim_size - 1 : i - 1;
        }
    };

    struct right : public hor_neighbor {
        using hor_neighbor::hor_neighbor;
        // TODO ist this fine if it is virtual?
        virtual __host__ __device__ size_t operator()(size_t i) const {
            // j is always i+1, expect when i is on the right side of the lattice
            // if i is on the right side of the lattice, j is i - (d - 1)
            // if i is one the right side of the lattice i % lat_dim = lat_dim - 1

            return (i % dim_size == dim_size - 1) ? i - (dim_size - 1) : i + 1;
        }
    };

    struct up : public vert_neighbor{
        using vert_neighbor::vert_neighbor;
        virtual __host__ __device__ size_t operator()(size_t i) const {
            // j is always i - d, except when i is on the upper bound of the lattice
            // if it is on the upper bound, j will be i + d(d-1)
            // if i is on the upper bound, i will be smaller than d

            // okay for rectangular lattices now both dimensions are relevant

            return (i < dim_size_x) ? i + dim_size_x * (dim_size_y - 1) : i - dim_size_x;
        }
    };

    struct down : public vert_neighbor {
        using vert_neighbor::vert_neighbor;
        virtual __host__ __device__ size_t operator()(size_t i) const {
            // j is always i + d, except when i is on the lower bound of the lattice
            // if it is on the lower bound, j will be i - d(d-1)
            // if i is on the lower bound, i will be larger than d * (d-1) - 1 = d*d - d - 1
            return (i >= dim_size_x * (dim_size_y -1 )) ? i - dim_size_x * (dim_size_y - 1) : i + dim_size_x;
        }
    };

    template<class State, class Deriv, class FunctorType>
    void force_calculation(State &x, Deriv &dxdt, double t, FunctorType functor) {
        BOOST_AUTO(start, thrust::make_zip_iterator(thrust::make_tuple(
                x.begin(),
                dxdt.begin() + n,
                thrust::make_permutation_iterator(
                        x.begin(),
                        thrust::make_transform_iterator(
                                thrust::counting_iterator<size_t>(0),
                                left(dim_size_x)           // for left the dim_size in x-direction is important
                        )
                ),
                thrust::make_permutation_iterator(
                        x.begin(),
                        thrust::make_transform_iterator(
                                thrust::counting_iterator<size_t>(0),
                                right(dim_size_x)          // for right the dim_size in x-direction is relevant
                        )
                ),
                thrust::make_permutation_iterator(
                        x.begin(),
                        thrust::make_transform_iterator(
                                thrust::counting_iterator<size_t>(0),
                                up(dim_size_x, dim_size_y)      // for up and down both dimension sizes are relevant
                        )
                ),
                thrust::make_permutation_iterator(
                        x.begin(),
                        thrust::make_transform_iterator(
                                thrust::counting_iterator<size_t>(0),
                                down(dim_size_x, dim_size_y)
                        )
                )
        )));
        thrust::for_each(start, start + n, functor);
    }

    template<class Stoch>
    void calc_diff(Stoch &theta, double t) {
        auto start = thrust::make_zip_iterator(thrust::make_tuple(theta.begin() + n, curand_states.begin()));
        // TODO does it work with start + n?
        thrust::for_each(start, start + n, curand(D));
    }

    template<class State, class Deriv>
    void calc_force(State &x, Deriv &dxdt, double t) {
        cout << "calc force dummy method of system is called" << endl;
    }

    template<class State>
    void map_state(State &x) {
        // Okay this maps the state to a certain interval, doesnt do anything for most systems, but for XY and dipol
        // system we map onto the intervals for the angles
    }

    // depricated constructors, use the one below
    System(map<Parameter, double>& paras) : dim_size_x((size_t)paras[Parameter::dim_size_x]),
                                            dim_size_y((size_t)paras[Parameter::dim_size_y]), eta(paras[Parameter::eta]), T(paras[Parameter::T]),
                                            n((size_t)paras[Parameter::dim_size_x] * (size_t)paras[Parameter::dim_size_y]),
                                            curand_random((bool) paras[Parameter::curand_random]), end_t(paras[Parameter::end_time]){
        D = (sqrt(2 * T * eta));
        cout << "System constructor from Enumeration type Map is called with eta = " << eta << "  T = " << T << endl;
        cout << "Using curand_random: " << curand_random << endl;
        if(curand_random) {
            curand_states = thrust::device_vector<curandState>(2 * n);  // TODO is this okay?
            auto sInit = thrust::make_zip_iterator(thrust::make_tuple(thrust::counting_iterator<int>(0), curand_states.begin()));
            // counting iterator counting from 0 to n. Every lattice site needs its own state i think
            // a state that corresponds to a different sequence. That the sequences are only one apart should not be relevant
            // now we call curand init with sequence numbers from the counting iterator
            chrono::microseconds mus = chrono::duration_cast<chrono::microseconds >(
                    chrono::system_clock::now().time_since_epoch()
            );
            // TODO are you still using this seed? Why is it not possible to get it to work with the default random
            // engine just advancing? Maybe just switch to curand?
            long seed = (mus.count() % 10000000) * 1000000000;
            thrust::for_each(sInit, sInit + 2*n, curand_setup(seed));      // now 2n, moare states should not be a problem and then we can use it also for the bbk random numbers?
        }
    }

    size_t get_dim_size_x() const{
        return dim_size_x;
    }

    virtual size_t get_Lx() const{
        return dim_size_x;
    }

    virtual size_t get_Ly() const{
        return dim_size_y;
    }

    size_t get_dim_size_y() const{
        return dim_size_y;
    }

    double get_cur_T() const{
        // If we Quench, we just have to keep this T up to date and we do that I think, then we can reuse this function
        return T;
    }
    template<class State>
    double calc_m(State &x) {
        double m = thrust::reduce(x.begin(), x.begin() + n, 0.0, thrust::plus<double>()) / (double) n;       // should work?
        return m;
    }

    template<class State>
    double calc_binder(State &x) {
        // No implementation here?
        return 0;
    }

    template<class State>
    void calc_xi(State& x, double& xix, double& xiy) {
        xix = 0;
        xiy = 0;
    }

    template <class T>
    struct sum_square_complex
    {
        __host__ __device__
        double operator()(const T& x, double y) const {
            return y + x.x * x.x + x.y * x.y;
        }
    };

    double get_eta() {
        return eta;
    }

    virtual double get_end_t() const {
        return end_t;
    }

    virtual string get_name() const {
        return "system";
    }

    virtual bool is_equilibrated() {
        return equilibrated;
    }

    virtual void set_equilibration(double t) {
        equilibrated = true;
    }

    virtual double get_quench_time() const {
        // ...dummy method so that we can call that function with observers that would usually need the quench time, but in certain usecases not
        return 0.0;
    }

    ~System() {
        //cudaDeviceSynchronize();
        cout << "System destroyed" << endl;
    }
};

struct System_OBC: virtual public System {
    using System::System;
    struct left : public hor_neighbor {
        using hor_neighbor::hor_neighbor;
        virtual __host__ __device__ size_t operator()(size_t i) const {
            // for open boundary conditions i would say if the first condition is true
            // we just yield i as neighbor of i and therefore no force
            return (i % dim_size == 0) ? i : i - 1;
        }
    };
    struct right : public hor_neighbor {
        using hor_neighbor::hor_neighbor;
        virtual __host__ __device__ size_t operator()(size_t i) const {
            return (i % dim_size == dim_size - 1) ? i : i + 1;
        }
    };
    struct up : public vert_neighbor{
        using vert_neighbor::vert_neighbor;
        virtual __host__ __device__ size_t operator()(size_t i) const {
            return (i < dim_size_x) ? i : i - dim_size_x;
        }
    };
    struct down : public vert_neighbor {
        using vert_neighbor::vert_neighbor;
        virtual __host__ __device__ size_t operator()(size_t i) const {
            return (i >= dim_size_x * (dim_size_y -1 )) ? i : i + dim_size_x;
        }
    };
};

struct quench : virtual public System {
    const double    T_start;       // Start-Temperture of the Quenching For example: 10
    const double    T_end;         // End-Temperature of the Quencheing. For example: 1
    double    s_eq_t;        // start equilibration time: Time the system gets to equilibrate at T_start
    const double    e_eq_t;        // end equilibration time: Time the system gets to equilibrate at T_end (after Quench)
    double          t_quench;            // total time the quench takes, is determined by the quench timescale tau and the temp difference
    double          end_quench_t;        // timepoint at which the Quenching ends = s_eq_t + t_quench
    const double    tau;

    void print_info() override {
        System::print_info();
        cout << "T_start = " << T_start << endl;
        cout << "T_end = " << T_end << endl;
        cout << "s_eq_t = " << s_eq_t << endl;
        cout << "e_eq_t = " << e_eq_t << endl;
        cout << "t_quench = " << t_quench << endl;
        cout << "end_quench_t = " << end_quench_t << endl;
        cout << "tau = " << tau << endl;
    }
    double linear_T(double t) {
        if(s_eq_t < t) {
            // if we are in the quench phase, we reduce T
            System::T = max(T_start - (t - s_eq_t)/tau, T_end);
        }
        return System::T;
    }

    template<class Stoch>
    void calc_diff(Stoch &theta, double t) {
        // I think this should work? We change the D of the System and then just calc the random numbers
        System::D = sqrt(2 * linear_T(t) * System::eta);
        System::calc_diff(theta, t);
        // this is a bit clunky but tbh i don't care to much at this moment anymore
        if(t > get_end_t()) {
            equilibrated = true;
        }
    }

public:
    quench(map<Parameter, double>& paras): System(paras),
                                        tau(paras[Parameter::tau]), T_end(paras[Parameter::end_temp]),
                                        s_eq_t(paras[Parameter::equil_time]),
                                        e_eq_t(paras[Parameter::equil_time_end]),
                                        T_start(paras[Parameter::starting_temp]) {
        // for what are we using end_quench_t again? for get_end_t, right and this is used in simulation.run to start
        // the step until, what do we do if we want to have this quench that dynamically equilibrates, quenches and stops?
        // the step until gets the end_t beforehand, the only other possibility to stop it is through the equilibration attribute
        // we can set s_eq_t through the observer when calling set_equilibration and therefore start the quench
        // after the time of t_quench (which we know beforehand) + s_eq_t + e_eq_t has passed, we return euqilibrated in is
        // equilibrated and therefore stop the stepper.
        if (s_eq_t == 0) {
            // ... we just set s_eq_t to be really large in this case becauce this case means we dynamically equilibrate?
            s_eq_t = 1e10;
        }
        t_quench = (get_quench_time());
        end_quench_t = t_quench + s_eq_t;
    }

    double get_quench_time() const override {
        // returns the time it takes to do the quench
        // in this system, we use a linear quench
        cout << "running get_quench_time:" << endl << "T_start = " << T_start << endl << "T_end = " << T_end << endl << "tau = " << tau << endl << endl;
        return (T_start - T_end) * tau;
    }

    double get_end_t() const override {
        // the total time are the two equilibriate times + the quench time
        return (double)(s_eq_t + e_eq_t + t_quench);
    }

    double get_end_quench_time() const  {
        return t_quench + s_eq_t;
    }

    void set_equilibration(double t) override {
        s_eq_t = t;
    }

    bool is_equilibrated() override{
        // the quench is not supposed to end when the system is equilibrated so we always return false
        return equilibrated;
    }
};

struct quench_nonlinear: virtual public quench {
    int gamma_val;
    void print_info() override {
        quench::print_info();
        cout << "gamma = " << gamma_val << endl;
    }

public:
    double nonlinear_T(double t) {
        if(s_eq_t < t) {
            // if we are in the quench phase, we reduce T
            System::T = max(T_start - pow((t - s_eq_t) / tau, gamma_val), T_end);
        }
        return System::T;
    }

    template<class Stoch>
    void calc_diff(Stoch &theta, double t) {
        // I think this should work? We change the D of the System and then just calc the random numbers
        System::D = sqrt(2 * nonlinear_T(t) * System::eta);
        System::calc_diff(theta, t);
        // this is a bit clunky but tbh i don't care to much at this moment anymore
        if(t > get_end_t()) {
            equilibrated = true;
        }
    }

    quench_nonlinear(map<Parameter, double>& paras): System(paras), quench(paras), gamma_val(paras[Parameter::gamma_exp]) {
        if (s_eq_t == 0) {
            // ... we just set s_eq_t to be really large in this case becauce this case means we dynamically equilibrate?
            s_eq_t = 1e10;
        }
        t_quench = (this->get_quench_time());
        end_quench_t = t_quench + s_eq_t;
    }

    double get_quench_time() const override {
        // returns the time it takes to do the quench
        // in this system, we use a linear quench
        // cout << "running get_quench_time:" << endl << "T_start = " << T_start << endl << "T_end = " << T_end << endl << "tau = " << tau << endl << endl;
        return pow(T_start - T_end, 1.0 / (double)gamma_val) * tau;
    }

};

class XY_model : virtual public System {
public:
    const double p_XY = 2;     // for the potential, will be 2 for bistable XY and 2.5 for silicon
    const double m = 1;     // prefactor for the interaction, important if we dont have the 2 pi periodicy

    double J;
    double h;

    struct xy_functor {
        const double Jx, Jy, h, eta, p_XY, m;

        xy_functor(const double eta, const double J, const double h, const double p_XY, const double m) :
        Jx(J), Jy(J), h(h), eta(eta), p_XY(p_XY), m(m) { }

        xy_functor(const double eta, const double Jx, const double Jy, const double h, const double p_XY, const double m) :
                Jx(Jx), Jy(Jy), h(h), eta(eta), p_XY(p_XY), m(m) { }

        template<class Tup>
        __host__ __device__ void operator()(Tup tup) {
            // Okay I think we have to think about where to % 2pi the system and I think i would like
            // to do it here since I can then easier switch between the models and do not have to adjust the stepper
            double q = thrust::get<0>( tup );
            double p = thrust::get<1>( tup );
            thrust::get<2>( tup ) = p;
            double q_left = thrust::get<4>(tup);
            double q_right = thrust::get<5>(tup);
            double q_up = thrust::get<6>(tup);
            double q_down = thrust::get<7>(tup);

            double interaction =  m * (
                        Jy * (sin(m * (q - q_up))   + sin(m * (q - q_down)))
                    +   Jx * (sin(m * (q - q_left)) + sin(m * (q - q_right)))
            );
            thrust::get<3>( tup ) = (-eta) * p                                                                                  // Friction
                                    + p_XY * h * sin(p_XY * q) // bistable potential
                                    - interaction;       // Interaction
        }
    };

    template <class T>
    struct sin_functor_thrust
    {
        sin_functor_thrust(double p): p(p) {}
        sin_functor_thrust(): p(1.0) {}
        double p;
        __host__ __device__
        T operator()(const T& x) const {
            return sin(p * x);
        }
    };
    struct xy_force {
        const double Jx, Jy, h, eta, p_XY, m;

        xy_force(const double eta, const double J, const double h, const double p_XY, const double m) :
                Jx(J), Jy(J), h(h), eta(eta), p_XY(p_XY), m(m) {
        }
        xy_force(const double eta, const double Jx, const double Jy, const double h, const double p_XY, const double m) :
                Jx(Jx), Jy(Jy), h(h), eta(eta), p_XY(p_XY), m(m) {
        }

        template<class Tup>
        __host__ __device__ void operator()(Tup tup) {
            // Okay I think we have to think about where to % 2pi the system and I think i would like
            // to do it here since I can then easier switch between the models and do not have to adjust the stepper
            double q = thrust::get<0>( tup );
            double q_left = thrust::get<2>(tup);
            double q_right = thrust::get<3>(tup);
            double q_up = thrust::get<4>(tup);
            double q_down = thrust::get<5>(tup);
            double interaction = m * (
                        Jy * (sin(m * (q - q_up))   + sin(m * (q - q_down)))
                    +   Jx * (sin(m * (q - q_left)) + sin(m * (q - q_right)))
            );
            thrust::get<1>( tup ) = p_XY * h * sin(p_XY * q) // bistable potential
                                    - interaction;       // Interaction
        }
    };

    struct map_functor {
        template<class value_type>
        __host__ __device__ value_type operator()(value_type q) {
            return positive_modulo(q, 2.0 * M_PI);
        }
    };

    struct init_functor {
        size_t n, dim_size_x;
        double equil_pos, range_min, range_max;

        init_functor(size_t n, double equil_pos, double range_min, double range_max, size_t dim_size_x): n(n),
        equil_pos(equil_pos), range_min(range_min), range_max(range_max), dim_size_x(dim_size_x) {}
        template<class State>
        void operator()(map<Parameter, double>& paras, State &x) {
            cout << "XY init_state is called with equil_pos = " << equil_pos << endl;
            if(paras[random_init] == 0.0) {
                // equilibrium initialization -> we are in XY model with p=2, meaning we have
                // our equilibria at pi/2 and 3pi/2, we initialize everything in the pi/2 minimum
                thrust::fill(x.begin(), x.begin()+n, equil_pos);
                // The measurments for systems with large h really do not want to leave this state, so we
                // flip every 10th dimer ... does this make a difference for the correlation function?
                // should we flip random dimers?

/*                for(int i = 0; i < n - 10; i += 10) {
                    x[i] *= (-1);   // is this again valid because I am accessing the gpu state?
                }*/
/*                int subsystem_size = (int)(paras[subsystem_Lx] * paras[subsystem_Ly]);
                int nr_random_numbers = subsystem_size / 10; // every 10th?*/
                vector<int> random_numbers = generateRandomIntegers(n / 10, n - 1);
                for(int flip_ind : random_numbers) {
                    x[flip_ind] *= (-1);
                }
/*                for(int subsys_nr = 0; subsys_nr < (int)paras[nr_subsystems]; subsys_nr++) {
                    int ind = subsys_nr * subsystem_size;
                    // ah shit it is again more complicated than this...
                    vector<int> randomNumbers = generateRandomIntegers(nr_random_numbers, subsystem_size);
                }*/

                if(paras[Parameter::J] < 0) {
                    cout << "initializing with chess trafo" << endl;
                    chess_trafo_rectangular(x, dim_size_x);
                }

            } else if (paras[random_init] == 1.0) {
                // random initialization
                double p_ampl = paras[p0];

                chrono::milliseconds ms = chrono::duration_cast<chrono::milliseconds >(
                        chrono::system_clock::now().time_since_epoch()
                );
                auto seed = ms.count() % 10000;
                thrust::counting_iterator<size_t> index_sequence_begin(seed * x.size());
                // theta is uniformly distributed between 0 and 2 pi
                thrust::transform(index_sequence_begin,
                                  index_sequence_begin + n,
                                  x.begin(),
                                  rand_uni_values(range_min, range_max));
                // impuls or angular velocity normal distributed around 0;
                thrust::transform(index_sequence_begin + n,
                                  index_sequence_begin + 2*n,
                                  x.begin() + n,
                                  rand_normal_values(p_ampl, 0, 1));

            } else {
                cout << "I am a bit confused, are the impulses not initalized to be zero?" << endl;
                thrust::fill(x.begin() + n, x.end(), 0.0);
            }
        }
    };

public:

    template<class State>
    void map_state(State &x) {
        // we map every q onto [-pi/2, pi/2]
        thrust::transform(x.begin(), x.begin() + n, x.begin(), map_functor());
    }

    template<class State>
    void init_state(map<Parameter, double>& paras, State &x) {
        init_functor(n, M_PI/2, 0, 2 * M_PI, dim_size_x)(paras, x);
    }

    void print_info() override{
        System::print_info();
        cout << "h = " << h << endl;
        cout << "J = " << J << endl;
    }

    template<class State, class Deriv>
    void calc_force(State &x, Deriv &dxdt, double t) {
        xy_force functor = xy_force(System::eta, J, h, p_XY, m);
        System::force_calculation(x, dxdt, t, functor);
    }

    XY_model(map<Parameter,double>& paras)
    : System(paras), J(paras[Parameter::J]), h(paras[alpha]), m(paras[Parameter::m]), p_XY(paras[Parameter::p]) {
        cout << "Xy model constructor is called "<< endl;
        print_info();
    }

    double get_h() {
        return h;
    }
};

class XY_Silicon: virtual public XY_model {
    // we only have to set the map functor new? And obviously the right p, m stuff
public:
    double p_XY = 2.57;           // does this work? Can i just redeclare them here?
    const double m = 2.0;

    struct map_functor {        // TODO Question if this already works without redefining map_state?
        template<class value_type>
        __host__ __device__ value_type operator()(value_type q) {
            return positive_modulo((q + M_PI/2), M_PI) - M_PI/2;
        }
    };

    template<class State>
    void map_state(State &x) {
        // we map every q onto [-pi/2, pi/2]
        thrust::transform(x.begin(), x.begin() + n, x.begin(), XY_Silicon::map_functor());
    }

    template<class State, class Deriv>
    void calc_force(State &x, Deriv &dxdt, double t) {
        XY_model::xy_force functor = XY_model::xy_force(System::eta, J, h, p_XY, m);
        force_calculation(x, dxdt, t, functor);
    }

    template<class State>
    void init_state(map<Parameter, double>& paras, State &x) {
        double equil_pos = (7.0 / 18.0) * M_PI;
        double range_min = - M_PI / 2;
        double range_max = M_PI / 2;
        init_functor(n, equil_pos, range_min, range_max, dim_size_x)(paras, x);
    }

    XY_Silicon(map<Parameter,double>& paras): System(paras), XY_model(paras) {
        if(paras[Parameter::p] > 0) {
            p_XY = paras[Parameter::p];
        }
    }

    void print_info() override {
        XY_model::print_info();
        cout << "p_XY = " << p_XY << endl;
        cout << "m = " << m << endl;
    }
};

class XY_Silicon_anisotrop: virtual public XY_Silicon {
public:
    const double Jx, Jy;
    XY_Silicon_anisotrop(map<Parameter,double>& paras): System(paras), XY_model(paras), XY_Silicon(paras), Jx(paras[Parameter::J]), Jy(paras[Parameter::Jy]) {}

    template<class State, class Deriv>
    void calc_force(State &x, Deriv &dxdt, double t) {
        XY_model::xy_force functor = XY_model::xy_force(System::eta, Jx, Jy, h, p_XY, m);
        force_calculation(x, dxdt, t, functor);
    }

};

class subsystems: virtual public System {
    // This one is supposed to simulate multiple small (isolated) systems in one run
    // the only thing we have to change is how the neighbors are selected?
protected:
    size_t n0;      // subsystem size
    const size_t Lx, Ly;  // subsystem dimensions


public:

    struct subsystem_running_index : thrust::unary_function<size_t, size_t> {
        const size_t dim_size_x, Lx, Ly;
        subsystem_running_index(size_t dimension_size_x, size_t Lx, size_t Ly):
                thrust::unary_function<size_t, size_t>(), dim_size_x(dimension_size_x), Lx(Lx), Ly(Ly){
        }
        __host__ __device__ size_t operator()(size_t i) const {
            // first we need to know in which system we are
            size_t sys_nr = i  /  (Lx * Ly);
            size_t pos_in_subsystem = i % (Lx * Ly);
            size_t row = pos_in_subsystem / Lx;
            size_t col = pos_in_subsystem % Lx;
            size_t index_in_whole_system = sys_nr * Lx + row * dim_size_x + col;
            return index_in_whole_system;
        }
    };

    struct running_chess_trafo_iterator : thrust::unary_function<thrust::tuple<double, size_t>, double> {
        size_t dim_size_x, Lx;
        running_chess_trafo_iterator(size_t dimension_size_x, size_t Lx):
                thrust::unary_function<thrust::tuple<double, size_t>, double>(),
                dim_size_x(dimension_size_x), Lx(Lx){
        }
        template<class Tup>
        __host__ __device__ double operator()(Tup tup) const {
            // first we need to know in which system we are
            if (dim_size_x % 2 == 1) {
                // if it is uneven we can just trafo every even
                if(thrust::get<1>(tup) % 2 == 0) {
                    return (-1) * thrust::get<0>(tup);
                }
            } else {
                // else we trafo in even rows even indices and in uneven rows
                // uneven indices
                int ind = thrust::get<1>(tup);
                int row_even = (ind / Lx) % 2;
                if(ind % 2 - row_even == 0) {
                    return (-1) * thrust::get<0>(tup) ;
                }
            }
            return  thrust::get<0>(tup);
        }
    };

    struct rearrange_ft_inds : thrust::unary_function<int, int> {
        int Ly, Lx;
        rearrange_ft_inds(int Lx, int Ly): Ly(Ly), Lx(Lx){
        }

        __host__ __device__ double operator()(int ind) const {
            int row = ind % Ly;
            int col = ind / Ly;
            return  row * Lx + col;
        }
    };

    subsystems(map<Parameter, double> paras): System(paras), Lx((size_t)paras[subsystem_Lx]), Ly((size_t)paras[subsystem_Ly]) {
        n0 = Ly * Lx;
    }

    size_t get_Lx() const override {
        return Lx;
    }

    size_t get_Ly() const override {
        return Ly;
    }

};

class subsystems_pbc : virtual public System, virtual public subsystems {
public:
    template<class State, class Deriv, class FunctorType>
    void force_calculation(State &x, Deriv &dxdt, double t, FunctorType functor) {
        BOOST_AUTO(start, thrust::make_zip_iterator(thrust::make_tuple(
                x.begin(),
                dxdt.begin() + n,
                thrust::make_permutation_iterator(
                        x.begin(),
                        thrust::make_transform_iterator(
                                thrust::counting_iterator<size_t>(0),
                                left(Lx)           // for left the dim_size in x-direction is important
                        )
                ),
                thrust::make_permutation_iterator(
                        x.begin(),
                        thrust::make_transform_iterator(
                                thrust::counting_iterator<size_t>(0),
                                right(Lx)          // for right the dim_size in x-direction is relevant
                        )
                ),
                thrust::make_permutation_iterator(
                        x.begin(),
                        thrust::make_transform_iterator(
                                thrust::counting_iterator<size_t>(0),
                                System::up(dim_size_x, Ly)      // for up and down both dimension sizes are relevant
                        )
                ),
                thrust::make_permutation_iterator(
                        x.begin(),
                        thrust::make_transform_iterator(
                                thrust::counting_iterator<size_t>(0),
                                System::down(dim_size_x, Ly)
                        )
                )
        )));
/*        auto left_inds = thrust::make_transform_iterator(thrust::counting_iterator<size_t>(0), left(Lx));
        auto right_inds = thrust::make_transform_iterator(thrust::counting_iterator<size_t>(0), right(Lx));
        auto up_inds = thrust::make_transform_iterator(thrust::counting_iterator<size_t>(0), System::up(dim_size_x, Ly));
        auto down_inds = thrust::make_transform_iterator(thrust::counting_iterator<size_t>(0), System::down(dim_size_x, Ly));
        for(int i = 0; i < dim_size_x * dim_size_y; i++) {
            cout << "pos: " << i << "  left: " << left_inds[i] << "  right: " << right_inds[i] << "  up: " << up_inds[i] << "  down: " << down_inds[i] << endl;
        }*/
        thrust::for_each(start, start + n, functor);
    }

    subsystems_pbc(map<Parameter, double> paras): System(paras), subsystems(paras) {
        cout << "initializing subsystems pbc with sizes: " << Lx << " x " << Ly << " = " << n0 << endl;
        if (Ly != dim_size_y) {
            cout << "dim_size_y has to be equal to L_y" << endl;
            exit(0);
        }
    }
};

class subsystems_obc : virtual public System_OBC, virtual public subsystems {
public:
    template<class State, class Deriv, class FunctorType>
    void force_calculation(State &x, Deriv &dxdt, double t, FunctorType functor) {
        BOOST_AUTO(start, thrust::make_zip_iterator(thrust::make_tuple(
                x.begin(),
                dxdt.begin() + n,
                thrust::make_permutation_iterator(
                        x.begin(),
                        thrust::make_transform_iterator(
                                thrust::counting_iterator<size_t>(0),
                                System_OBC::left(Lx)           // for left the dim_size in x-direction is important
                        )
                ),
                thrust::make_permutation_iterator(
                        x.begin(),
                        thrust::make_transform_iterator(
                                thrust::counting_iterator<size_t>(0),
                                System_OBC::right(Lx)          // for right the dim_size in x-direction is relevant
                        )
                ),
                thrust::make_permutation_iterator(
                        x.begin(),
                        thrust::make_transform_iterator(
                                thrust::counting_iterator<size_t>(0),
                                System_OBC::up(dim_size_x, Ly)      // for up and down both dimension sizes are relevant
                        )
                ),
                thrust::make_permutation_iterator(
                        x.begin(),
                        thrust::make_transform_iterator(
                                thrust::counting_iterator<size_t>(0),
                                System_OBC::down(dim_size_x, Ly)
                        )
                )
        )));
/*        auto left_inds = thrust::make_transform_iterator(thrust::counting_iterator<size_t>(0), left(Lx));
        auto right_inds = thrust::make_transform_iterator(thrust::counting_iterator<size_t>(0), right(Lx));
        auto up_inds = thrust::make_transform_iterator(thrust::counting_iterator<size_t>(0), System::up(dim_size_x, Ly));
        auto down_inds = thrust::make_transform_iterator(thrust::counting_iterator<size_t>(0), System::down(dim_size_x, Ly));
        for(int i = 0; i < dim_size_x * dim_size_y; i++) {
            cout << "pos: " << i << "  left: " << left_inds[i] << "  right: " << right_inds[i] << "  up: " << up_inds[i] << "  down: " << down_inds[i] << endl;
        }*/
        thrust::for_each(start, start + n, functor);
    }

    subsystems_obc(map<Parameter, double> paras): System_OBC(paras), System(paras), subsystems(paras) {
        cout << "initializing subsystems obc with sizes: " << Lx << " x " << Ly << " = " << n0 << endl;
        if (Ly != dim_size_y) {
            cout << "dim_size_y has to be equal to L_y" << endl;
            exit(0);
        }
    }
};

struct XY_silicon_anisotrop_subsystems : public subsystems_pbc, public XY_Silicon_anisotrop {
public:
    XY_silicon_anisotrop_subsystems(map<Parameter, double> paras) : XY_model(paras),
                                                            subsystems(paras),
                                                            subsystems_pbc(paras),
                                                            XY_Silicon(paras),
                                                            XY_Silicon_anisotrop(paras),
                                                            System(paras) {}
    template<class State, class Deriv>
    void calc_force(State &x, Deriv &dxdt, double t) {
        xy_force functor = xy_force(System::eta, Jx, Jy, h, XY_Silicon::p_XY, XY_Silicon::m);
        subsystems_pbc::force_calculation(x, dxdt, t, functor);
    }

    struct segment_functor {
        size_t seg_len;
        segment_functor(size_t seg_len): seg_len(seg_len){}
        __host__ __device__ int operator()(int ind) const {
            return ind / (seg_len);
        }
    };

    struct avg_functor {
        double Lx, Ly;
        avg_functor(double Lx, double Ly): Lx(Lx), Ly(Ly){}
        __host__ __device__ double operator()(double m) const {
            return m / (Lx * Ly);
        }
    };

    template<class State>
    double calc_binder(State& x) {
        // We return the average over the subsystems we have?
        // so we need to extract the subsystems and calc m, average at the end
        // actually we should implement this for the normal XY model and make it useable for everything
        // ahh difficult becauso polymorphism doenst work and we inherit from two classes
        // hmmm decisions...
        // we also cannot really reuse the methods from our lattice calculations because we want to use
        // gpu methods
        // i think for now we will only implement it for this class...
        bool gpu = true;
        int nr_subsystems = dim_size_x / Lx;
        thrust::host_vector<double> m_vec{};
        // damn is this already it for the cell?
        auto cell = thrust::make_permutation_iterator(
                x.begin(),
                thrust::make_transform_iterator(
                        thrust::counting_iterator<size_t>(0),
                        subsystem_running_index(dim_size_x, Lx, Ly)
                )
        );
/*        auto subsystem_inds = thrust::make_transform_iterator(
                thrust::counting_iterator<size_t>(0),
                subsystem_running_index(dim_size_x, Lx, Ly)
        );*/
/*        for (int j = 0; j < dim_size_x * Ly; j++) {
            cout << "pos in sub: " << j << " pos in whole: " << subsystem_inds[j] << " value:" << cell[j] << endl;
        }*/
        // doing chess trafo
        auto tuple = thrust::make_zip_iterator(thrust::make_tuple(cell, thrust::counting_iterator<size_t>(0)));
        //BOOST_AUTO(tuple, thrust::make_zip_iterator(thrust::make_tuple(cell, thrust::counting_iterator<size_t>(0))));
        /*thrust::for_each(tuple, tuple + (dim_size_x * Ly), running_chess_trafo(dim_size_x, Lx));*/
        auto cell_trafo = thrust::make_transform_iterator(tuple, running_chess_trafo_iterator(dim_size_x, Lx));
/*        for (int j = 0; j < dim_size_x * Ly; j++) {
            cout << "pos in sys: " << j  << " value: " << cell_trafo[j] << endl;
        }
        for (int j = 0; j < dim_size_x * Ly; j++) {
            if (j % dim_size_x == 0) {
                cout << endl;
            }
            cout << cell_trafo[j] << ", ";

        }*/

        if (gpu) {
            // before we used reduce by key, we extracted the sinus of the cell by using transform reduce,
            // now we should transform beforehand i guess
            auto sin_cell = thrust::make_transform_iterator(cell_trafo, sin_functor_thrust<double>(XY_Silicon::p_XY / 2.0));
            thrust::device_vector<double> m_vec_gpu(nr_subsystems);
            //auto segment_functor = [this] __device__ (int ind) {return ind / (Lx * Ly);};
            auto segment_keys = thrust::make_transform_iterator(thrust::counting_iterator<int>(0), segment_functor(Lx * Ly));
/*            for (int j = 0; j < dim_size_x * Ly; j++) {
                if (j % dim_size_x == 0) {
                    cout << endl;
                }
                cout << segment_keys[j] << ", ";
            }*/
            thrust::reduce_by_key(segment_keys,        // is it slow because of this this?
                                  segment_keys + (nr_subsystems * Lx * Ly),
                                  sin_cell,
                                  thrust::make_discard_iterator(),
                                  m_vec_gpu.begin());

            // we would now have to copy the m values back to the host to advance with the rest of the code
            // other possibility is to just perform the rest on the GPU aswell
            // I am not sure what will timewise be the better option, I guess simpler would be copying for now
            // we need to average the U_Ls
            //auto avg_functor = [this] __device__ (double m) {return m/(double)(Lx * Ly);};
/*            cout << endl << endl << "m_vec_gpu [";
            for(int i = 0; i < nr_subsystems; i++) {
                cout << m_vec_gpu[i] << ", ";
            }
            cout << "]" << endl;*/ // from those values it seems he does everything correct so the error has to be earlier
            thrust::transform(m_vec_gpu.begin(), m_vec_gpu.end(), m_vec_gpu.begin(), avg_functor((double)Lx, (double)Ly));
            m_vec = thrust::host_vector<double>(m_vec_gpu);
        } else {
            for(int i = 0; i < nr_subsystems; i++) {
                // for every subsystem we need to extract it
                double m;
                if (gpu) {
                    m = thrust::transform_reduce(cell_trafo + i * (Lx * Ly), cell_trafo + (i+1) * (Lx * Ly),
                                                 sin_functor_thrust<double>(XY_Silicon::p_XY /  2.0), 0.0, thrust::plus<double>()) / ((double) (Lx * Ly));
                    // cout << "m =" << m << endl;
                } else {
                    vector<double> cell(Lx * Ly);
                    extract_cell(x, cell, i, Lx, Ly, dim_size_x);
                    // could work right? :)
                    // chess trafo
                    chess_trafo_rectangular(cell, Lx);
                    // calc m
                    m = transform_reduce(cell.begin(), cell.end(), 0.0, plus<double>(), sin_functor(XY_Silicon::p_XY /  2.0)) / ((double) (Lx * Ly));
                }
                m_vec.push_back(m);
            }
        }

        double m_L2 = std::transform_reduce(m_vec.begin(), m_vec.end(),
                                            0.0, // initial value for the reduction (sum)
                                            std::plus<double>(), // transformation (square)
                                            [](double m) -> double { return m * m; });
        // cout << "Dividing by " << m_vec.size() << endl;
        // exit(0);
        m_L2 /= (double)m_vec.size();
        double m_L4 = std::transform_reduce(m_vec.begin(), m_vec.end(),
                                            0.0, // initial value for the reduction (sum)
                                            std::plus<>(), // transformation (square)
                                            [](double m) { return (pow(m, 4)); });
        m_L4 /= (double)m_vec.size();

        double cum = m_L4 / (m_L2 * m_L2);
        // cout << "cum = " << cum << endl << endl;
        return cum;
    }

    template<class State>
    double calc_m(State& x) {
        vector<double> m_vec = calc_m_vec(x);
        double m = std::reduce(m_vec.begin(), m_vec.end()) / (double)m_vec.size();

        return m;
    }

    template<class State>
    thrust::host_vector<double> calc_m_vec(State& x) {
        // We return the average over the subsystems we have?
        // so we need to extract the subsystems and calc m, average at the end
        // actually we should implement this for the normal XY model and make it useable for everything
        // ahh difficult becauso polymorphism doenst work and we inherit from two classes
        // hmmm decisions...
        // we also cannot really reuse the methods from our lattice calculations because we want to use
        // gpu methods
        // i think for now we will only implement it for this class...
        bool gpu = true;
        int nr_subsystems = dim_size_x / Lx;
        thrust::host_vector<double> m_vec{};
        // damn is this already it for the cell?
        auto cell = thrust::make_permutation_iterator(
                x.begin(),
                thrust::make_transform_iterator(
                        thrust::counting_iterator<size_t>(0),
                        subsystem_running_index(dim_size_x, Lx, Ly)
                )
        );

        // doing chess trafo
        auto tuple = thrust::make_zip_iterator(thrust::make_tuple(cell, thrust::counting_iterator<size_t>(0)));
        //BOOST_AUTO(tuple, thrust::make_zip_iterator(thrust::make_tuple(cell, thrust::counting_iterator<size_t>(0))));
        auto cell_trafo = thrust::make_transform_iterator(tuple, running_chess_trafo_iterator(dim_size_x, Lx));

        if (gpu) {
            // before we used reduce by key, we extracted the sinus of the cell by using transform reduce,
            // now we should transform beforehand i guess
            auto sin_cell = thrust::make_transform_iterator(cell_trafo, sin_functor_thrust<double>(XY_Silicon::p_XY / 2.0));
            thrust::device_vector<double> m_vec_gpu(nr_subsystems);
            //auto segment_functor = [this] __device__ (int ind) {return ind / (Lx * Ly);};
            auto segment_keys = thrust::make_transform_iterator(thrust::counting_iterator<int>(0), segment_functor(Lx * Ly));
            thrust::reduce_by_key(segment_keys,        // is it slow because of this this?
                                  segment_keys + (nr_subsystems * Lx * Ly),
                                  sin_cell,
                                  thrust::make_discard_iterator(),
                                  m_vec_gpu.begin());

            thrust::transform(m_vec_gpu.begin(), m_vec_gpu.end(), m_vec_gpu.begin(), avg_functor((double)Lx, (double)Ly));
            m_vec = thrust::host_vector<double>(m_vec_gpu);
        } else {
            for(int i = 0; i < nr_subsystems; i++) {
                // for every subsystem we need to extract it
                double m;
                if (gpu) {
                    m = thrust::transform_reduce(cell_trafo + i * (Lx * Ly), cell_trafo + (i+1) * (Lx * Ly),
                                                 sin_functor_thrust<double>(XY_Silicon::p_XY /  2.0), 0.0, thrust::plus<double>()) / ((double) (Lx * Ly));
                    // cout << "m =" << m << endl;
                } else {
                    vector<double> cell(Lx * Ly);
                    extract_cell(x, cell, i, Lx, Ly, dim_size_x);
                    // could work right? :)
                    // chess trafo
                    chess_trafo_rectangular(cell, Lx);
                    // calc m
                    m = transform_reduce(cell.begin(), cell.end(), 0.0, plus<double>(), sin_functor(XY_Silicon::p_XY /  2.0)) / ((double) (Lx * Ly));
                }
                m_vec.push_back(m);
            }
        }

        return m_vec;
    }

    template<class State>
    double calc_chi(State& x) {
        thrust::host_vector<double> m_vec = calc_m_vec(x);
        double m_L2 = std::transform_reduce(m_vec.begin(), m_vec.end(),
                                            0.0, // initial value for the reduction (sum)
                                            std::plus<double>(),
                                            [](double m_val) -> double { return m_val * m_val; }) /
                                                    (double)(m_vec.size() * m_vec.size());

        return m_L2;
    }

    struct cos_real_to_complex {
        template<class Tup>
        __host__ __device__ void operator()(Tup tup) const {
            thrust::get<0>(tup).x = cos(thrust::get<1>(tup));
            thrust::get<0>(tup).y = 0.0;
        }
    };

    struct sin_real_to_complex {
        sin_real_to_complex(double p): p(p) {}
        sin_real_to_complex(): p(1.0) {}
        double p;
        template<class Tup>
        __host__ __device__ void operator()(Tup tup) const {
            thrust::get<0>(tup).x = sin(p * thrust::get<1>(tup));
            thrust::get<0>(tup).y = 0.0;
        }
    };

    template<class State>
    void calc_xi(State& x, double& xix, double& xiy) {
        thrust::host_vector<double> ft_squared_k;
        thrust::host_vector<double> ft_squared_l;
        get_ft_vectors(x, ft_squared_k, ft_squared_l);
        fit_ft(ft_squared_k, ft_squared_l, xix, xiy);
    }

    void fit_ft(thrust::host_vector<double> &ft_squared_k, thrust::host_vector<double> &ft_squared_l, double& xix, double& xiy) {
        bool cut_zero_impuls = true;
        bool cut_around_peak = true;
        double *ft_k_fit;
        double *ft_l_fit;

        vector<double> kx = get_frequencies_fftw_order(Lx);
        vector<double> ky = get_frequencies_fftw_order(Ly);
        if(cut_zero_impuls) {
            // since we now know how pointers work, we know that we can just increase the value of ft_squared_k
            // which is the memory adress of the first value in the array by 1 to increase the adress value by the size of a double
            // this way we have a pointer that points to the second value of ft_squared_k
            ft_k_fit = &ft_squared_k[1];
            ft_l_fit = &ft_squared_l[1];
            // remove first value form ks
            kx.erase(kx.begin());
            ky.erase(ky.begin());
            // size for fitting
        } else {
            ft_k_fit = &ft_squared_k[0];
            ft_l_fit = &ft_squared_l[0];
        }
/*        cout << "ft_l_fit BEFORE cut around peak:" << endl;
        print_array(ft_l_fit, ky.size());
        cout << "" << endl;
        cout << "ky BEFORE cut around peak:" << endl;
        print_vector(ky);
        cout << "" << endl;*/
        if(cut_around_peak) {
            auto cut_data_k = cut_data_around_peak(kx, ft_k_fit);   // It doesnt need more computation power if i asign it to a pair beforehand or?
            kx = cut_data_k.first;
            ft_k_fit = cut_data_k.second;

            auto cut_data_l = cut_data_around_peak(ky, ft_l_fit);   // It doesnt need more computation power if i asign it to a pair beforehand or?
            ky = cut_data_l.first;
            ft_l_fit = cut_data_l.second;
        }

        // TODO improve the fit, have a thingy that makes it grounded and maybe only fit the peak ? For now ok
// We now call another function if we want to cut with offset, this is tedious to change but I dont know if
// we ever will so I won't bother for now
/*        cout << "ft_l_fit after cut around peak:" << endl;
        print_array(ft_l_fit, ky.size());
        cout << "ky AFTER cut around peak:" << endl;
        print_vector(ky);*/
        Singleton_timer::set_endpoint("Xi summing and averaging");
        Singleton_timer::set_startpoint("Xi Fitting");
        Eigen::VectorXd paras_x = fit_offset_lorentz_peak(kx, ft_k_fit);
        Eigen::VectorXd paras_y = fit_offset_lorentz_peak(ky, ft_l_fit);
        Singleton_timer::set_endpoint("Xi Fitting");
        xix = abs(paras_x(1));
        xiy = abs(paras_y(1));

        Singleton_timer::set_endpoint("Xi calculation Evaluation of FFTs");
    }

    template<class State>
    void get_ft_vectors(const State &x, thrust::host_vector<double> &ft_squared_k,
                   thrust::host_vector<double> &ft_squared_l) const {
        cufftHandle plan;
        Singleton_timer::set_startpoint("Xi calculation FFTs");
        using dim_t = array<int, 2>;
        dim_t fft_size = {(int) Ly, (int) Lx};
        // data is x
        int nr_batches = (int) (dim_size_x / Lx);     // so if I understand correctly the batch size is the number of multiple
        int batch_size = fft_size[0] * fft_size[1];
        // ffts running at the same time. so since I want to do a fft for every subsystem, my batch size will
// be the number of subsystems?

        // sum vector for alter, this will be the squared ft
        thrust::device_vector<double> sum(batch_size);
        cufftCreate(&plan);
        cufftPlanMany(&plan, fft_size.size(), fft_size.data(),
                      nullptr, 1, 0,
                      nullptr, 1, 0,
                      CUFFT_C2C, nr_batches);

        auto cell = thrust::make_permutation_iterator(
                x.begin(),
                thrust::make_transform_iterator(
                        thrust::counting_iterator<size_t>(0),
                        subsystem_running_index(dim_size_x, Lx, Ly)
                )
        );
        BOOST_AUTO(tuple, thrust::make_zip_iterator(thrust::make_tuple(cell, thrust::counting_iterator<size_t>(0))));
        auto cell_trafo = thrust::make_transform_iterator(tuple, running_chess_trafo_iterator(dim_size_x, Lx)); // so cell trafo is the transformed running_chess_it(cell, ind) which is a double again

        // from cell trafo we need to somehow create a pointer to complex<double>
// workaraound -> device vector erstellen
        thrust::device_vector<cufftComplex> input_vector(dim_size_x * Ly), output_vector(dim_size_x * Ly);

        // fill input vector with the real part being cos of theta
        BOOST_AUTO(input_tuple, thrust::make_zip_iterator(thrust::make_tuple(input_vector.begin(), cell_trafo)));
        thrust::for_each(input_tuple, input_tuple + dim_size_x * Ly, cos_real_to_complex());

        // execute FFT
        cufftExecC2C(plan, thrust::raw_pointer_cast(input_vector.data()),
                     thrust::raw_pointer_cast(output_vector.data()), CUFFT_FORWARD);
        cudaDeviceSynchronize();

        //thrust::reduce_by_key(segment_keys, segment_keys + nr_batches * batch_size,
//                      output_vector.begin(), thrust::make_discard_iterator(), sum.begin())
// reduce by key doesnt work here as easy as expected. also since the number of batches is usally small
// for the correlation length calculation, this wont we the problem?
        for(int i = 0; i < nr_batches; i++) {
            thrust::transform(output_vector.begin() + i * batch_size,
                              output_vector.begin() + (i + 1) * batch_size,
                              sum.begin(), sum.begin(), sum_square_complex<cufftComplex>());
        }

        // Now everything again for the sin part?
        input_tuple = thrust::make_zip_iterator(thrust::make_tuple(input_vector.begin(), cell_trafo));
        thrust::for_each(input_tuple, input_tuple + dim_size_x * Ly, sin_real_to_complex());    // fill input with sin

        // execute fft
        cufftExecC2C(plan, thrust::raw_pointer_cast(input_vector.data()),
                     thrust::raw_pointer_cast(output_vector.data()), CUFFT_FORWARD); // fft into output vector
        cudaDeviceSynchronize();

        // sum up
        for(int i = 0; i < nr_batches; i++) {
            thrust::transform(output_vector.begin() + i * batch_size, output_vector.begin() + (i + 1) * batch_size,
                              sum.begin(), sum.begin(), sum_square_complex<cufftComplex>());
        }
        Singleton_timer::set_endpoint("Xi calculation FFTs");
        Singleton_timer::set_startpoint("Xi calculation Evaluation of FFTs");
        Singleton_timer::set_startpoint("Xi summing and averaging");
        thrust::transform(sum.begin(), sum.end(), sum.begin(), thrustDivideBy((double)nr_batches));
        thrust::device_vector<double> ft_squared_k_gpu(Lx);
        thrust::device_vector<double> ft_squared_l_gpu(Ly);

        auto segment_keys_l_sum = thrust::make_transform_iterator(thrust::counting_iterator<int>(0), segment_functor(Lx));     // segments i am summing are Ly long
        auto segment_keys_k_sum = thrust::make_transform_iterator(thrust::counting_iterator<int>(0), segment_functor(Ly));     // segments j are Ly long


        thrust::reduce_by_key(segment_keys_l_sum, segment_keys_l_sum + batch_size,
                              sum.begin(), thrust::make_discard_iterator(), ft_squared_l_gpu.begin());

        auto rearranged_inds = thrust::make_transform_iterator(thrust::make_counting_iterator<int>(0), rearrange_ft_inds(
                Lx, Ly));
        auto rearranged_sum = thrust::make_permutation_iterator(sum.begin(), rearranged_inds);


        thrust::reduce_by_key(segment_keys_k_sum, segment_keys_k_sum + batch_size,
                              rearranged_sum, thrust::make_discard_iterator(), ft_squared_k_gpu.begin());


        // now the entries have to be averaged
        thrust::transform(ft_squared_k_gpu.begin(), ft_squared_k_gpu.end(), ft_squared_k_gpu.begin(), thrustDivideBy((double)pow(
                Ly, 4)));
        thrust::transform(ft_squared_l_gpu.begin(), ft_squared_l_gpu.end(), ft_squared_l_gpu.begin(), thrustDivideBy((double)pow(
                Lx, 4)));// cut zero impuls, I think we can always do that and it won't be much of a difference?
        cufftDestroy(plan);

        ft_squared_k = thrust::host_vector<double>(ft_squared_k_gpu);
        ft_squared_l = thrust::host_vector<double>(ft_squared_l_gpu);
    }

    template<class State>
    void calc_ft(State& x, double* ft_k, double* ft_l) {
        using dim_t = std::array<int, 2>;
        dim_t fft_size = {(int)Ly, (int)Lx};
        // data is x
        int nr_batches = (int) (dim_size_x / Lx);     // so if I understand correctly the batch size is the number of multiple
        int batch_size = fft_size[0] * fft_size[1];
        // ffts running at the same time. so since I want to do a fft for every subsystem, my batch size will
        // be the number of subsystems?

        // sum vector for alter, this will be the squared ft
        thrust::device_vector<double> sum(Lx * Ly);

        cufftHandle plan;
        cufftCreate(&plan);
        cufftPlanMany(&plan, fft_size.size(), fft_size.data(),
                      nullptr, 1, 0,
                      nullptr, 1, 0,
                      CUFFT_C2C, nr_batches);

        auto cell = thrust::make_permutation_iterator(
                x.begin(),
                thrust::make_transform_iterator(
                        thrust::counting_iterator<size_t>(0),
                        subsystem_running_index(dim_size_x, Lx, Ly)
                )
        );
        BOOST_AUTO(tuple, thrust::make_zip_iterator(thrust::make_tuple(cell, thrust::counting_iterator<size_t>(0))));
        auto cell_trafo = thrust::make_transform_iterator(tuple, running_chess_trafo_iterator(dim_size_x, Lx));

        // from cell trafo we need to somehow create a pointer to complex<double>
        // workaraound -> device vector erstellen
        thrust::device_vector<cufftComplex> input_vector(dim_size_x * Ly), output_vector(dim_size_x * Ly);

        // fill input vector with the real part being cos of theta
        BOOST_AUTO(input_tuple, thrust::make_zip_iterator(thrust::make_tuple(input_vector.begin(), cell_trafo)));
        thrust::for_each(input_tuple, input_tuple + dim_size_x * Ly, cos_real_to_complex());

        // execute FFT
        cufftExecC2C(plan, thrust::raw_pointer_cast(input_vector.data()),
                     thrust::raw_pointer_cast(output_vector.data()), CUFFT_FORWARD);
        cudaDeviceSynchronize();

        for(int i = 0; i < nr_batches; i++) {
            thrust::transform(output_vector.begin() + i * batch_size,
                              output_vector.begin() + (i + 1) * batch_size,
                              sum.begin(), sum.begin(), sum_square_complex<cufftComplex>());
        }

        // Now everything again for the sin part?
        input_tuple = thrust::make_zip_iterator(thrust::make_tuple(input_vector.begin(), cell_trafo));
        thrust::for_each(input_tuple, input_tuple + dim_size_x * Ly, sin_real_to_complex());    // fill input with sin

        // execute fft
        cufftExecC2C(plan, thrust::raw_pointer_cast(input_vector.data()),
                     thrust::raw_pointer_cast(output_vector.data()), CUFFT_FORWARD); // fft into output vector
        cudaDeviceSynchronize();

        // sum up
        for(int i = 0; i < nr_batches; i++) {
            thrust::transform(output_vector.begin() + i * batch_size, output_vector.begin() + (i + 1) * batch_size,
                              sum.begin(), sum.begin(), sum_square_complex<cufftComplex>());
        }

        thrust::transform(sum.begin(), sum.end(), sum.begin(), thrustDivideBy((double)nr_batches));
        thrust::device_vector<double> ft_squared_k_gpu(Lx);
        thrust::device_vector<double> ft_squared_l_gpu(Ly);

        auto segment_keys_l_sum = thrust::make_transform_iterator(thrust::counting_iterator<int>(0), segment_functor(Lx));     // segments i am summing are Ly long
        auto segment_keys_k_sum = thrust::make_transform_iterator(thrust::counting_iterator<int>(0), segment_functor(Ly));     // segments j are Ly long


        thrust::reduce_by_key(segment_keys_l_sum, segment_keys_l_sum + batch_size,
                              sum.begin(), thrust::make_discard_iterator(), ft_squared_l_gpu.begin());

        auto rearranged_inds = thrust::make_transform_iterator(thrust::make_counting_iterator<int>(0), rearrange_ft_inds(
                Lx, Ly));
        auto rearranged_sum = thrust::make_permutation_iterator(sum.begin(), rearranged_inds);


        thrust::reduce_by_key(segment_keys_k_sum, segment_keys_k_sum + batch_size,
                              rearranged_sum, thrust::make_discard_iterator(), ft_squared_k_gpu.begin());


        // now the entries have to be averaged
        thrust::transform(ft_squared_k_gpu.begin(), ft_squared_k_gpu.end(), ft_squared_k_gpu.begin(), thrustDivideBy((double)pow(
                Ly, 4)));
        thrust::transform(ft_squared_l_gpu.begin(), ft_squared_l_gpu.end(), ft_squared_l_gpu.begin(), thrustDivideBy((double)pow(
                Lx, 4)));// cut zero impuls, I think we can always do that and it won't be much of a difference?
        cufftDestroy(plan);

        auto ft_squared_k = thrust::host_vector<double>(ft_squared_k_gpu);
        auto ft_squared_l = thrust::host_vector<double>(ft_squared_l_gpu);

        // TODO oh boy does this even work? why were you working with arrays anyway?
        ft_k = &ft_squared_k[0];
        ft_l = &ft_squared_l[0];

        cufftDestroy(plan);
    }

    string get_name() const override {
        return "XY_silicon_anisotrop_subsystems";
    }
};

struct XY_silicon_anisotrop_subsystems_obc : public subsystems_obc, virtual public XY_silicon_anisotrop_subsystems{
public:
    XY_silicon_anisotrop_subsystems_obc(map<Parameter, double> paras) :
                                                                    XY_silicon_anisotrop_subsystems(paras),
                                                                    XY_model(paras),
                                                                    subsystems(paras),
                                                                    subsystems_obc(paras),
                                                                    XY_Silicon(paras),
                                                                    System_OBC(paras),
                                                                    System(paras){}

    template<class State, class Deriv>
    void calc_force(State &x, Deriv &dxdt, double t) {
        xy_force functor = xy_force(System::eta, Jx, Jy, h, XY_Silicon::p_XY, XY_Silicon::m);
        subsystems_obc::force_calculation(x, dxdt, t, functor);
    }

    string get_name() const override {
        return "XY_silicon_anisotrop_subsystems_obc";
    }
};

struct XY_silicon_anisotrop_subsystems_quench : public XY_silicon_anisotrop_subsystems, public quench {
public:
    void print_info() override {
        XY_silicon_anisotrop_subsystems::print_info();
        quench::print_info();
    }

    XY_silicon_anisotrop_subsystems_quench(map<Parameter, double> paras) : quench(paras),
    XY_silicon_anisotrop_subsystems(paras),
    XY_Silicon(paras),
    XY_model(paras),
    subsystems(paras),
    System(paras){
        cout << "XY_silicon_anisotrop_subsystems_quench system constructed";
    }

    string get_name() const override {
        return "XY_silicon_anisotrop_subsystems_quench";
    }
};

struct XY_silicon_anisotrop_subsystems_nonlinear_quench: public XY_silicon_anisotrop_subsystems, public quench_nonlinear{
public:
    void print_info() override {
        XY_silicon_anisotrop_subsystems::print_info();
        quench_nonlinear::print_info();
    }

    XY_silicon_anisotrop_subsystems_nonlinear_quench(map<Parameter, double> paras):
            quench(paras),
            quench_nonlinear(paras),
            XY_silicon_anisotrop_subsystems(paras),
            XY_Silicon(paras),
            XY_model(paras),
            subsystems(paras),
            System(paras){
        cout << "XY_silicon_anisotrop_subsystems_nonlinear_quench system constructed";
    }

    string get_name() const override {
        return "XY_silicon_anisotrop_subsystems_nonlinear_quench";
    }
};

struct XY_silicon_anisotrop_subsystems_quench_h : public XY_silicon_anisotrop_subsystems, public quench {
    double h_start;
    double h_end;

    void linear_h(double t) {
        // If I just overwrite the XY_model::h and call this thing everytime I want to call the force
        if(s_eq_t < t) {
            // if we are in the quench phase, we reduce T
            XY_model::h = max(h_start - (h_start / quench::T_start) * (t - s_eq_t)/tau, h_end);
        }
    }
public:
    template<class State, class Deriv>
    void calc_force(State &x, Deriv &dxdt, double t) {
        linear_h(t);
        XY_silicon_anisotrop_subsystems::calc_force(x, dxdt, t);
    }

    void print_info() override {
        XY_silicon_anisotrop_subsystems::print_info();
        quench::print_info();
    }

    XY_silicon_anisotrop_subsystems_quench_h(map<Parameter, double> paras) : quench(paras),
                                                                           XY_silicon_anisotrop_subsystems(paras),
                                                                           XY_Silicon(paras),
                                                                           XY_model(paras),
                                                                           subsystems(paras),
                                                                           System(paras){
        h_start = XY_model::h;
        h_end = (quench::T_end - quench::T_start) * quench::tau;
        cout << "XY_silicon_anisotrop_subsystems_quench system constructed";
    }

    string get_name() const override {
        return "XY_silicon_anisotrop_subsystems_quench_h";
    }
};

struct XY_silicon_anisotrop_subsystems_quench_obc : public XY_silicon_anisotrop_subsystems_obc, public quench {
public:
    void print_info() override {
        XY_silicon_anisotrop_subsystems::print_info();
        quench::print_info();
    }

    XY_silicon_anisotrop_subsystems_quench_obc(map<Parameter, double> paras) : quench(paras),
                                                                           XY_silicon_anisotrop_subsystems_obc(paras),
                                                                           XY_silicon_anisotrop_subsystems(paras),
                                                                           XY_Silicon(paras),
                                                                           XY_model(paras),
                                                                           subsystems(paras),
                                                                           System_OBC(paras),
                                                                           System(paras){
        cout << "XY_silicon_anisotrop_subsystems_quench system constructed";
    }
    string get_name() const override {
        return "XY_silicon_anisotrop_subsystems_quench_obc";
    }
};
#endif //CUDAPROJECT_SYSTEMS_CUDA_CUH