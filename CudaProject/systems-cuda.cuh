//
// Created by andi on 30.04.23.
//

#ifndef CUDAPROJECT_SYSTEMS_CUDA_CUH
#define CUDAPROJECT_SYSTEMS_CUDA_CUH

#include "main.cuh"
// #include "parameters.cuh"
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
#include <curand_kernel.h>

using namespace std;


struct System {
public:
    const size_t n;

    // needed to "advance the random engine" and generate different random numbers for the different steps
    size_t step_nr;
    double T;
    const double eta;
    double D;
    // possibility of different sizes in the different directions
    const size_t dim_size_x;
    const size_t dim_size_y;
    thrust::device_vector<curandState> curand_states;

    checkpoint_timer timer {{}};           // checkpoint timer with no names, for now only total time

    template<class State>
    void init_state(map<Parameter, double>& paras, State &x) {
        cout << "standard init state called" << endl;
    }

    virtual void print_info() {
        cout << "Printing system info:" << endl;
        cout << "step_nr = " << step_nr << endl;
        cout << "T = " << T << endl;
        cout << "eta = " << eta << endl;
        cout << "D = " << D << endl;
        cout << "dim_size_x = " << dim_size_x << endl;
        cout << "dim_size_y = " << dim_size_y << endl;
        cout << "n = " << n << endl;
    }

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


    template<class State, class Deriv, class Stoch>
    void operator()(const State &x, Deriv &dxdt, Stoch &theta, double t) {}

    struct Functor {
        Functor(){}

        virtual void operator()() {}

    };

    template<class State, class Deriv, class FunctorType>
    void derivative_calculation(State &x, Deriv &dxdt, double t, FunctorType functor) {
        BOOST_AUTO(start, thrust::make_zip_iterator(thrust::make_tuple(
                x.begin(),
                x.begin() + n,
                dxdt.begin(),
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
        step_nr++;
    }

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
        step_nr++;
    }

    template<class Stoch>
    void calc_diff_thrust(Stoch &theta, double t) {
        // TODO actually t is never needed here with the current architecture, but i am too lazy to fix that
        // as it will probably improve nothing
        chrono::microseconds mus = chrono::duration_cast<chrono::microseconds >(
                chrono::system_clock::now().time_since_epoch()
        );
        // TODO are you still using this seed? Why is it not possible to get it to work with the default random
        // engine just advancing? Maybe just switch to curand?
        long seed = (mus.count() % 10000000) * 1000000000;
        thrust::counting_iterator<size_t> index_sequence_begin(seed);

        thrust::transform(index_sequence_begin,
                          index_sequence_begin + n,
                          theta.begin() + n,
                          rand(D));

    }
    template<class Stoch>
    void calc_diff(Stoch &theta, double t) {
        auto start = thrust::make_zip_iterator(thrust::make_tuple(theta.begin() + n, curand_states.begin()));
        // TODO does it work with start + n?
        thrust::for_each(start, start + n, curand(D));
    }

    template<class Stoch>
    void calc_diff_bbk(Stoch &theta, double dt) {
        auto start = thrust::make_zip_iterator(thrust::make_tuple(theta.begin(), curand_states.begin()));
        thrust::for_each(start, start + (2 * n), curand(1.0));      // first n are n_1, second n are n_2
        // okay i want the first n entries of theta to be sqrt(2T / eta) zeta_2, the second n entries will be sqrt(2Teta) zeta_1
        // the problem is that i need th stepsize for this, which i find a bit unelegent, but whatever
        // so this line is supposed to transform the n_2 that sit in the second part of theta to zeta_2
        thrust::transform(theta.begin(), theta.begin() + n, theta.begin() + n, theta.begin() + n, zeta_2<double>(dt, eta, T));

        // the second one is doable with a lambda, i just need to recall how to write those...
        double pref = sqrt(2.0 * T * eta * tau_2(dt, eta));
        thrust::transform(theta.begin(), theta.begin() + n, theta.begin(), zeta_1<double>(pref));
    }

    template <class value_type>
    __host__ __device__
    value_type static tau_1(value_type dt, value_type eta) {
        return 1.0 / eta * (1.0 - exp(- eta * dt));
    }

    template <class value_type>
    __host__ __device__
    value_type static tau_2(value_type dt, value_type eta) {
        return 1.0 / (2.0 * eta) * (1.0 - exp(-2.0 * eta * dt));
    }

    template <class T>
    struct zeta_1
    {
        T pref;
        zeta_1(T dt, T eta, T temp) {
            // TODO the thing is those tau_1 and tau_2 are static as soon as we switched to a constant stepsize..
            pref = sqrt(2 * temp * eta * tau_2(dt, eta));
        }
        zeta_1(T pref): pref(pref) {}
        __host__ __device__
        T operator()(const T& x) const {
            return pref * x;
        }
    };

    template <class T>
    struct zeta_2
    {
        T dt, Tau_1, Tau_2, pref;
        zeta_2(T dt, T eta, T temp): dt(dt){
            // TODO the thing is those tau_1 and tau_2 are static as soon as we switched to a constant stepsize..
            Tau_1 = tau_1(dt, eta);
            Tau_2 = tau_2(dt, eta);
            pref = sqrt(2.0 * temp * eta);
        }
        zeta_2(T dt, T eta, T temp, T Tau_1, T Tau_2): dt(dt), Tau_1(Tau_1), Tau_2(Tau_2){
            pref = sqrt(2.0 * temp / eta);
        }
        __host__ __device__
        T operator()(const T& x, const T& y) const {
            return pref * (((Tau_1 - Tau_2) / sqrt(Tau_2)) * x + sqrt(dt - pow(Tau_1, 2) / Tau_2) * y);
        }
    };

    template<class State>
    void map_state(State &x) {
        // Okay this maps the state to a certain interval, doesnt do anything for most systems, but for XY and dipol
        // system we map onto the intervals for the angles
    }

    template <class T>
    struct var
    {
        T mean;
        var(T mean): mean(mean) {}

        __host__ __device__
        T operator()(const T& x) const {
            return pow(mean - x, 2);
        }
    };

    // depricated constructors, use the one below
    System(size_t step_nr, const double eta, const double T, const size_t lat_dim) : step_nr(step_nr), dim_size_x(lat_dim), dim_size_y(lat_dim), n(lat_dim * lat_dim), eta(eta), T(T), D(sqrt(2 * T * eta)) {
        cout << "System constructor is called with eta = " << eta << "  T = " << T << endl;
    }
    System(map<string, double>& paras) : step_nr((size_t)paras["step_nr"]), dim_size_x((size_t)paras["lat_dim"]), dim_size_y((size_t)paras["lat_dim"]), eta(paras["eta"]), T(paras["T"]),
                                         n((size_t)paras["lat_dim"] * (size_t)paras["lat_dim"]) {
        D = (sqrt(2 * T * eta));
        cout << "System constructor from Map is called with eta = " << eta << "  T = " << T << endl;
    }

    System(map<Parameter, double>& paras) : step_nr((size_t)paras[Parameter::step_nr]), dim_size_x((size_t)paras[Parameter::dim_size_x]),
                                            dim_size_y((size_t)paras[Parameter::dim_size_y]), eta(paras[Parameter::eta]), T(paras[Parameter::T]),
                                            n((size_t)paras[Parameter::dim_size_x] * (size_t)paras[Parameter::dim_size_y]) {
        D = (sqrt(2 * T * eta));
        cout << "System constructor from Enumeration type Map is called with eta = " << eta << "  T = " << T << endl;
        curand_states = thrust::device_vector<curandState>(2 * n);  // TODO is this okay?
        // counting iterator counting from 0 to n. Every lattice site needs its own state i think
        // a state that corresponds to a different sequence. That the sequences are only one apart should not be relevant
        auto sInit = thrust::make_zip_iterator(thrust::make_tuple(thrust::counting_iterator<int>(0), curand_states.begin()));
        // now we call curand init with sequence numbers from the counting iterator
        thrust::for_each(sInit, sInit + 2*n, curand_setup(step_nr));      // now 2n, moare states should not be a problem and then we can use it also for the bbk random numbers?

    }

    size_t get_dim_size_x() const{
        return dim_size_x;
    }

    size_t get_dim_size_y() const{
        return dim_size_y;
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

    template <class T>
    struct kinetic_energy
    {
        __host__ __device__
        T operator()(const T& x, const T& y) const {
            return 0.5 * x * y;
        }
    };

    template<class State>
    double calc_kinetic_energy(State &x) {
        double E_kin = 0.5 * thrust::transform_reduce(x.begin() + n, x.end(), square<double>(), 0.0, thrust::plus<double>());
        return E_kin;
    }

    template<class State>
    double calc_f_mm(State &x) {
        double m2 = thrust::transform_reduce(x.begin(), x.begin() + n, square<double>(), 0.0, thrust::plus<double>()) / (double) n;
        double m = thrust::reduce(x.begin(), x.begin() + n, 0.0, thrust::plus<double>()) / (double) n;       // should work?
        double f_mm = (double)n * (m2 / (m * m) - 1);
        return f_mm;
    }



    template<class State, class Functor>
    double calc_f_me(State &x, Functor functor) {
        // ah okay damn, the energy depends on the system and the interaction potential
        // how do i solve this the easiest? do i hand over a functor that knows the potential?
        // ah damn it gets pretty complicated if I want to do this the right way... because of the neighbors
        // it will look like in the universal opterations function
        thrust::device_vector<double> e(n);
        thrust::device_vector<double> me(n);
        BOOST_AUTO(start, thrust::make_zip_iterator(thrust::make_tuple(
                x.begin(),
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
                ),
                e.begin()
        )));
        thrust::for_each(start, start + n, functor);        // potential energy
        // TODO check if this is correct.
        thrust::device_vector<double> e_p(n);
        thrust::transform(x.begin() + n, x.end(), x.begin() + n, e_p.begin(), kinetic_energy<double>());    // kinetic energy
        thrust::transform(e.begin(), e.end(), e_p.begin(), e.begin(), thrust::plus<double>());               // add both into e
        thrust::transform(x.begin(), x.begin() + n, e.begin(), me.begin(), thrust::multiplies<double>());
        double me_avg = thrust::reduce(me.begin(), me.end(), 0.0, thrust::plus<double>()) / (double) n;
        double m = thrust::reduce(x.begin(), x.begin() + n, 0.0, thrust::plus<double>()) / (double) n;       // should work?
        double e_avg = thrust::reduce(e.begin(), e.end(), 0.0, thrust::plus<double>()) / (double) n;
        double f_me = (double)n * (me_avg / (m * e_avg) - 1);
        return f_me;
    }

    template<class State>
    double calc_f_me(State &x) {
        return 0;
    }

    template<class State>
    double calc_m(State &x) {
        double m = thrust::reduce(x.begin(), x.begin() + n, 0.0, thrust::plus<double>()) / (double) n;       // should work?
        return m;
    }

    double test() {
        return 1.0;
    }

    size_t get_step_nr() {
        return step_nr;
    }

    double get_eta() {
        return eta;
    }
};


struct NNN_System: public System {
    using System::System;       // inherit constructor

    struct up_right : public up, public right{
        up_right(size_t dimension_size_x, size_t dimension_size_y):
                right(dimension_size_x), up(dimension_size_x, dimension_size_y){}
        __host__ __device__ size_t operator()(size_t i) const override{
            // up right is just the upper neighbor of the right neighbor
            return right::operator()(up::operator()(i));
        }
    };

    struct up_left : public up, public left{
        up_left(size_t dimension_size_x, size_t dimension_size_y):
                left(dimension_size_x), up(dimension_size_x, dimension_size_y){}
        __host__ __device__ size_t operator()(size_t i) const override{
            // up right is just the upper neighbor of the right neighbor
            return left::operator()(up::operator()(i));
        }
    };
    struct down_right : public down, public right{
        down_right(size_t dimension_size_x, size_t dimension_size_y):
                right(dimension_size_x), down(dimension_size_x, dimension_size_y){}
        __host__ __device__ size_t operator()(size_t i) const override{
            // down right is just the downper neighbor of the right neighbor
            return right::operator()(down::operator()(i));
        }
    };
    struct down_left : public down, public left{
        down_left(size_t dimension_size_x, size_t dimension_size_y):
                left(dimension_size_x), down(dimension_size_x, dimension_size_y){}
        __host__ __device__ size_t operator()(size_t i) const override{
            // down left is just the downper neighbor of the left neighbor
            return left::operator()(down::operator()(i));
        }
    };

    template<class State, class Deriv, class FunctorType>
    void derivative_calculation(State &x, Deriv &dxdt, double t, FunctorType functor) {
        BOOST_AUTO(start, thrust::make_zip_iterator(thrust::make_tuple(
                x.begin(),
                x.begin() + n,
                dxdt.begin(),
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
                ),
                // okay tuple can only have 10 entries so we make another tuple here?
                thrust::make_zip_iterator(thrust::make_tuple(thrust::make_permutation_iterator(
                                           x.begin(),
                                           thrust::make_transform_iterator(
                                                   thrust::counting_iterator<size_t>(0),
                                                   down_right(dim_size_x, dim_size_y)
                                           )
                                           ),
                                           thrust::make_permutation_iterator(
                                                   x.begin(),
                                                   thrust::make_transform_iterator(
                                                           thrust::counting_iterator<size_t>(0),
                                                           down_left(dim_size_x, dim_size_y)
                                                   )
                                           ),
                                           thrust::make_permutation_iterator(
                                                   x.begin(),
                                                   thrust::make_transform_iterator(
                                                           thrust::counting_iterator<size_t>(0),
                                                           up_right(dim_size_x, dim_size_y)
                                                   )
                                           ),
                                           thrust::make_permutation_iterator(
                                                   x.begin(),
                                                   thrust::make_transform_iterator(
                                                           thrust::counting_iterator<size_t>(0),
                                                           up_left(dim_size_x, dim_size_y)
                                                   )
                                           )))
        )));
        thrust::for_each(start, start + n, functor);
        step_nr++;
    }

};


struct constant_bath : public System {
    using left = typename System::left;
    using right = typename System::right;
    using up = typename System::up;
    using down = typename System::down;
public:
    const double alpha;
    const double beta;
    const double J;
    // systemsize should probably be a template argument?
    // In contrast to the brownian system we need one more parameter for the timescale of the cooling down
    // the current temperature

    void print_info() override {
        System::print_info();
    }

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
        timer.set_startpoint(functor_point);

        bath_functor functor = bath_functor(System::eta, alpha, beta, J);

        this->derivative_calculation(x, dxdt, t, functor);
        timer.set_endpoint(functor_point);
    }

    template<class Stoch>
    void calc_diff(Stoch &theta, double t) {
        timer.set_startpoint(rng);
        System::calc_diff(theta, t);
        timer.set_endpoint(rng);
    }
public:
    constant_bath(const double T, const double eta, const double alpha, const double beta, const double J, const size_t lat_dim, const int init_step=0)
            : System(init_step, eta, T, lat_dim), alpha(alpha), beta(beta), J(J) {
        cout << "constant Bath System is constructed" << endl;
    }

    double get_cur_T() const{
        return System::T;
    }
    ~constant_bath() {
        cout << "Bath System is destroyed" << endl;
    }
};



struct coulomb_interaction : public System {
    using left = typename System::left;
    using right = typename System::right;
    using up = typename System::up;
    using down = typename System::down;
public:
    const double alpha;
    const double beta;
    const double J;


    void print_info() override {
        System::print_info();
        cout << "alpha = " << alpha << endl;
        cout << "beta = " << beta << endl;
        cout << "J = " << J << endl;
    }
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

    struct energy_functor {

    };

public:
    template<class State, class Deriv>
    void calc_drift(const State &x, Deriv &dxdt, double t) {
        coulomb_functor functor = coulomb_functor(System::eta, alpha, beta, J);
        this->derivative_calculation(x, dxdt, t, functor);
    }



    coulomb_interaction(const double T, const double eta, const double alpha, const double beta, const double J, const size_t lat_dim, const int init_step=0)
            : System(init_step, eta, T, lat_dim), alpha(alpha),
            beta(beta), J(J) {
    }
    coulomb_interaction(map<string, double>& paras)
            : System(paras), alpha(paras["alpha"]), beta(paras["beta"]), J(paras["J"]) {
    }
    coulomb_interaction(map<Parameter, double>& paras) : System(paras), alpha(paras[Parameter::alpha]), beta(paras[Parameter::beta]), J(paras[Parameter::J]) {}

};


class coulomb_constant : public coulomb_interaction {
public:
    coulomb_constant(const double T, const double eta, const double alpha, const double beta, const double J, const size_t lat_dim, const int init_step=0)
    : coulomb_interaction(T, eta, alpha, beta, J, lat_dim, init_step) {

    }
    coulomb_constant(map<string, double>& paras)
            : coulomb_interaction(paras) {
    }

    coulomb_constant(map<Parameter, double>& paras)
            : coulomb_interaction(paras) {
    }

    template<class Stoch>
    void calc_diff(Stoch &theta, double t) {
        System::calc_diff(theta, t);
    }
};



struct gpu_bath : public coulomb_interaction {
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
        System::T = max(T_start - t/tau, T_end);
        return System::T;
    }

    double linear_T(double t) {
        if(s_eq_t < t && t < end_quench_t) {
            // if we are in the quench phase, we reduce T
            System::T = T_start - (t - s_eq_t)/tau;
        }
        return System::T;
    }

    template<class Stoch>
    void calc_diff(Stoch &theta, double t) {
        // I think this should work? We change the D of the System and then just calc the random numbers
        System::D = sqrt(2 * linear_T(t) * System::eta);
        System::calc_diff(theta, t);
    }
public:
    gpu_bath(const double T, const double T_end, const double eta, const double alpha, const double beta, const double J, const double tau, const size_t lat_dim, size_t init_step = 0, double eq_t = 30)
            : coulomb_interaction(T, eta, alpha, beta, J, lat_dim, init_step), T_start(T),
              tau(tau), T_end(T_end), s_eq_t(eq_t), e_eq_t(eq_t) {
        t_quench = (get_quench_time());
        end_quench_t = t_quench + s_eq_t;       // end of quench is the quench time + equilibrate at beginning
    }

    // TODO Constructor with map
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


struct quench : virtual public System {
    const double    T_start;       // Start-Temperture of the Quenching For example: 10
    const double    T_end;         // End-Temperature of the Quencheing. For example: 1
    const double    s_eq_t;        // start equilibration time: Time the system gets to equilibrate at T_start
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
        if(s_eq_t < t && t < end_quench_t) {
            // if we are in the quench phase, we reduce T
            System::T = T_start - (t - s_eq_t)/tau;
        }
        return System::T;
    }

    template<class Stoch>
    void calc_diff(Stoch &theta, double t) {
        // I think this should work? We change the D of the System and then just calc the random numbers
        System::D = sqrt(2 * linear_T(t) * System::eta);
        System::calc_diff(theta, t);
    }

public:
   quench(const double T, const double T_end, const double eta, const double tau, const size_t lat_dim, size_t init_step = 0, double eq_t = 30)
            : System(init_step, eta, T, lat_dim), T_start(T), tau(tau), T_end(T_end), s_eq_t(eq_t), e_eq_t(eq_t) {
        t_quench = (get_quench_time());
        end_quench_t = t_quench + s_eq_t;       // end of quench is the quench time + equilibrate at beginning
    }
    quench(map<string, double>& paras): System(paras),
            tau(paras["tau"]), T_end(paras["end_T"]), s_eq_t(paras["t_eq"]), e_eq_t(paras["t_eq"]), T_start(paras["starting_T"]) {
        t_quench = (get_quench_time());
        end_quench_t = t_quench + s_eq_t;
    }

    quench(map<Parameter, double>& paras): System(paras),
                                        tau(paras[Parameter::tau]), T_end(paras[Parameter::end_temp]), s_eq_t(paras[Parameter::equil_time]),
                                        e_eq_t(paras[Parameter::equil_time]), T_start(paras[Parameter::starting_temp]) {
        t_quench = (get_quench_time());
        end_quench_t = t_quench + s_eq_t;
    }

    double get_quench_time() {
        // returns the time it takes to do the quench
        // in this system, we use a linear quench
        cout << "running get_quench_time:" << endl << "T_start = " << T_start << endl << "T_end = " << T_end << endl << "tau = " << tau << endl << endl;
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


struct anisotropic_coulomb_interaction : virtual public System {
    using left = typename System::left;
    using right = typename System::right;
    using up = typename System::up;
    using down = typename System::down;
public:
    const double alpha;
    const double beta;
    const double Jx, Jy;

    void print_info() override {
        System::print_info();
        cout << "alpha = " << alpha << endl;
        cout << "beta = " << beta << endl;
        cout << "Jx = " << Jx << endl;
        cout << "Jy = " << Jy << endl;
    }

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

    struct potential_energy_functor {
        const double alpha, beta, Jx, Jy;

        potential_energy_functor(const double alpha,
                            const double beta, const double Jx, const double Jy) :
                            alpha(alpha), beta(beta), Jx(Jx), Jy(Jy){ }

        template<class Tup>
        __host__ __device__ void operator()(Tup tup) {
            double q = thrust::get<0>( tup );
            double q_left =     thrust::get<1>(tup);
            double q_right =    thrust::get<2>(tup);
            double q_up =       thrust::get<3>(tup);
            double q_down =     thrust::get<4>(tup);

            // now here the logic for the energy
            double E_dw = 0.5 * (alpha * pow(q, 4) - beta * pow(q, 2));
            double E_interaction = Jx * (1.0 / sqrt(1.0 + (q-q_left) * (q-q_left))
                            + 1.0 / sqrt(1.0 + (q-q_right) * (q-q_right)))
                            + Jy * (1.0 / sqrt(1.0 + (q-q_up) * (q-q_up))
                                    + 1.0 / sqrt(1.0 + (q-q_down) * (q-q_down)));
            thrust::get<5>(tup) = E_dw + E_interaction;
        }
    };


public:
    template<class State, class Deriv>
    void calc_drift(const State &x, Deriv &dxdt, double t) {
        ani_coulomb_functor functor = ani_coulomb_functor(System::eta, alpha, beta, Jx, Jy);
        this->derivative_calculation(x, dxdt, t, functor);
    }

    anisotropic_coulomb_interaction(const double T, const double eta, const double alpha, const double beta, const double Jx, const double Jy, const size_t lat_dim, const int init_step=0)
            : System(init_step, eta, T, lat_dim), alpha(alpha),
              beta(beta), Jx(Jx), Jy(Jy) {

    }
    anisotropic_coulomb_interaction(map<string, double>& paras)
            : System(paras), alpha(paras["alpha"]), beta(paras["beta"]), Jx(paras["J"]), Jy(paras["Jy"]) {
    }
    anisotropic_coulomb_interaction(map<Parameter, double>& paras)
            : System(paras), alpha(paras[Parameter::alpha]), beta(paras[Parameter::beta]), Jx(paras[Parameter::J]), Jy(paras[Parameter::Jy]) {
    }

    template<class State>
    double calc_f_me(State &x) {
        auto functor = potential_energy_functor(alpha, beta, Jx, Jy);
        return System::calc_f_me(x, functor);
    }

};



class anisotropic_coulomb_constant : public anisotropic_coulomb_interaction {
public:
    anisotropic_coulomb_constant(const double T, const double eta, const double alpha, const double beta, const double Jx, const double Jy, const size_t lat_dim, const int init_step=0)
            : System(init_step, eta, T, lat_dim), anisotropic_coulomb_interaction(T, eta, alpha, beta, Jx, Jy, lat_dim, init_step) {
        cout << "creating anisotropic_coulomb_constant system" << endl;
    }

    anisotropic_coulomb_constant(map<string, double>& paras)
            : System(paras), anisotropic_coulomb_interaction(paras) {
        cout << "creating anisotropic_coulomb_constant system" << endl;
    }

    anisotropic_coulomb_constant(map<Parameter, double>& paras)
            : System(paras), anisotropic_coulomb_interaction(paras) {
        cout << "creating anisotropic_coulomb_constant system" << endl;
    }


    template<class Stoch>
    void calc_diff(Stoch &theta, double t) {
        // cout << "calc_diff is called" << endl;
        System::calc_diff(theta, t);
    }
};


class anisotropic_coulomb_quench: public anisotropic_coulomb_interaction, public quench {
public:
    anisotropic_coulomb_quench(const double T, const double T_end, const double eta, const double alpha,
                               const double beta, const double Jx,const double Jy, const double tau, const size_t lat_dim,
                               size_t init_step = 0, double eq_t = 30)
            : anisotropic_coulomb_interaction(T, eta, alpha, beta, Jx, Jy, lat_dim, init_step),
                    quench(T, T_end, eta, tau, lat_dim, init_step, eq_t),
                            System(init_step, eta, T, lat_dim){

    }
    anisotropic_coulomb_quench(map<string,double>& paras)
            : anisotropic_coulomb_interaction(paras), quench(paras), System(paras){

    }

    anisotropic_coulomb_quench(map<Parameter,double>& paras)
            : anisotropic_coulomb_interaction(paras), quench(paras), System(paras){

    }

    void print_info() override {
        quench::print_info();
        anisotropic_coulomb_interaction::print_info();
    }

};

template <class value_type>
__host__ __device__
inline value_type positive_modulo(value_type i, value_type n) {
    return fmod((fmod(i, n) + n), n);
}

class XY_model : virtual public System {
    using left = typename System::left;
    using right = typename System::right;
    using up = typename System::up;
    using down = typename System::down;
    const double p_XY = 2;     // for the potential, will be 2 for bistable XY and 2.5 for silicon
    const double m = 1;     // prefactor for the interaction, important if we dont have the 2 pi periodicy
protected:

    double J;
    double h;

    struct xy_functor {
        const double J, h, eta, p_XY, m;

        xy_functor(const double eta, const double J, const double h, const double p_XY, const double m) :
        J(J), h(h), eta(eta), p_XY(p_XY), m(m) { }

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

            double interaction = J * (
                    +   sin(m * (q - q_up))   + sin(m * (q - q_down))
                    +   sin(m * (q - q_left)) + sin(m * (q - q_right))
            );
            thrust::get<3>( tup ) = (-eta) * p                                                                                  // Friction
                                    + p_XY * h * sin(p_XY * q) // bistable potential
                                    - interaction;       // Interaction
        }
    };

    struct xy_force {
        const double J, h, eta, p_XY, m;

        xy_force(const double eta, const double J, const double h, const double p_XY, const double m) :
                J(J), h(h), eta(eta), p_XY(p_XY), m(m) {
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
            double interaction = J * (
                    +   sin(m * (q - q_up))   + sin(m * (q - q_down))
                    +   sin(m * (q - q_left)) + sin(m * (q - q_right))
            );
            thrust::get<1>( tup ) = p_XY * h * sin(p_XY * q) // bistable potential
                                    - interaction;       // Interaction
        }
    };

    struct potential_energy_functor {
        const double J, h;

        potential_energy_functor(const double J, const double h) : J(J), h(h) { }

        template<class Tup>
        __host__ __device__ void operator()(Tup tup) {
            double q = thrust::get<0>( tup );
            double q_left =     thrust::get<1>(tup);
            double q_right =    thrust::get<2>(tup);
            double q_up =       thrust::get<3>(tup);
            double q_down =     thrust::get<4>(tup);

            // now here the logic for the energy
            double E_pot = h * cos(2 * q);
            double E_interaction = J * (cos(q - q_left) + cos(q - q_right) + cos(q - q_up) + cos(q - q_down));
            thrust::get<5>(tup) = E_pot + E_interaction;
        }
    };

    struct map_functor {
        template<class value_type>
        __host__ __device__ value_type operator()(value_type q) {
            return positive_modulo(q, 2.0 * M_PI);
        }
    };

    struct init_functor {
        size_t n;
        double equil_pos, range_min, range_max;

        init_functor(size_t n, double equil_pos, double range_min, double range_max): n(n), equil_pos(equil_pos), range_min(range_min), range_max(range_max) {}
        template<class State>
        void operator()(map<Parameter, double>& paras, State &x) {
            cout << "XY init_state is called with equil_pos = " << equil_pos << endl;
            if(paras[random_init] == 0.0) {
                // equilibrium initialization -> we are in XY model with p=2, meaning we have
                // our equilibria at pi/2 and 3pi/2, we initialize everything in the pi/2 minimum
                thrust::fill(x.begin(), x.begin()+n, equil_pos);

            } else {
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
        init_functor(n, M_PI/2, 0, 2 * M_PI)(paras, x);
    }

    void print_info() override{
        System::print_info();
        cout << "h = " << h << endl;
        cout << "J = " << J << endl;
    }

    template<class State, class Deriv>
    void calc_drift(State &x, Deriv &dxdt, double t) {
        xy_functor functor = xy_functor(System::eta, J, h, p_XY, m);
        System::derivative_calculation(x, dxdt, t, functor);
    }

    template<class State, class Deriv>
    void calc_force(State &x, Deriv &dxdt, double t) {
        xy_force functor = xy_force(System::eta, J, h, p_XY, m);
        System::force_calculation(x, dxdt, t, functor);
    }


    template<class State>
    double calc_f_me(State &x) {
        auto functor = potential_energy_functor(J, h);
        return System::calc_f_me(x, functor);
    }

    XY_model(map<Parameter,double>& paras)
    : System(paras), J(paras[Parameter::J]), h(paras[alpha]), m(paras[Parameter::m]), p_XY(paras[Parameter::p]) {
        cout << "Xy model constructor is called "<< endl;
        print_info();
    }
};


class XY_quench: public XY_model, public quench {
public:
    XY_quench(map<Parameter, double>& paras): XY_model(paras), quench(paras), System(paras) {
        cout << "XY Quench System is constructed" << endl;
    }
    void print_info() override {
        quench::print_info();
        XY_model::print_info();
    }

};


class XY_Silicon: virtual public XY_model {
    // we only have to set the map functor new? And obviously the right p, m stuff
    const double p_XY = 2.57;           // does this work? Can i just redeclare them here?
    const double m = 2.0;

    struct map_functor {        // TODO Question if this already works without redefining map_state?
        template<class value_type>
        __host__ __device__ value_type operator()(value_type q) {
            return positive_modulo((q + M_PI/2), M_PI) - M_PI/2;
        }
    };

    template<class State, class Deriv>
    void calc_drift(State &x, Deriv &dxdt, double t) {
        cout << "Using p_XY = " << p_XY << endl;
        XY_model::xy_functor functor = XY_model::xy_functor(System::eta, J, h, p_XY, m);
        System::derivative_calculation(x, dxdt, t, functor);
    }

    template<class State>
    void init_state(map<Parameter, double>& paras, State &x) {
        double equil_pos = (7.0 / 18.0) * M_PI;
        double range_min = - M_PI / 2;
        double range_max = M_PI / 2;
        init_functor(n, equil_pos, range_min, range_max)(paras, x);
    }

    XY_Silicon(map<Parameter,double>& paras): System(paras), XY_model(paras) {}
};

#define DIAG_DIST 11.18033988749895
#define DIAG_PREF 1.4
#define P_COS   2.5714285714285716

class dipol_interaction: public NNN_System {
    double p_mom;   // q * l dipole moment
    double h;   // strength of anisotropy potential

public:

    struct dipol_functor {
        double p_mom;
        double h;
        double eta;
        dipol_functor(double eta, double p_mom, double h): eta(eta), p_mom(p_mom), h(h) {}

        template<class Tup>
        __host__ __device__ void operator()(Tup tup) {
            double q = thrust::get<0>( tup );
            double p = thrust::get<1>( tup );
            double q_left = thrust::get<4>(tup);
            double q_right = thrust::get<5>(tup);
            double q_up = thrust::get<6>(tup);
            double q_down = thrust::get<7>(tup);

            // thrust::get<8>(tup) should get us the tuple with the 4 NNN, so we get the value with get<i>
            double q_down_right = thrust::get<0>(thrust::get<8>(tup));
            double q_down_left = thrust::get<1>(thrust::get<8>(tup));
            double q_up_right = thrust::get<2>(thrust::get<8>(tup));
            double q_up_left = thrust::get<0>(thrust::get<8>(tup));
            double interaction = (DIAG_PREF * cos(q) * sin(q_down_right) + sin(q) * cos(q_down_right) +
                            DIAG_PREF * cos(q) * sin(q_down_left) + sin(q) * cos(q_down_left) +
                            DIAG_PREF * cos(q) * sin(q_up_right) + sin(q) * cos(q_up_right) +
                            DIAG_PREF * cos(q) * sin(q_up_left) + sin(q) * cos(q_up_left)
                            )  / DIAG_DIST;                         // NNN Interaction

            interaction += 1.0 / 8.0 * (
                    2 * cos(q) * sin(q_right) + sin(q) * cos(q_right) +         // right neighbor
                    2 * cos(q) * sin(q_left) + sin(q) * cos(q_left)           // left neighbor
            );
            interaction += - (
                     - cos(q) * sin(q_up) + sin(q) * cos(q_up) +         // up neighbor
                     - cos(q) * sin(q_down) + sin(q) * cos(q_down));           // down neighbor

            thrust::get<2>( tup ) = p;
            thrust::get<3>( tup ) = (-eta) * p                                                                                  // Friction
                                       + P_COS * h * sin(P_COS * q) // bistable potential
                                       - pow(p_mom, 2) * interaction;       // Interaction
        }
    };

    struct potential_energy_functor {
        const double J, h;

        potential_energy_functor(const double J, const double h) : J(J), h(h) { }

        template<class Tup>
        __host__ __device__ void operator()(Tup tup) {
            double q = thrust::get<0>( tup );
            double q_left =     thrust::get<1>(tup);
            double q_right =    thrust::get<2>(tup);
            double q_up =       thrust::get<3>(tup);
            double q_down =     thrust::get<4>(tup);

            // now here the logic for the energy
            double E_pot = h * cos(P_COS * q);
            // TODO not implemented
            double E_interaction = 0;
            thrust::get<5>(tup) = E_pot + E_interaction;
        }
    };

    struct map_functor {
        template<class value_type>
        __host__ __device__ value_type operator()(value_type q) {
            return positive_modulo((q + M_PI/2), M_PI) - M_PI/2;
        }
    };

    void print_info() override {
        System::print_info();
        cout << "h = " << h << endl;
        cout << "p_mom = " << p_mom << endl;
    }

    template<class State, class Deriv>
    void calc_drift(State &x, Deriv &dxdt, double t) {
        dipol_functor functor = dipol_functor(System::eta, p_mom, h);
        NNN_System::derivative_calculation(x, dxdt, t, functor);
    }
    template<class State>
    void map_state(State &x) {
        // maybe we don't have to redeclare map state again all the time?
        // we map every q onto [-pi/2, pi/2]
        thrust::transform(x.begin(), x.begin() + n, x.begin(), map_functor());
    }


    template<class State>
    double calc_f_me(State &x) {
        auto functor = potential_energy_functor(J, h);
        return System::calc_f_me(x, functor);
    }

    dipol_interaction(map<Parameter,double>& paras)
    : NNN_System(paras), p_mom(paras[Parameter::J]), h(paras[alpha]) {
        cout << "dipol model constructor is called "<< endl;
        print_info();
    }

};


struct quadratic_trapped_lattice : public System {
    using rand = typename System::rand;
    using left = typename System::left;
    using right = typename System::right;
    using up = typename System::up;
    using down = typename System::down;
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

    struct linear_force {
        const double J, alpha;

        linear_force(const double J, const double alpha) :
                J(J), alpha(alpha) {
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
            double interaction = J * ((q - q_left) + (q - q_right) + (q - q_up) + (q - q_down));

            thrust::get<1>( tup ) =  -alpha * (2 * q)  // linear force
                                    + interaction;       // Interaction
        }
    };

    template<class State, class Deriv>
    void calc_drift(const State &x, Deriv &dxdt, double t) {
        timer.set_startpoint(functor_point);

        harmonic_trap_functor functor = harmonic_trap_functor(System::eta, alpha, J);

        System::derivative_calculation(x, dxdt, t, functor);
        timer.set_endpoint(functor_point);
    }

    template<class State, class Deriv>
    void calc_force(State &x, Deriv &dxdt, double t) {
        linear_force functor = linear_force(J, alpha);
        System::force_calculation(x, dxdt, t, functor);
    }

    template<class Stoch>
    void calc_diff(Stoch &theta, double t) {
        timer.set_startpoint(rng);
        System::calc_diff(theta, t);
        timer.set_endpoint(rng);
    }

public:
    quadratic_trapped_lattice(const double T, const double eta, const double alpha, const double J, const size_t lat_dim, const int init_step=0)
            : System(init_step, eta, T, lat_dim), alpha(alpha), J(J) {
        cout << "Lattice quadratic Trap System is constructed" << endl;
    }

    quadratic_trapped_lattice(map<Parameter, double> paras): System(paras), alpha(paras[Parameter::alpha]), J(paras[Parameter::J]) {
        cout << "Lattice quadratic Trap System from parameter map is constructed" << endl;
    }

    double get_cur_T() const{
        return System::T;
    }
};


struct chain : virtual public System {
    using left = typename System::left;
    using right = typename System::right;

    template<class State, class Deriv, class FunctorType>
    void derivative_calculation(State &x, Deriv &dxdt, double t, FunctorType functor) {
        BOOST_AUTO(start, thrust::make_zip_iterator(thrust::make_tuple(
                x.begin(),
                x.begin() + n,
                dxdt.begin(),
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
                thrust::counting_iterator<size_t>(0)
        )));
        thrust::for_each(start, start + n, functor);
        step_nr++;
    }

    chain(map<Parameter, double> paras): System(paras) {}

};

struct XY_pairs : public chain, public XY_model {
    double J;
    double h;
    struct XY_pair_functor {
        const double J, h, eta;

        XY_pair_functor(const double J, const double h, const double eta): J(J), h(h), eta(eta) {

        }

        template<class Tup>
        __host__ __device__ void operator()(Tup tup) {
            double q = thrust::get<0>( tup );
            double p = thrust::get<1>( tup );
            thrust::get<2>( tup ) = p;
            size_t i = thrust::get<6>(tup);
            // depending on i either the left or the right neighbor is my partner
            // even -> right .... uneven -> left
            double q_partner = 0;
            if (i % 2 == 0) {   // 5 is right
                q_partner = thrust::get<5>(tup);
            } else {            // 4 is left
                q_partner = thrust::get<4>(tup);
            }

            thrust::get<3>( tup ) = (-eta) * p                         // Friction
                                    + h * sin(2 * q)                // on site potential
                                    - J * sin(q - q_partner);       // interaction
        }
    };

public:
    XY_pairs(map<Parameter,double>& paras) : chain(paras), XY_model(paras), System(paras),
    J(paras[Parameter::J]), h(paras[Parameter::alpha]) {
        cout << "XY_pairs system is constructed" << endl;
    }

    template<class State, class Deriv>
    void calc_drift(State &x, Deriv &dxdt, double t) {
        XY_pair_functor functor = XY_pair_functor(J, h, chain::eta);
        chain::derivative_calculation(x, dxdt, t, functor);
    }
};

struct quadratic_chain : public chain {
public:
    const double J;

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


public:
    template<class State, class Deriv>
    void calc_drift(const State &x, Deriv &dxdt, double t) {
        functor Functor = functor(J, eta);
        chain::derivative_calculation(x, dxdt, t, Functor);
    }
    quadratic_chain(map<Parameter, double> paras): chain(paras), J(paras[Parameter::J]) {}
};



struct gpu_oscillator_chain : chain{
public:
    const double alpha;
    // no tau, constant temp
    // no beta, just harmonic potential
    // no interaction
    // systemsize should probably be a template argument?
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

public:
    gpu_oscillator_chain(map<Parameter, double> paras): chain(paras), alpha(paras[Parameter::alpha]) {
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

#endif //CUDAPROJECT_SYSTEMS_CUDA_CUH