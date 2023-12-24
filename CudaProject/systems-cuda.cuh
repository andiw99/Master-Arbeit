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
#include <cuda_runtime.h>
#include <cufft.h>


// #include "cufft_utils.h"

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
    bool curand_random = false;

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
    void calc_diff_curand(Stoch &theta, double t) {
        auto start = thrust::make_zip_iterator(thrust::make_tuple(theta.begin() + n, curand_states.begin()));
        // TODO does it work with start + n?
        thrust::for_each(start, start + n, curand(D));
    }

    template<class Stoch>
    void calc_diff(Stoch &theta, double t) {
        if(curand_random) {
            calc_diff_curand(theta, t);
        } else {
            calc_diff_thrust(theta, t);
        }
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
                                            n((size_t)paras[Parameter::dim_size_x] * (size_t)paras[Parameter::dim_size_y]),
                                            curand_random((bool) paras[Parameter::curand_random]){
        D = (sqrt(2 * T * eta));
        cout << "System constructor from Enumeration type Map is called with eta = " << eta << "  T = " << T << endl;
        cout << "Using curand_random: " << curand_random << endl;
        if(curand_random) {
            curand_states = thrust::device_vector<curandState>(2 * n);  // TODO is this okay?
            auto sInit = thrust::make_zip_iterator(thrust::make_tuple(thrust::counting_iterator<int>(0), curand_states.begin()));
            // counting iterator counting from 0 to n. Every lattice site needs its own state i think
            // a state that corresponds to a different sequence. That the sequences are only one apart should not be relevant
            // now we call curand init with sequence numbers from the counting iterator
            thrust::for_each(sInit, sInit + 2*n, curand_setup(step_nr));      // now 2n, moare states should not be a problem and then we can use it also for the bbk random numbers?
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
    struct sum_square_complex
    {
        __host__ __device__
        double operator()(const T& x, double y) const {
            return y + x.x * x.x + x.y * x.y;
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

    double test() {
        return 1.0;
    }

    size_t get_step_nr() {
        return step_nr;
    }

    double get_eta() {
        return eta;
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

struct System_anitemp: virtual public System {
public:
    thrust::device_vector<double> D;
    struct curand
    {
        double mu, sigma;

        __host__ __device__
        curand(double mu = 0.0, double sigma = 1.0) : mu(mu), sigma(sigma) {};

        template<class Tuple>
        __device__
        float operator()(Tuple tup) const
        {
            curandState local_state = thrust::get<1>(tup);
            // in this class the D is location dempendent
            thrust::get<0>(tup) = thrust::get<2>(tup) * curand_normal(&local_state);
            thrust::get<1>(tup) = local_state;
        }
    };
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
                          rand(System::D));

    }

    template<class Stoch>
    void calc_diff_curand(Stoch &theta, Stoch& D_vec, double t) {
        auto start = thrust::make_zip_iterator(thrust::make_tuple(theta.begin() + n, curand_states.begin(), D_vec.begin()));
        // TODO does it work with start + n?
        thrust::for_each(start, start + n, curand());
    }

    template<class Stoch>
    void calc_diff(Stoch &theta, Stoch &D_vec, double t) {
        if(curand_random) {
            calc_diff_curand(theta, D_vec, t);
        } else {
            calc_diff_thrust(theta, D_vec, t);
        }
    }

    System_anitemp(map<Parameter, double>& paras) : System(paras) {
        D = thrust::device_vector<double>(n);
        thrust::fill(D.begin(), D.end(), System::D);
    }
    void print_info() override {
        System::print_info();
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

struct quench_ani : virtual public System_anitemp, virtual public quench {
    void print_info() override {
        System_anitemp::print_info();
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

    template<class Stoch, class FunctorType>
    void calc_diff(Stoch &theta, double t, FunctorType functor) {
        // I think this should work? We change the D of the System and then just calc the random numbers
        // okay here we need the logic on how to calculate the D for the different stuffs, probably with a functor
        // that is then implemented in the child classes?
        auto tuple = thrust::make_zip_iterator(thrust::make_tuple(D.begin(), thrust::counting_iterator<size_t>(0)));
        thrust::for_each(tuple, tuple + n, functor);
        System_anitemp::calc_diff(theta, D, t);
    }

public:
    quench_ani(map<Parameter, double>& paras): System_anitemp(paras), quench(paras), System(paras) {
        t_quench = (get_quench_time());
        end_quench_t = t_quench + s_eq_t;
    }

};

struct quench_left_right : virtual public quench_ani {
    double T_diff;
public:
    struct  quench_functor  {
        double s_eq_t, T_end, T_start, T_diff, tau, t, eta;
        size_t Lx;
        quench_functor(double s_eq_t, double T_end, double T_start, double T_diff, double tau,
                       size_t Lx, double t, double eta):
        s_eq_t(s_eq_t), T_end(T_end), T_start(T_start), T_diff(T_diff), tau(tau), t(t), eta(eta), Lx(Lx) {}

        template<class Tup>
        __host__ __device__
        void operator()(Tup tup) {
                // we have a tuple of D and a counting iterator so we know the place in the lattice
                // how can we model the temperature?
                if(s_eq_t <= t) {
                    // we need to know how far we are from the left end of our subsystem
                    // it will just be i % Lx?
                    int lx = thrust::get<1>(tup) % Lx;
                    double T = max((T_start + lx * (T_diff / (double)Lx)) - (t - s_eq_t) / tau, T_end);
                    thrust::get<0>(tup) = sqrt(2 * eta * T);
                }
        }
    };
    template<class Stoch>
    void calc_diff(Stoch&theta, double t) {
        quench_functor functor(quench::s_eq_t, quench::T_end, quench::T_start, T_diff, quench::tau, get_Lx(),
                               t, eta);
        quench_ani::calc_diff(theta, t, functor);
    }

    quench_left_right(map<Parameter, double>& paras) : quench_ani(paras), quench(paras),
    System(paras), System_anitemp(paras) {
        T_diff = 0.5 * (T_start - T_end);
        // we will initialize the system with the anisotrope temperature
        auto tuple = thrust::make_zip_iterator(thrust::make_tuple(D.begin(), thrust::counting_iterator<size_t>(0)));
        // the time that will enter the functor will be the s_eq_t so that the if statement works
        quench_functor functor(quench::s_eq_t, quench::T_end, quench::T_start, T_diff, quench::tau, get_Lx(),
                               quench::s_eq_t, eta);
        thrust::for_each(tuple, tuple + n, functor);
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

    struct potential_energy_functor {
        const double Jx, Jy, h;

        potential_energy_functor(const double J, const double h) : Jx(J), Jy(J), h(h) { }
        potential_energy_functor(const double Jx, const double Jy, const double h) : Jx(Jx), Jy(Jy), h(h) { }

        template<class Tup>
        __host__ __device__ void operator()(Tup tup) {
            double q = thrust::get<0>( tup );
            double q_left =     thrust::get<1>(tup);
            double q_right =    thrust::get<2>(tup);
            double q_up =       thrust::get<3>(tup);
            double q_down =     thrust::get<4>(tup);

            // now here the logic for the energy
            double E_pot = h * cos(2 * q);
            double E_interaction = Jx * (cos(q - q_left) + cos(q - q_right)) + Jy * (cos(q - q_up) + cos(q - q_down));
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
        size_t n, dim_size_x;
        double equil_pos, range_min, range_max;

        init_functor(size_t n, double equil_pos, double range_min, double range_max, size_t dim_size_x): n(n), equil_pos(equil_pos), range_min(range_min), range_max(range_max), dim_size_x(dim_size_x) {}
        template<class State>
        void operator()(map<Parameter, double>& paras, State &x) {
            cout << "XY init_state is called with equil_pos = " << equil_pos << endl;
            if(paras[random_init] == 0.0) {
                // equilibrium initialization -> we are in XY model with p=2, meaning we have
                // our equilibria at pi/2 and 3pi/2, we initialize everything in the pi/2 minimum
                thrust::fill(x.begin(), x.begin()+n, equil_pos);
                if(paras[Parameter::J] < 0) {
                    cout << "initializing with chess trafo" << endl;
                    chess_trafo_rectangular(x, dim_size_x);
                }

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
        init_functor(n, M_PI/2, 0, 2 * M_PI, dim_size_x)(paras, x);
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
public:
    const double p_XY = 2.57;           // does this work? Can i just redeclare them here?
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
    void calc_drift(State &x, Deriv &dxdt, double t) {
        XY_model::xy_functor functor = XY_model::xy_functor(System::eta, J, h, p_XY, m);
        System::derivative_calculation(x, dxdt, t, functor);
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

    XY_Silicon(map<Parameter,double>& paras): System(paras), XY_model(paras) {}

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
    void calc_drift(State &x, Deriv &dxdt, double t) {
        XY_model::xy_functor functor = XY_model::xy_functor(System::eta, Jx, Jy, h, p_XY, m);
        System::derivative_calculation(x, dxdt, t, functor);
    }

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

    struct subsystem_index : thrust::unary_function<size_t, size_t> {
        const size_t dim_size_x, Lx, subsystem_ind;
        subsystem_index(size_t dimension_size_x, size_t Lx, size_t subsystem_ind):
                thrust::unary_function<size_t, size_t>(), dim_size_x(dimension_size_x), Lx(Lx), subsystem_ind(subsystem_ind){
        }
        virtual __host__ __device__ size_t operator()(size_t i) const {
            // first we need to know the row of the index
            int row = i / Lx;
            int col = i % Lx;
            int index_in_whole_system = row * dim_size_x + subsystem_ind * Lx + col;
            return index_in_whole_system;
        }
    };
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

    struct running_chess_trafo : thrust::unary_function<size_t, size_t> {
        size_t dim_size_x, Lx;
        running_chess_trafo(size_t dimension_size_x, size_t Lx):
                thrust::unary_function<size_t, size_t>(),
                        dim_size_x(dimension_size_x), Lx(Lx){
        }
        template<class Tup>
        __host__ __device__ void operator()(Tup tup) const {
            // first we need to know in which system we are
            if (dim_size_x % 2 == 1) {
                // if it is uneven we can just trafo every even
                if(thrust::get<1>(tup) % 2 == 0) {
                    thrust::get<0>(tup) = (-1) * thrust::get<0>(tup);
                }
            } else {
                // else we trafo in even rows even indices and in uneven rows
                // uneven indices
                int ind = thrust::get<1>(tup);
                int row_even = (ind / Lx) % 2;
                if(ind % 2 - row_even == 0) {
                    thrust::get<0>(tup) = (-1) * thrust::get<0>(tup) ;
                }
            }
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

    struct real_to_complex {
        template<class Tup>
        __host__ __device__ void operator()(Tup tup) const {
            thrust::get<0>(tup).x = thrust::get<1>(tup);
            thrust::get<0>(tup).y = 0.0;
            printf("%.4f   ", thrust::get<0>(tup).x);
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
        step_nr++;
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
        step_nr++;
    }

    subsystems_obc(map<Parameter, double> paras): System_OBC(paras), System(paras), subsystems(paras) {
        cout << "initializing subsystems pbc with sizes: " << Lx << " x " << Ly << " = " << n0 << endl;
        if (Ly != dim_size_y) {
            cout << "dim_size_y has to be equal to L_y" << endl;
            exit(0);
        }
    }
};

struct XY_subsystems : public subsystems_pbc, public XY_model {
public:
    XY_subsystems(map<Parameter, double> paras): subsystems_pbc(paras), subsystems(paras), XY_model(paras), System(paras) {}

    template<class State, class Deriv>
    void calc_force(State &x, Deriv &dxdt, double t) {
        xy_force functor = xy_force(System::eta, J, h, XY_model::p_XY, XY_model::m);
        subsystems_pbc::force_calculation(x, dxdt, t, functor);
    }
};

struct XY_silicon_subsystems : public subsystems_pbc, public XY_Silicon {
public:
    XY_silicon_subsystems(map<Parameter, double> paras) : XY_model(paras),
                                                          subsystems(paras),
                                                          subsystems_pbc(paras),
                                                          XY_Silicon(paras),
                                                          System(paras) {}
    template<class State, class Deriv>
    void calc_force(State &x, Deriv &dxdt, double t) {
        xy_force functor = xy_force(System::eta, J, h, XY_Silicon::p_XY, XY_Silicon::m);
        subsystems_pbc::force_calculation(x, dxdt, t, functor);
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
        vector<double> m_vec{};
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
        BOOST_AUTO(tuple, thrust::make_zip_iterator(thrust::make_tuple(cell, thrust::counting_iterator<size_t>(0))));
        /*thrust::for_each(tuple, tuple + (dim_size_x * Ly), running_chess_trafo(dim_size_x, Lx));*/
        auto cell_trafo = thrust::make_transform_iterator(tuple, running_chess_trafo_iterator(dim_size_x, Lx));
/*        for (int j = 0; j < dim_size_x * Ly; j++) {
            cout << "pos in sub: " << j << " pos in whole: " << subsystem_inds[j] << " value:" << cell_trafo[j] << endl;
        }*/
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
        double m_L2 = std::transform_reduce(m_vec.begin(), m_vec.end(),
                                            0.0, // initial value for the reduction (sum)
                                            std::plus<double>(), // transformation (square)
                                            [](double m) -> double { return m * m; });
        m_L2 /= m_vec.size();
        double m_L4 = std::transform_reduce(m_vec.begin(), m_vec.end(),
                                            0.0, // initial value for the reduction (sum)
                                            std::plus<>(), // transformation (square)
                                            [](double m) { return (pow(m, 4)); });
        m_L4 /= m_vec.size();

        double cum = m_L4 / (m_L2 * m_L2);

        return cum;
    }

    struct cos_real_to_complex {
        template<class Tup>
        __host__ __device__ void operator()(Tup tup) const {
            thrust::get<0>(tup).x = cos(thrust::get<1>(tup));
            thrust::get<0>(tup).y = 0.0;
        }
    };

    struct sin_real_to_complex {
        template<class Tup>
        __host__ __device__ void operator()(Tup tup) const {
            thrust::get<0>(tup).x = sin(thrust::get<1>(tup));
            thrust::get<0>(tup).y = 0.0;
        }
    };

    template<class State>
    void calc_xi(State& x, double& xix, double& xiy) {
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

        /*cout << endl << endl;
        thrust::host_vector<cufftComplex> host_output(output_vector);
        cout << "cos transform fft output: " << endl;
        for(int i = 0; i < n; i++) {
            cout << host_output[i].x << "  " << host_output[i].y << endl;
        }*/

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

/*        host_output = thrust::host_vector<cufftComplex>(output_vector);
        thrust::host_vector<cufftComplex> host_input(input_vector);
        cout << endl << endl;
        cout << "sin transform fft output: " << endl;
        for(int i = 0; i < n; i++) {
            if (i % batch_size == 0) {
                cout << endl;
            }
            cout << *(cell_trafo + i) << "   " << host_input[i].x << "   " << host_input[i].y << "   " << host_output[i].x << "  " << host_output[i].y << endl;
        }*/

         // sum up
        for(int i = 0; i < nr_batches; i++) {
            thrust::transform(output_vector.begin() + i * batch_size, output_vector.begin() + (i + 1) * batch_size,
                              sum.begin(), sum.begin(), sum_square_complex<cufftComplex>());
        }

        thrust::host_vector<double> host_sum(sum);

        for(int i = 0; i < batch_size; i++) {
            host_sum[i] /= nr_batches;
        }

/*        cout << endl << endl;
        for(int i = 0; i < batch_size; i++) {
            cout << host_sum[i] << endl;
        }*/

        double* ft_squared_k = new double[Lx];
        double* ft_squared_l = new double[Ly];
        for(int i = 0; i < Lx; i++) {
            ft_squared_k[i] = 0.0;
        }
        // cout << endl;
        for(int i = 0; i < Ly; i++) {
            ft_squared_l[i] = 0.0;
        }
        for(int i = 0; i < Lx; i++) {
            for(int j = 0; j < Ly; j++) {
                int k_ind = j * Lx + i;
                ft_squared_k[i] += host_sum[k_ind];
                ft_squared_l[j] += host_sum[k_ind];
            }
        }
        // cout << endl;
        for(int i = 0; i < Lx; i++) {
            ft_squared_k[i] /= Ly;
            // cout << ft_squared_k[i] << "  ";
        }
        //cout << endl;
        for(int i = 0; i < Ly; i++) {
            ft_squared_l[i] /= Lx;
            // cout << ft_squared_l[i] << "  ";
        }
        // cout << endl;
        auto kx = get_frequencies_fftw_order(Lx);
        auto ky = get_frequencies_fftw_order(Ly);

        Eigen::VectorXd paras_x = fit_lorentz_peak(kx, ft_squared_k);
        Eigen::VectorXd paras_y = fit_lorentz_peak(ky, ft_squared_l);


        xix = paras_x(1);
        xiy = paras_y(1);

        // cufftDestroy(plan);
        // cudaDeviceReset();
        //cudaFree(input_pointer);
        //cudaFree(output_pointer);
        delete[] ft_squared_k;
        delete[] ft_squared_l;
        cufftDestroy(plan);
        // cudaDeviceReset();

        cout << "xix " << xix << "  xiy " << xiy << endl;
        //exit(0);
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

        thrust::host_vector<double> host_sum(sum);

        for(int i = 0; i < batch_size; i++) {
            host_sum[i] /= nr_batches;
        }

        for(int i = 0; i < Lx; i++) {
            ft_k[i] = 0.0;
        }
        // cout << endl;
        for(int i = 0; i < Ly; i++) {
            ft_l[i] = 0.0;
        }
        for(int i = 0; i < Lx; i++) {
            for(int j = 0; j < Ly; j++) {
                int k_ind = j * Lx + i;
                ft_k[i] += host_sum[k_ind];
                ft_l[j] += host_sum[k_ind];
            }
        }
        // cout << endl;
        for(int i = 0; i < Lx; i++) {
            ft_k[i] /= Ly;
        }
        //cout << endl;
        for(int i = 0; i < Ly; i++) {
            ft_l[i] /= Lx;
        }

        cufftDestroy(plan);

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

};

struct XY_silicon_subsystems_quench : public XY_silicon_subsystems, public quench {
public:
    void print_info() override {
        XY_silicon_subsystems::print_info();
        quench::print_info();
    }

    XY_silicon_subsystems_quench(map<Parameter, double> paras) : quench(paras),
                                                                 XY_silicon_subsystems(paras),
                                                                 subsystems(paras),
                                                                 XY_model(paras),
                                                                 System(paras){
        cout << "XY_silicon_subsystems_quench system constructed";
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

};

struct XY_silicon_anisotrop_subsystems_quench_LR_obc : public XY_silicon_anisotrop_subsystems_obc, public quench_left_right {
public:
    void print_info() override {
        XY_silicon_anisotrop_subsystems::print_info();
        quench::print_info();
    }

    XY_silicon_anisotrop_subsystems_quench_LR_obc(map<Parameter, double> paras) :   quench(paras),
                                                                                    quench_ani(paras),
                                                                                    quench_left_right(paras),
                                                                                    XY_silicon_anisotrop_subsystems_obc(paras),
                                                                                    XY_silicon_anisotrop_subsystems(paras),
                                                                                    XY_Silicon(paras),
                                                                                    XY_model(paras),
                                                                                    subsystems(paras),
                                                                                    System_OBC(paras),
                                                                                    System_anitemp(paras),
                                                                                    System(paras){
        cout << "XY_silicon_anisotrop_subsystems_quench system constructed";
    }

};

#endif //CUDAPROJECT_SYSTEMS_CUDA_CUH