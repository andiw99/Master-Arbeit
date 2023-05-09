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

using namespace std;


template <size_t lat_dim>
struct gpu_bath {
public:
    const double T;
    const double alpha;
    const double eta;
    const double tau;
    const double beta;
    const double J;
    // systemsize should probably be a template argument?
    const size_t n;
    // needed to "advance the random engine" and generate different random numbers for the different steps
    size_t step_nr;
    // In contrast to the brownian system we need one more parameter for the timescale of the cooling down
    // the current temperature
    double T_t;
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
            /*


            // We now have 4 additional iterators in the tuple 4 -> 8
            double q = thrust::get<0>( tup );
            double p = thrust::get<1>( tup );
            // this line has to be the dgl for x
            // dqdt =
            thrust::get<2>( tup ) = p;
            // this the one for p
            // dpdt=

            // this is the one depending on the whole q stuff. We can use a normal calculation here right?
            thrust::get<3>( tup ) = (-eta) * p - dVdq(q, q_left, q_right, q_up, q_down);
            */
            double q = thrust::get<0>( tup );
            double p = thrust::get<1>( tup );
            // this line has to be the dgl for x
            // dqdt =
            thrust::get<2>( tup ) = p;
            // neighbors
            double q_left = thrust::get<4>(tup);
            double q_right = thrust::get<5>(tup);
            double q_up = thrust::get<6>(tup);
            double q_down = thrust::get<7>(tup);

            // this the one for p
            // dpdt=
            // I think the problem is that we call another function and this is not legal
            // i guess I will just use an inline implementation for my first implementation
            thrust::get<3>( tup ) = (-eta) * p                                                                                  // Friction
                                    - alpha * (2 * q * q * q - beta * q)                                                        // double well potential
                                    - J * (sin(q - q_left) + sin(q - q_right) + sin(q - q_up) + sin(q - q_down));       // Interaction
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

    double linear_T(double t) {
        // parametrisierung für die Temperatur
        // linearer Abfall
        T_t = max(T - t/tau, 1.0);
        return T_t;
    }

    template<class State, class Deriv, class Stoch>
    void operator()(const State &x, Deriv &dxdt, Stoch &theta, double t) {
        thrust::counting_iterator<size_t> index_sequence_begin(step_nr * n);
        if(step_nr == 0) {
            // TODO will this already be initialized to zero without this statement?
            thrust::fill(theta.begin(), theta.begin() + n, 0);
        }
        // One thing is that now T has to change, which we will realize like we have for the cpu case
        // In the cpu case i had a starting T and one that was changing, but why was that again?
        // i think so that we were able to get T(t), because otherwise your only option is to decrease T in every step
        // but then you don't know which time it is. But is that important?
        // We also have the step_nr but actually I think I will try to just decrement T in every step rn i don't see
        // the disadvantage
        // The disadvantage is that we don't have dt available here, so i will work with the cur_T attribute
        double D = sqrt(2 * linear_T(t) * eta);
        // ok so this should be done again, the linear_T computation is performed on the cpu i think which is good
        // since the Temperature is the same for all lattice sites and so it is only one computation
        // generate the random values
        thrust::transform(index_sequence_begin,
                          index_sequence_begin + n,
                          theta.begin() + n,
                          rand(D));
        /*
        cout << "theta[n-1] = " << theta[n-1] << endl;
        cout << "theta[n] = " << theta[n] << endl;
        cout << "theta[n+1] = " << theta[n+1] << endl;
        */

        // Problem now is that i have to generate iterators that iterate over all neighbors
        // We have an example of how to do this but i actually do not understand it
        // the most reasonable thing to do would be to just do it analogous and hope that it works
        // i only need the q values of the neighbors, so i get 4 additional iterators, one for each neighbor
        BOOST_AUTO(start, thrust::make_zip_iterator(thrust::make_tuple(
                // x begin has all q values
                x.begin(),
                // x begin + n has all p values
                x.begin() + n,
                // dxdt begin has all dqdt values
                dxdt.begin(),
                // dxdt begin + n has all dpdt values
                dxdt.begin() + n,
                // left neighbor q values
                // TODO i don't know what exactly is happening here but
                // i understand what effect it has or supposed to have
                thrust::make_permutation_iterator(
                        x.begin(),
                        thrust::make_transform_iterator(
                                thrust::counting_iterator<size_t>(0),
                                left()
                        )
                ),
                // right neighbor q values
                thrust::make_permutation_iterator(
                        x.begin(),
                        thrust::make_transform_iterator(
                                thrust::counting_iterator<size_t>(0),
                                right()
                        )
                ),
                // TODO this looks horibly bloated, but safest, fastest, easiest
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
        // bath_functor(const double eta, const double alpha,
        //                     const double beta, const double J)
        // apply bath_functor()() for each tuple between start and start + n
        thrust::for_each(start, start + n, bath_functor(eta, alpha, beta, J));

/*        cout << "x[1] = " << x[1] << endl;
        // the problem is here, actually dxdt[0] = x[1] should be
        // somehow dxdt is not set but i dont really get why
        cout << "dxdt[0] = " << dxdt[0] << endl;*/

        step_nr++;
    }
public:
    gpu_bath(const double T, const double eta, const double alpha, const double beta, const double J, const double tau, size_t init_step = 0)
            : T(T), step_nr(init_step), n(lat_dim * lat_dim), T_t(T), eta(eta), alpha(alpha), beta(beta), J(J), tau(tau) {
    }

    double get_cur_T() const{
        return T_t;
    }

    size_t get_lattice_dim() const{
        return lat_dim;
    }
};


template <size_t lat_dim>
struct constant_bath {
public:
    const double alpha;
    const double eta;
    const double beta;
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

            thrust::get<3>( tup ) = (-eta) * p                                                                                  // Friction
                                    - alpha * (2 * q * q * q - beta * q)                                                        // double well potential
                                    - J * (sin(q - q_left) + sin(q - q_right) + sin(q - q_up) + sin(q - q_down));       // Interaction
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
    struct left : thrust::unary_function<size_t, size_t> {
        __host__ __device__ size_t operator()(size_t i) const {
            return (i % lat_dim == 0) ? i + lat_dim - 1 : i - 1;
        }
    };

    struct right : thrust::unary_function<size_t, size_t> {
        __host__ __device__ size_t operator()(size_t i) const {
            return (i % lat_dim == lat_dim - 1) ? i - (lat_dim - 1) : i + 1;
        }
    };

    struct up : thrust::unary_function<size_t, size_t> {
        __host__ __device__ size_t operator()(size_t i) const {
            return (i < lat_dim) ? i + lat_dim * (lat_dim - 1) : i - lat_dim;
        }
    };

    struct down : thrust::unary_function<size_t, size_t> {
        __host__ __device__ size_t operator()(size_t i) const {
            return (i >= lat_dim * (lat_dim -1 )) ? i - lat_dim * (lat_dim - 1) : i + lat_dim;
        }
    };


    template<class State, class Deriv, class Stoch>
    void operator()(const State &x, Deriv &dxdt, Stoch &theta, double t) {
        thrust::counting_iterator<size_t> index_sequence_begin(step_nr * n);
        if(step_nr == 0) {
            // TODO will this already be initialized to zero without this statement?
            thrust::fill(theta.begin(), theta.begin() + n, 0);
        }

        thrust::transform(index_sequence_begin,
                          index_sequence_begin + n,
                          theta.begin() + n,
                          rand(D));
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
        thrust::for_each(start, start + n, bath_functor(eta, alpha, beta, J));
        step_nr++;
        // cout << "x[1] = " << x[1] << endl;
        // the problem is here, actually dxdt[0] = x[1] should be
        // somehow dxdt is not set but i dont really get why
        // cout << "dxdt[0] = " << dxdt[0] << endl;
    }
public:
    constant_bath(const double T, const double eta, const double alpha, const double beta, const double J)
            : T(T), step_nr(0), n(lat_dim * lat_dim), eta(eta), alpha(alpha), beta(beta), J(J), D(sqrt(2 * T * eta)) {
    }

    double get_cur_T() const{
        return T;
    }

    size_t get_lattice_dim() const{
        return lat_dim;
    }
};


template <size_t lat_dim>
struct gpu_oscillator_chain {
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
    // no neighbors


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



typedef vector<double> state_type;

struct brownian_particel {
    const double eta, T;
    static inline auto normal_dice = bind(normal_distribution<double>(0, 1), default_random_engine(13));

    brownian_particel(const double eta, const double T)
            : eta(eta), T(T) { }

    /*
     * call of the system needs to take in the theta state type
     */
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
