//
// Created by andi on 21.04.23.
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
#include <fstream>
#include <filesystem>



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
};


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
    gpu_bath(const double T, const double eta, const double alpha, const double beta, const double J, const double tau)
            : T(T), step_nr(0), n(lat_dim * lat_dim), T_t(T), eta(eta), alpha(alpha), beta(beta), J(J), tau(tau) {
    }

    double get_cur_T() const{
        return T_t;
    }

    size_t get_lattice_dim() const{
        return lat_dim;
    }
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
        count = 0;
    }

    template<class State, class System>
    void operator()(const System sys, const State &x , double t ) {
        // so i think if statements are not to costly on c++
        if (count % write_every == 0) {
            // now we write the state to the file
            // i call a function here that would not be available for all Systems, not to clean and a bit fishy
            // but i guess it works for now
            size_t lat_dim = sys.get_lattice_dim();
            double T = sys.get_cur_T();

            // writing, but only q for the moment
            file << "t, " << t << ",";
            for(int i = 0; i < lat_dim * lat_dim; i++) {
                    file << x[i] << ",";

            }
            // for the last we write T? I mean thats how i did it, right?
            // Temperatur saven
            file << T;

            // Zeilenumbruch
            file << "\n";
        }
        count++;
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
class euler_mayurama_stepper {
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
        sys(x, dxdt, theta, t);
        // this should set the correct values for dxdt and theta so that they can be applied in apply_em
        // can we print them here?

//        cout << "x[0] = " << x[0] << ", " << "x[1] = " << x[1] << endl;
//        cout << "dxdt[0] = " << dxdt[0] << ", " << "dxdt[1] = " << dxdt[1] << endl;
//        cout << "theta[0] = " << theta[0] << ", " << "theta[1] = " << theta[1] << endl;

        // for the update we need to define a type of a function, but i don't really understand the syntax
        // okay so the plan is to have a new name for the templated struct apply_em that is a member of operations
        // so that we can later more clearly instantiate the struct apply_em with the correct template
        typedef typename operations::template apply_em<time_type> apply_em;
        // this is the operation that i want to apply on every lattice site
        algebra::for_each(x, x, dxdt, theta, apply_em(dt));
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
                                                                                {}

/*    euler_mayurama_stepper(size_t N) : N(N) {
        Observer = new observer();
    }*/

    // now for the memory allocation. I am not sure if this ever changes for me but the implementation should not harm
    // on the other hand, my class doesn't have any saved state_types that could potentially be resized
    // so we skip this for now
private:
    // the system size, but what is N if i am in 2 dimensions? probably n * n. But how can i initialize the inital
    // values for all sites? not just by "stat_type(N)"
    // i still don't fully understand what the state_type is, is it just (x, p) for one lattice site? Or is it the
    // the state of the whole system? In the former, N would correspond with the dimensionality of my DGL system
    // in the latter i would not really know how to initialize the initial states independent of the dimensionality
    // of my problem
    const size_t N;
    state_type dxdt, theta;
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

int main() {
    // We try out the code for the brownian motion i would say
    // But we cannot use our old class system I think because there the whole system is already on a lattice
    // But we can quickly write another system i guess
    const int steps = 1000000;
    const double dt = 0.005;
    const double eta2 = 5;
    const double T = 70;
    const double J = 50;
    const double alpha = 5;
    const double beta = 10;
    const double tau = 100;
    const double eta = 5;
    const int nr_save_values = 32;
    size_t write_every = steps / nr_save_values;
    const size_t lattice_dim = 100;
    // system size
    const size_t n = lattice_dim * lattice_dim;
    // DGLs per lattice site
    const size_t N = 2;
    
    cout << "Starting Simulation for a " << lattice_dim << " by " << lattice_dim << " lattice for " << steps << "steps." << endl;
    // last time i didnt have to specify the dimensionality i think (in terms of (x, p) )
    const double D = T / eta2;
    double theo_msd = 2 * D * dt * steps;
    double theo_mu = 0;
    double mu = 0;
    double msd = 0;

/*    observer base_obs = observer();

    euler_mayurama_stepper<state_type, container_algebra, default_operations> stepper(2, &base_obs);
    brownian_particel system(eta2, T);

    // I guess we directly have to check whether the distribution parameters are still the same.
    const size_t runs = 10;



    for(size_t i = 0; i < runs; i++) {
        // init the inital values for every run
        state_type x(2, 0.0);
        for( size_t n=0 ; n<steps ; ++n ) {
            stepper.do_step(system, x, dt, n*dt);
            // cout << n*dt << " ";
            // cout << x[0] << " " << x[1] << endl;
        }
        // add to msd and mu
        mu += x[0];
        msd += x[0] * x[0];
    }
    // averaging
    mu /= runs;
    msd /= runs;

    cout << "mu = " << mu << endl;
    cout << "msd = " << msd << "   theo value msd = " << theo_msd << endl;*/

    // file stuff
    string storage_root = "../../../Generated content/GPU Test/";
    string dir_name = create_tree_name(eta, T, dt, n, alpha, beta, J, tau,
                                       storage_root);
    string name = dir_name + "/0";
    ofstream file;
    ofstream parafile;
    file.open(name + ".csv");
    parafile.open(name + ".txt");

    // init observer

    bath_observer Obs(file, write_every);

    typedef thrust::device_vector<double> gpu_state_type;
    euler_mayurama_stepper<gpu_state_type, thrust_algebra, thrust_operations > gpu_stepper(N * n);

    // initialize the system...
    // gpu_bath(const double T, const double eta, const double alpha, const double beta, const double J, const double tau)
    gpu_bath<lattice_dim> gpu_system(T, eta, alpha, beta, J, tau);
    gpu_state_type x(N * n, 0.0);
    /*


    // We initialize a system of size 50000... what does that even mean?
    gpu_brownian_system gpu_system(eta, T, n);

    // so this state type only has N=2 values, which are set to 0. Maybe we need to initialize N * n values?
    gpu_state_type x(N * n, 0.0);
      */

    auto start = chrono::high_resolution_clock::now();
    double t = 0;
    for( size_t i=0 ; i<steps ; ++i ) {
        gpu_stepper.do_step(gpu_system, x, dt, t);
        // ... why don't we just write here? wouldn't that be faster?
        // i guess we do that
        Obs(gpu_system, x, t);
        t += dt;
        // cout << n*dt << " ";
        // cout << x[0] << " " << x[1] << endl;
    }

    write_parameters(parafile, eta, T, dt, n, alpha, beta, J, tau);

    file.close();
    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::milliseconds>(end - start);
    cout << "execution took " << duration.count() << "ms, meaning " <<
         duration.count() * 1000/(steps) << "ms per 1000 steps." << endl;
    cout << "for a " << lattice_dim << " by " << lattice_dim << " lattice." << endl;
    // print this shit
    // TODO we could use this reduction stuff to compute the moments
    mu = 0;
    msd = 0;
    for(int i = 0; i <= n; i++) {
        mu += x[i];
        msd += x[i] * x[i];
    }
    mu /= n;
    msd /= n;

    cout << "mu = " << mu << endl;
    cout << "msd = " << msd << "   theo value msd = " << theo_msd << endl;


    // add to msd and mu
    return 0;
}
