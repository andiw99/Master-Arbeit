//
// Created by andi on 30.04.23.
//

#ifndef CUDAPROJECT_MAIN_CUH
#define CUDAPROJECT_MAIN_CUH
//

#include "../LearningProject/Header/Helpfunctions and Classes.h"
#include <functional>
#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <boost/typeof/typeof.hpp>
#include <thrust/iterator/counting_iterator.h>
#include <thrust/random.h>
#include <thrust/reduce.h>
#include <thrust/sequence.h>
#include <thrust/execution_policy.h>
#include <map>

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

    int get_elapsed_time() const {
        auto endtime = chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
                endtime - starttime);
        int duration_count = (int)duration.count();
        return duration_count;
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

struct constant_bath_timer : public timer {
    long rng_generation_time = 0;
    long functor_time = 0;
    chrono::time_point<chrono::high_resolution_clock> rng_start;
    chrono::time_point<chrono::high_resolution_clock> rng_end;
    chrono::time_point<chrono::high_resolution_clock> functor_start;
    chrono::time_point<chrono::high_resolution_clock> functor_end;
public:
    constant_bath_timer() : timer() {
        cout << "timer is constructed" << endl;
        rng_start = chrono::high_resolution_clock ::now();
        rng_end = chrono::high_resolution_clock ::now();
        functor_start = chrono::high_resolution_clock ::now();
        functor_end = chrono::high_resolution_clock ::now();
    }

    void set_rng_start() {
        rng_start = chrono::high_resolution_clock ::now();
    }

    void set_rng_end() {
        rng_end = chrono::high_resolution_clock ::now();
        rng_generation_time += std::chrono::duration_cast<std::chrono::microseconds>(
                rng_end - rng_start).count();
    }

    void set_functor_start() {
        functor_start = chrono::high_resolution_clock ::now();
    }

    void set_functor_end() {
        functor_end = chrono::high_resolution_clock ::now();
        functor_time += std::chrono::duration_cast<std::chrono::microseconds>(
                functor_end - functor_start).count();
    }

    ~constant_bath_timer() {
        cout << "RNG took " << (double)rng_generation_time * 0.001 << " ms" << endl;
        cout << "Functor executions took " << (double)functor_time * 0.001 << " ms" << endl;
    }
};

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

template <class State>
void fill_init_values(State &state, float x0, float p0, int run = 0, double mu=0, double sigma=1) {

    int n = state.size()/2;

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

template <class state_type>
class state_initializer {
public:
    map<string, double> paras;
    state_initializer(map<string, double>& paras): paras(paras) {

    }
    virtual void init_state(state_type& x) {
        cout << "Dummy initializer is called" << endl;
    }

};

template <class state_type>
class random_initializer : public state_initializer<state_type>  {
    double x0;
    double p0;
    size_t lattice_dim;
public:
    random_initializer(map<string, double>& paras): state_initializer<state_type> (paras) {
        const double beta = paras["beta"];
        x0 = paras["x0"] * sqrt(beta / 2.0);
        p0 = paras["p0"] * sqrt(beta / 2.0);
        lattice_dim = (size_t) paras["lat_dim"];
    }
    void init_state(state_type& x) {
        chrono::milliseconds ms = chrono::duration_cast<chrono::milliseconds >(
                chrono::system_clock::now().time_since_epoch()
        );
        fill_init_values<state_type>(x, (float) x0, (float) p0, ms.count() % 10000);
    }
};

template <class state_type>
class memory_initializer : public state_initializer<state_type> {
    fs::path root;
    int count = 0;
    size_t lattice_dim;

public:
    memory_initializer(map<string, double>& paras, fs::path& root): state_initializer<state_type> (paras), root(root) {
        lattice_dim = (size_t) paras["lat_dim"];
    }
    void init_state(state_type& x) {
        // since the initialization sometimes depend on the parameters
        // like here T, i have to create the initializer everytime i call repeat for the first time
        double T = this->paras["T"];
        size_t n = lattice_dim * lattice_dim;
        cout << "checking for initial state in folder..." << endl;
        // listing the temp folders that are already inside
        vector<fs::path> temp_paths = list_dir_paths(root);
        // we check every folder name for the value and use the one that is closest to our actual temp
        string closest_T = findClosestDir(temp_paths, T);
        cout << "clostest folder to " << T << " already existing is " << closest_T << endl;
        if(closest_T != "None") {
            // now we list every csv file and take the one that is closest to i
            vector<fs::path> csv_files = list_csv_files(closest_T);
            string closest_i = findClosestStem(csv_files, count);
            cout << "clostest index to " << count << " already existing is " << closest_i << endl;
            if(closest_i != "None") {
                fs::path pre_name = closest_i;
                cout << "Trying to read " << pre_name << endl;
                ifstream pre_file = safe_read(pre_name, true);
                // and we are ready to read in the last file?
                double prev_T;
                double prev_t;
                vector<double> pre_lattice = readDoubleValuesAt(pre_file, -1, prev_T, prev_t);
                for(int i = 0; i < n; i++) {
                    x[i] = pre_lattice[i];
                }
            }
        }

        // increment count
        count++;
    }
};

template <class state_type>
class equilibrium_initializer : public state_initializer<state_type>  {
    size_t n;
    double J;
    double beta;
public:
    equilibrium_initializer(map<string, double> & paras) : state_initializer<state_type> (paras) {
        n = (size_t)paras["n"];
        J = paras["J"];
        beta = paras["beta"];
    }

    void init_state(state_type& x) {
        cout << get_name() << " is called" << endl;
        thrust::fill(x.begin(), x.begin() + n, sqrt(beta / 2.0));
        if (J < 0) {
            chess_trafo(x, (size_t)sqrt(n));
        }
        thrust::fill(x.begin() + n, x.begin() + 2*n, 0);
    }

    string get_name() {
        return "equilibrium initializer";
    }
};

template <class state_type>
state_initializer<state_type>* create_state_initializer(int random, map<string, double>& paras, fs::path& root) {
    // again some fishy function to create different derived classes based on the value of an parameter
    // don't know better, dont care...
    if(random == -1) {
        // memory initializer case
        return new memory_initializer<state_type> (paras, root);
    } else if (random == 1) {
        // case for the random initializer
        return new random_initializer<state_type> (paras);
    } else if (random == 0) {
        return new equilibrium_initializer<state_type> (paras);
    }
}




#endif //CUDAPROJECT_MAIN_CUH
