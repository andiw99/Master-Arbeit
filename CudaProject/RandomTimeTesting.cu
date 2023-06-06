//
// Created by weitze73 on 06.06.23.
//
#include <iostream>
#include "main.cuh"
#include "systems.cuh"
#include "parameters.cuh"


/**************************************************/
/* RANDOM NUMBERS GENERATION STRUCTS AND FUNCTION */
/**************************************************/
template<typename T>
struct rand_01 {
    __host__ T operator()(T& VecElem) const {


        return (T)rand() / RAND_MAX; }
};

template<typename T>
struct rand_01_for_each {
    int n = 0;
    __host__ __device__ void operator()(T& VecElem) {
        thrust::default_random_engine rng;
        thrust::uniform_real_distribution<float> dist(0, 1);
        n++;
        rng.discard(n);

        VecElem = dist(rng) / RAND_MAX;

    }
};

struct curand_slow{
    float a, b;
    curandState s;
    int seed;

    __host__ __device__
    curand_slow(float _a = 0.f, float _b = 1.f, int seed=0) : a(_a), b(_b), seed(seed) {

    };

    __device__
    float operator()(const unsigned int n)
    {
        curand_init(seed, n, 0, &s);
        return curand_normal(&s);
    }
};

struct curand_fast{
    float a, b;
    curandState s;
    int seed;

    __host__ __device__
    curand_fast(float _a = 0.f, float _b = 1.f, int seed=0) : a(_a), b(_b), seed(seed) {

    };

    template<class Tuple>
    __device__ void operator()(Tuple tup) const
    {
        curandState local_state;
        local_state = thrust::get<1>(tup);
        thrust::get<0>(tup) = curand_normal(&local_state);
    }
};
struct curand_fast_setup
{
    using init_tuple = thrust::tuple<int, curandState &>;
    __device__
    void operator()(init_tuple t){
        curandState s;
        int id = thrust::get<0>(t);
        curand_init(0, id, 0, &s);
        thrust::get<1>(t) = s;
    }
};

template<typename T>
__host__ T rand_01_fcn() { return ((T)rand() / RAND_MAX); }

struct prg
{
    float a, b;

    __host__ __device__
    prg(float _a = 0.f, float _b = 1.f) : a(_a), b(_b) {};

    __host__ __device__
    float operator()(const unsigned int n) const
    {
        thrust::default_random_engine rng;
        thrust::uniform_real_distribution<float> dist(a, b);
        rng.discard(n);

        return dist(rng);
    }
};

/********/
/* MAIN */
/********/
int main() {


    const int N = 2 << 18;
    cout << N;
    //const int N = 64;

    const int numIters = 50;

    thrust::host_vector<double>     h_v1(N);
    thrust::device_vector<double>     h_v2(N);
    thrust::host_vector<double>     h_v3(N);
    thrust::host_vector<double>     h_v4(N);
    thrust::device_vector<double>     h_v5(N);
    thrust::device_vector<double>     h_v6(N);

    printf("N = %d\n", N);

    {
        string name = "transform on host";
        checkpoint_timer timer{{name}};
        for (int k = 0; k < numIters; k++) {
            timer.set_startpoint(name);
            thrust::transform(thrust::host, h_v1.begin(), h_v1.end(), h_v1.begin(), rand_01<double>());
            timer.set_endpoint(name);
        }
    }

    {
        string name = "transform on device";
        checkpoint_timer timer{{name}};
        for (int k = 0; k < numIters; k++) {
            timer.set_startpoint(name);
            thrust::counting_iterator<unsigned int> index_sequence_begin(0);
            thrust::transform(index_sequence_begin,
                              index_sequence_begin + N,
                              h_v2.begin(),
                              prg(0.f, 1.f));
            timer.set_endpoint(name);
        }
    }

    {
        string name = "for each";
        checkpoint_timer timer{{name}};
        for (int k = 0; k < numIters; k++) {
            timer.set_startpoint(name);
            thrust::for_each(h_v3.begin(), h_v3.end(), rand_01_for_each<double>());
            timer.set_endpoint(name);
        }
    }

    //std::cout << "Values generated: " << std::endl;
    //for (int k = 0; k < N; k++)
    //  std::cout << h_v3[k] << " : ";
    //std::cout << std::endl;
    {
        string name = "generate";
        checkpoint_timer timer{{name}};
        for (int k = 0; k < numIters; k++) {
            timer.set_startpoint(name);
            thrust::generate(h_v4.begin(), h_v4.end(), rand_01_fcn<double>);
            timer.set_endpoint(name);
        }

    }

    {
        string name = "curand";
        checkpoint_timer timer{{name}};
        for (int k = 0; k < numIters; k++) {
            timer.set_startpoint(name);
            thrust::counting_iterator<unsigned int> index_sequence_begin(0);
            thrust::transform(index_sequence_begin, index_sequence_begin + N, h_v5.begin(), curand_slow());
            timer.set_endpoint(name);
        }
        cout << "Sample random values for " << name << endl;
        for(int i = 0; i < 20; i++) {
            cout << h_v5[i] << endl;
        }
        double mu = 0;
        double msd = 0;
        for(int i = 0; i < N; i++) {
            mu += h_v5[i];
            msd += h_v5[i] * h_v5[i];
        }
        cout << "mean: " << mu/N << endl;
        cout << "msd: " << msd/N << endl;


    }

    {
        string name = "curand fast";
        checkpoint_timer timer{{name}};
        thrust::device_vector<curandState> s1(N);
        auto sInit = thrust::make_zip_iterator(thrust::make_tuple(thrust::counting_iterator<int>(0), s1.begin()));

        thrust::for_each(sInit, sInit + N, curand_fast_setup());

        auto start = thrust::make_zip_iterator(thrust::make_tuple(h_v6.begin(), s1.begin()));
        auto end = thrust::make_zip_iterator(thrust::make_tuple(h_v6.end(), s1.end()));
        for (int k = 0; k < numIters; k++) {
            timer.set_startpoint(name);
            thrust::for_each(start, end, curand_fast());
            timer.set_endpoint(name);
        }

        cout << "Sample random values for " << name << endl;
        for(int i = 0; i < 20; i++) {
            cout << h_v5[i] << endl;
        }
        double mu = 0;
        double msd = 0;
        for(int i = 0; i < N; i++) {
            mu += h_v6[i];
            msd += h_v6[i] * h_v6[i];
        }
        cout << "mean: " << mu/N << endl;
        cout << "msd: " << msd/N << endl;


    }
    //std::cout << "Values generated: " << std::endl;
    //for (int k = 0; k < N; k++)
    //  std::cout << h_v4[k] << " : ";
    //std::cout << std::endl;

    //std::cout << "Values generated: " << std::endl;
    //for (int k = 0; k < N * 2; k++)
    //  std::cout << h_v[k] << " : ";
    //std::cout << std::endl;

    return 0;
}