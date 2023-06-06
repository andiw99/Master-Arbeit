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
    __host__ T operator()(T& VecElem) const { return (T)rand() / RAND_MAX; }
};

template<typename T>
struct rand_01_for_each {
    __host__ void operator()(T& VecElem) const { VecElem = (T)rand() / RAND_MAX; }
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

    thrust::device_vector<double>     h_v1(N);
    thrust::device_vector<double>     h_v2(N);
    thrust::device_vector<double>     h_v3(N);
    thrust::device_vector<double>     h_v4(N);

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