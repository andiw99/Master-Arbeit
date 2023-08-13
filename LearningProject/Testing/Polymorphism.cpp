//
// Created by andi on 30.05.23.
//

#include "Polymorphism.h"
#include <iostream>
#include <chrono>
#include <stdio.h>
#include <vector>
#include <numeric>

using namespace std;

template <size_t lat_dim>
struct Parent {
public:
    struct Substruct {
    public:
        Substruct() {}
        void operator() () {
            cout << "Do i work? " << lat_dim << endl;
        }
    };

    struct Functor {
    public:
        Functor() {}
        virtual void operator()(){
            cout << "Dummy Functor is called" << endl;
        }
    };

    template <class FunctorType>
    void universalOperation(FunctorType functor) {
        functor();
    }

    void normalOp() {
        cout << "Just a normal fucking funciton" << endl;
    }

};

template <size_t lat_dim>
struct Child : public Parent<lat_dim> {
    using Substruct = typename Parent<lat_dim>::Substruct;
    using Functor = typename Parent<lat_dim>::Functor;
public:
    Child() {}

    void someMethod() {
        Substruct()();
    }

    struct ChildFunctor : public Functor {
    public:
        int id;
        ChildFunctor(int id) : id(id) {}


        void operator()(){
            cout << "Child Functor is called " << id <<  endl;
        }
    };

    void operator()() {
        ChildFunctor functor = ChildFunctor(6);
        this->universalOperation(functor);
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


class A{
public:
    virtual void run() {
        cout << "running the virutal function?" << endl;
    }
    void repeat() {
        this->run();
    }
};

class B : public A {
    void run() override {
        cout << "calling the right shit fam" << endl;
    }
};

int main(int argc, char* argv[]) {

    Child<5> child;

    child.someMethod();
    child();

    child.normalOp();

    size_t k = 0;
    unsigned int k2 = 0;
    k--;
    k2--;
    cout << k << endl;
    cout << k2 << endl;


    int count = 0;
    {
        timer timer;
        for(int i = 0; i < 2000000; i++) {
            if (count < i) {
                count++;
            } else {
                count--;
            }
        }

    }

    if (argc >= 2) {
        // Access and print the command-line arguments
        for (int i = 1; i < argc; ++i) {
            std::cout << "Argument " << i << ": " << argv[i] << std::endl;
        }
    } else {
        std::cout << "No arguments provided." << std::endl;
    }

    double fckn_testdouble = 0.01;
    printf("%f", fckn_testdouble);
    cout << "hello?" << endl;
    cout << time(NULL) % 10000 << endl;
#include <chrono>

// ...

    using namespace std::chrono;
    milliseconds ms = duration_cast< milliseconds >(
            system_clock::now().time_since_epoch()
    );

    cout << ms.count() % 100000 << endl;
    int a = 5;
    int b = 2;
    cout << "integer division" << endl;
    cout << a / b << endl;



    std::vector<int> input = {1, 2, 3, 4, 5};

    int sumOfSquares = std::transform_reduce(input.begin(), input.end(),

                                             0, // initial value for the reduction (sum)
                                              std::plus<int>(), // transformation (square)
                                             [](int b) { return b * b; }); // reduction (sum)
    std::vector<double> input2 = {1, 2, 3, 4, 5};

    double dsumOfSquares = std::transform_reduce(input2.begin(), input2.end(),

                                             0, // initial value for the reduction (sum)
                                             std::plus<int>(), // transformation (square)
                                             [](int b) { return b * b; }); // reduction (sum)

    std::cout << "Sum of squares: " << sumOfSquares << std::endl;

    B B_class{};
    B_class.repeat();

    return 0;

}
