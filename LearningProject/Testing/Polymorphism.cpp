//
// Created by andi on 30.05.23.
//

#include "Polymorphism.h"
#include <iostream>
#include <chrono>

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


    return 0;
}
