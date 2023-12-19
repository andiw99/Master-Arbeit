//
// Created by andi on 30.05.23.
//

#include "Polymorphism.h"
#include <iostream>
#include <chrono>
#include <stdio.h>
#include <vector>
#include <numeric>
#include <boost/asio/ip/host_name.hpp>
#include "../Header/Helpfunctions and Classes.h"

using namespace std;


struct ParentClass {
public:
    struct parentFunctor {
        void operator()(double x) {
            cout << "parent Functor" << endl;
        }
    };
};

struct SpecializedParentClass: virtual public ParentClass {
public:
    struct parentFunctor {
        void operator()(double x) {
            cout << "specialized Functor called" << endl;
        }
    };
};

struct directChildClass : public ParentClass {
public:
    void operation(double x) {
        parentFunctor()(x);
    }
};

struct specializedChildClass : public SpecializedParentClass {
public:
    void operation(double x) {
        parentFunctor()(x);
    }
};


struct usedClass: public specializedChildClass {
};

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
    void run() {
        cout << "running the virutal function?" << endl;
    }
    void repeat() {
        this->run();
    }
};

class B : public A {
public:
    B(int a) {
        cout << "B constructor called with " << a << endl;
    }
    void run() {
        cout << "calling the right shit fam" << endl;
    }
};

class b : public B {
public:
    b(int x, int y) : B(x) {
        cout << "b constructor called with " << x << "   " << y << endl;
    }
};

class c : public B {
public:
    c(int x, double y) : B(x) {
        cout << "b constructor called with " << x << "   " << y << endl;
    }
};


class C : public A {
public:
    C(double a) {
        cout << "C constructor called with " << a << endl;
    }
    void run() {
        cout << "calling the right shit fam" << endl;
    }
};


class ParentTemplateA {
public:
    template <class value_type>
    void universalOP(value_type x) {
        cout << "universal parentTemplate method" << endl;
    }
};

class DerivedTemplateB : virtual public ParentTemplateA {
public:
    virtual void specificOP() {
        cout << "DerivedTemplate B specific OP" << endl;
        universalOP(1.0);
    }
};

class DerivedTemplateC : virtual public ParentTemplateA {
public:
    template <class value_type>
    void universalOP(value_type x) {
        cout << "universal OP override" << endl;
    }
};

class DerivedTemplateBC : public DerivedTemplateB, public DerivedTemplateC {
};

class ParentA {
public:
    virtual void universalOP(double x) {
        cout << "universal parent method" << endl;
    }
};

class DerivedB : virtual public ParentA {
public:
    virtual void specificOP() {
        cout << "Derived B specific OP" << endl;
        universalOP(1.0);
    }
};

class DerivedC : virtual public ParentA {
public:
    void universalOP(double x) override {
        cout << "universal OP override" << endl;
    }
};

class DerivedBC : public DerivedB, public DerivedC {
};

class Container {
public:
    vector<A*> observers;

    template <class obsver>
    void register_observer(obsver* obs) {
        cout << "Registering observers" << endl;
        obs->run();
        observers.push_back(obs);
        observers[0]->run();
    }
};

class DerivedContainer: public Container {

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
        for (int i = 0; i < 2000000; i++) {
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
    milliseconds ms = duration_cast<milliseconds>(
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

    A *base = new A();
    DerivedContainer container;
    container.register_observer(base);

    int rows = 3;
    int cols = 4;

    double **matrix = new double *[rows];
    for (int i = 0; i < rows; ++i) {
        matrix[i] = new double[cols];
    }

    // Initializing values
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            matrix[i][j] = i * cols + j;
        }
    }

    // Accessing and printing values
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < cols; ++j) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }

    // Deallocating memory
    for (int i = 0; i < rows; ++i) {
        delete[] matrix[i];
    }
    delete[] matrix;


    DerivedBC derived{};
    derived.specificOP();

    DerivedTemplateBC derivedTemplated{};
    derivedTemplated.specificOP();


    const auto host_name = boost::asio::ip::host_name();
    cout << host_name << endl;

    auto timepoint = chrono::system_clock::now().time_since_epoch();
    auto hour = duration_cast<hours>(timepoint) % 24;
    auto minute = duration_cast<minutes>(timepoint) % 60;
    auto second = duration_cast<seconds>(timepoint) % 60;
    auto millisecond = duration_cast<milliseconds>(timepoint) % 1000;

    cout << "current time: " << hour.count() << ":" << minute.count() << ":" << second.count() << ":"
         << millisecond.count() << endl;
    cout << "current time: " << get_current_time() << endl;

    std::vector<double> data = {1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8, 9.9};

    // Define the range of the subset using iterators
    auto subset_begin = data.begin() + 2;  // Iterator pointing to the third element
    auto subset_end = data.begin() + 6;    // Iterator pointing to the seventh element

    // Create a set of pointers to the subset

    vector<int> inds = {3, 6, 7};

    std::vector<double*> subset_pointers(3);
    int j = 0;
    for (auto ind : inds) {
        subset_pointers[j] = (&(*(data.begin() + ind)));
        j++;
    }

    // Print the values in the subset using the pointers
    std::cout << "Subset values: ";
    for (int i= 0; i < subset_pointers.size(); i++) {
        std::cout << *(subset_pointers[i]) << " ";
        *(subset_pointers[i]) = 1.0;
    }

    cout << endl << "Values:" << endl;
    for(auto thing : data) {
        cout << thing << "  " ;
    }

    std::cout << std::endl;

    usedClass childClass;
    childClass.operation(2.0);

}
