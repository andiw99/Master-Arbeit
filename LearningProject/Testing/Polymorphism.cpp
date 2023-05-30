//
// Created by andi on 30.05.23.
//

#include "Polymorphism.h"
#include <iostream>

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

    void universalOperation(Functor *functor) {
        functor->operator()();
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
        ChildFunctor() {}


        void operator()(){
            cout << "Child Functor is called" << endl;
        }
    };

    void specifiedOperation() {
        ChildFunctor functor = ChildFunctor();
        this->universalOperation(&functor);
    }
};



int main() {

    Child<5> child;

    child.someMethod();
    child.specifiedOperation();

    return 0;
}
