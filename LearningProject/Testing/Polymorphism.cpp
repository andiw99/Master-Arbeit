//
// Created by andi on 30.05.23.
//

#include "Polymorphism.h"
#include <iostream>

using namespace std;

struct Parent {
public:
    struct Substruct {
    public:
        Substruct() {}
        void operator() () {
            cout << "Do i work?" << endl;
        }
    };
};

struct Child : public Parent {
public:
    Child() {}

    void someMethod() {
        Substruct()();
    }
};



int main() {

    Child child = Child();

    child.someMethod();

    return 0;
}
