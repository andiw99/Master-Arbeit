//
// Created by andi on 06.04.23.
//

#include "Testclass.h"
#include <iostream>

class B {
public:
    virtual void foo() {
        std::cout << "B::foo()" << std::endl;
    }
};

class C : public B {
public:
    void foo() override {
        std::cout << "C::foo()" << std::endl;
    }
};

class A {
public:
    B* b; // using a pointer to B
    A(class B *b_val) {
        b = b_val;
    }
    A(){}
};

int main() {
    A a{};
    C c;
    a.b = &c; // assigning the address of c to the pointer to B in A
    a.b->foo(); // calls C::foo() due to dynamic binding

    C c2;
    A a2{&c2};
    a2.b->foo();
    return 0;
}