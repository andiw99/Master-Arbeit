//
// Created by andi on 06.04.23.
//

#include "Testclass.h"
#include <iostream>

class parent {
public:
    virtual void foo() {
        std::cout << "parent::foo()" << std::endl;
    }
};

class child : public parent {
public:
    void foo() override {
        std::cout << "child::foo()" << std::endl;
    }
};

class pointer_container {
public:
    parent* b; // using a pointer to B
    pointer_container(class parent *b_val) {
        b = b_val;
    }
    pointer_container(){}
};


class container {
public:
    parent member;
    container(parent member) : member(member) {}

};

class Shape {
public:
    void draw() {
        std::cout << "Drawing a shape" << std::endl;
    };
};

class Circle : public Shape {
public:
    void draw() {
        std::cout << "Drawing a circle" << std::endl;
    }
};

class Square : public Shape {
public:
    void draw() {
        std::cout << "Drawing a square" << std::endl;
    }
};

int main() {
    child c;
    pointer_container a{&c};
    c.foo();
    a.b = &c; // assigning the address of c to the pointer to B in A
    a.b->foo(); // calls C::foo() due to dynamic binding

    parent* c2;

    c2.foo();
    pointer_container a2{&c2};
    a2.b->foo();

    child c3;
    container con(c3);
    con.member.foo();

    Shape shapes[] = {Circle(), Square() };

    for (auto shape : shapes) {
        shape.draw();
    }

    return 0;
}