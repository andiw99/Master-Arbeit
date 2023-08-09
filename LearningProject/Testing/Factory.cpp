//
// Created by andi on 09.08.23.
//

#include "Factory.h"
#include <iostream>
using namespace std;

template <int N>
class MyTemplateClass {
public:
    MyTemplateClass() {
        // Constructor code
        cout << "Constructor of Class1 is called" << endl;
    }
    // Other class methods
};

template <int N>
class MyTemplateClass2 {
public:
    MyTemplateClass2() {
        // Constructor code
        cout << "Constructor of Class2 is called" << endl;
    }
    // Other class methods
};


template <int N>
class MyTemplateClassFactory {
public:
        template <template<size_t> class T>
        static T<N> create() {
            return T<N>();
        }

/*        template<> MyTemplateClass<N> create() {
            cout << "Create for 1 is called" << endl;
        }
        template<> MyTemplateClass2<N> create() {
            cout << "Create2 is called" << endl;
        }*/
};


template<>
template<>
MyTemplateClass<5> MyTemplateClassFactory<5>::create<MyTemplateClass2<5>>() {

}

template< size_t lat_dim, template<size_t> class templateClass>
void routine() {
        templateClass<lat_dim> instance = MyTemplateClassFactory<lat_dim>::template create<templateClass>();
        }


int main() {
    MyTemplateClass<5> instance1 =   MyTemplateClassFactory<5>::create<MyTemplateClass>();
    MyTemplateClass2<10> instance2 = MyTemplateClassFactory<10>::create<MyTemplateClass2>();

    // Use instances

    routine<10, MyTemplateClass2>();

    return 0;
}
