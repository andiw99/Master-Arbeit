//
// Created by weitze73 on 09.03.23.
//


#include <array>
#include <iostream>
#include <boost/numeric/odeint.hpp>

template <size_t n>
std::array<std::array<int, n>, n> createQuadraticArray() {
    std::array<std::array<int, n>, n> arr{};
    for (size_t i = 0; i < n; i++) {
        for (size_t j = 0; j < n; j++) {
            arr[i][j] = 0;
        }
    }
    return arr;
}

// wir benötigen doch ein 3D array weil wir q und q' speichern müssen
// bzw müssen die einträge vom array 2D Vektoren sein
std::array<std::array<std::array<double, 2>, 200>, 200> createQuadraticArray2(int n) {
    std::array<std::array<std::array<double, 2>, 200>, 200> arr{};
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            // Anfangsbedingungen q' soll hier mal 1 sein
            arr[i][j][0] = 0;
            arr[i][j][1] = 1;
        }
    }
    return arr;
}

/* The type of container used to hold the state vector */
typedef std::array< double, 2 > state_type;

const double gam = 0.15;

// The rhs of x' = f(x)
// man übergibt vektor x = (q, q') und vektor x' = (q', q'')
// erste Gleichung ist nur x'[0] = x[1]
// zweite Gleichung ist dann die Gleichung für harmonischen Oszillator: q'' = -q - gam * q'
void harmonic_oscillator( const state_type &x , state_type &dxdt , const double /* t */ )
{
    dxdt[0] = x[1];
    dxdt[1] = -x[0] - gam*x[1];
}

int main() {
    // initialize parameters
    int n = 10;         // size of the lattice (one dimension)


    auto arr = createQuadraticArray2(n);
    for(int i = 0; i < 5; i++) {
        std::size_t steps = boost::numeric::odeint::integrate(harmonic_oscillator, arr[i][0], 0.0, 10.0, 0.1);
        std::cout << arr[i][0][0] << arr[i][0][1] << std::endl;
    }




    return 0;
}
