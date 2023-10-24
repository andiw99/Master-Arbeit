//
// Created by andi on 11.08.23.
//

#include "Simulation.h"
#include "configs.h"

int main() {
    typedef sin_functor transformation_functor;
    simulation sim = simulation<transformation_functor>(root);
    sim.run(calcs);

    return 0;
}