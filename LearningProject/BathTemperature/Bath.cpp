//
// Created by weitze73 on 12.04.23.
//

#include "Bath.h"
#include "../Helpfunctions and Classes.h"
#include "../Systems.h"
#include "../Solvers.h"



int main() {
    // So the plan is just to change the temperature dependency from eps(T) to theta(T)
    // We should be able to do this by just defining another lattice_system class in Functions and Classes.h
    // which has the desired behaviour
    // we then have to rewrite the search grid function to take in the kind of system that should be tried out by the
    // parameters

    string storage_root = "../../Generated content/Cooling Bath/Test/";

    // we need to make sure that we cross the phase transition
    // so we have to select our parameters wisely
    // We choose a starting Temperature and an ending temperature
    // at best they are equally distant to the phase transition
    // the problem is that we dont exactly know where the phase transition takes place
    // i would say we give T_start and T_end and tau and the maximum stepsize and a desired number of steps
    // then we calculate the stepsize from the number of steps


    // TODO maybe replace dt with tmax
    vector<double> eta_values{5};
    vector<double> T_start{70};
    vector<double> T_end{20};
    vector<int> steps_values{100000};
    vector<int> n_values{100};
    vector<double> alpha_values{5};
    vector<double> beta_values{10};
    vector<double> J_values{50};
    vector<double> tau_values{10};

    double max_dt = 0.005;

    search_grid_bath<cooling_bath>(eta_values, T_start, T_end, steps_values, n_values, alpha_values, beta_values,
                                   J_values, tau_values, storage_root, max_dt);

//    search_grid<cooling_bath>(eta_values, T_values, dt_values, steps_values, n_values, alpha_values, beta_values, J_values, tau_values,
  //                                       storage_root);

    return 0;
}
