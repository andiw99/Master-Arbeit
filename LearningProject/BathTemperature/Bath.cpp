//
// Created by weitze73 on 12.04.23.
//

#include "Bath.h"
#include "../Helpfunctions and Classes.h"
#include "../Systems.h"
#include "../Solvers.h"

double ind_value(vector<double> paras, int ind) {
    return paras[ind % paras.size()];
}


int ind_value(vector<int> paras, int ind) {
    return paras[ind % paras.size()];
}

template <typename sys>
void search_grid_bath(vector<double> eta_values, vector<double> T_start, vector<double> T_end,
                      vector<int> steps_values, vector<int> n_values,
                      vector<double> alpha_values, vector<double> beta_values, vector<double> J_values,
                      vector<double> tau_values,
                      string storage_root, double starting_t = 0, double max_dt = 0.005) {
    // first we find the vector with the largest size
    int configs = max({eta_values.size(), T_start.size(), T_end.size(), steps_values.size(), n_values.size(),
                   alpha_values.size(), beta_values.size(), J_values.size(), tau_values.size()});
    cout << configs << " configs" << endl;

    for(int i = 0; i < configs; i++) {
        // we expect that the i-th entry of every vektor belongs to one config
        // but if there are not enough values we just use the first one or random?
        // we know calculate the step size for this run
        int steps = steps_values[i % steps_values.size()];
        double tau = ind_value(tau_values, i);
        double dt = (ind_value(T_start, i)- ind_value(T_end, i)) * tau
                / steps;
        // check whether dt is to small
        if(dt > max_dt) {
            dt = max_dt;
            // set the steps for this step
            steps = (T_start[i % T_start.size()] - T_end[i % T_end.size()]) * tau / dt;
            cout << "Had to rescale steps. New number of steps: " << steps << endl;
        }

        init_and_run<sys>(ind_value(eta_values, i), ind_value(T_start, i), dt, steps,
                          ind_value(n_values, i), ind_value(alpha_values, i),
                          ind_value(beta_values, i), ind_value(J_values, i), tau,
                          storage_root, starting_t);
        cout << "run " << i+1 << " / " << configs << endl;
    }
}



int main() {
    // So the plan is just to change the temperature dependency from eps(T) to theta(T)
    // We should be able to do this by just defining another lattice_system class in Functions and Classes.h
    // which has the desired behaviour
    // we then have to rewrite the search grid function to take in the kind of system that should be tried out by the
    // parameters

    string storage_root = "../../Generated content/Cooling Bath/Comparison/";

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
    vector<int> n_values{50};
    vector<double> alpha_values{5};
    vector<double> beta_values{10};
    vector<double> J_values{50};
    vector<double> tau_values{1};

    double max_dt = 0.005;

    search_grid_bath<cooling_bath>(eta_values, T_start, T_end, steps_values, n_values, alpha_values, beta_values,
                                   J_values, tau_values, storage_root, max_dt);

//    search_grid<cooling_bath>(eta_values, T_values, dt_values, steps_values, n_values, alpha_values, beta_values, J_values, tau_values,
  //                                       storage_root);

    return 0;
}
