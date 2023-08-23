//
// Created by andi on 11.08.23.
//

#include "Simulation.h"


int main() {
    // possible Calcs: BinderCumulant, CorrLength, SecondMomentCorr
    vector<Calc> calcs = {Calc::BinderCumulant, Calc::CorrLength};
    for (auto calc : calcs) {
        cout << (int)calc << endl;
    }
    fs::path root = "../../../Generated content/Defense/Binder Cumulant/";
    simulation sim = simulation(root);
    sim.run(calcs);

    return 0;
}