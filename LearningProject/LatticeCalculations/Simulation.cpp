//
// Created by andi on 11.08.23.
//

#include "Simulation.h"


int main() {
    // possible Calcs: BinderCumulant, CorrLength, SecondMomentCorr
    vector<Calc> calcs = {Calc::SecondMomentCorr};
    for (auto calc : calcs) {
        cout << (int)calc << endl;
    }
    fs::path root = "../../Generated content/AA/AA Quench";
    simulation sim = simulation(root);
    sim.run(calcs);

    return 0;
}