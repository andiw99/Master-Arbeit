//
// Created by andi on 11.08.23.
//

#include "Simulation.h"


int main() {
    // possible Calcs: BinderCumulant, CorrLength, SecondMomentCorr, StructFact
    vector<Calc> calcs = {Calc::BinderCumulant};
    for (auto calc : calcs) {
        cout << (int)calc << endl;
    }
    fs::path root = "../../../Generated content/Defense/eta=0.6/";
    simulation sim = simulation(root);
    sim.run(calcs);

    return 0;
}