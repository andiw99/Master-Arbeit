//
// Created by andi on 25.08.23.
//

#ifndef LEARNINGPROJECT_CONFIGS_H
#define LEARNINGPROJECT_CONFIGS_H

enum class Calc {
    BinderCumulant,
    CorrLength,
    SecondMomentCorr,
    StructFact,
};
// possible Calcs: BinderCumulant, CorrLength, SecondMomentCorr, StructFact
vector<Calc> calcs = {Calc::StructFact};
// fs::path root = "../../../Generated content/Trash/New/Overdamped Quenching 2/";
fs::path root = "../../../Generated content/AA/AA Quench Copy";

map<string, double> StructFactConfig {
        {"cell_L", 40}
};


#endif //LEARNINGPROJECT_CONFIGS_H
