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
// vector<Calc> calcs = {Calc::BinderCumulant, Calc::CorrLength};
vector<Calc> calcs = {Calc::StructFact};
// fs::path root = "../../../Generated content/Trash/New/Overdamped Quenching 2/";
fs::path root = "../../../Generated content/Defense/Real Quench";

map<string, double> StructFactConfig {
        {"cell_L", 64}
};

map<string, double> CorrLengthHandlerConfig {
        {"starting_k", 4},
        {"nr_Ls", 18}
};

map<string, double> BinderHandlerConfig {
        {"starting_k", 4},
        {"nr_Ls", 18}
};


#endif //LEARNINGPROJECT_CONFIGS_H
