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
    QuenchProcess
};
// possible Calcs: BinderCumulant, CorrLength, SecondMomentCorr, StructFact
vector<Calc> calcs = {Calc::BinderCumulant, Calc::CorrLength, Calc::StructFact};
// vector<Calc> calcs = {Calc::StructFact, Calc::QuenchProcess};
//fs::path root = "../../../Generated content/Trash/New/Overdamped Quenching 2/";
fs::path root = "../../../Generated content/Testing Rectangular/Relaxation";

map<string, double> StructFactConfig {
        {"cell_L", 128},
        {"cutup", 8}
};

map<string, double> CorrLengthHandlerConfig {
        {"starting_k", 8},
        {"nr_Ls", 10}
};

map<string, double> BinderHandlerConfig {
        {"starting_k", 10},
        {"nr_Ls", 16}
};
map<string, double> QuenchProcessHandlerConfig {
        {"cell_L", 128},
        {"cutup", 4}
};

#endif //LEARNINGPROJECT_CONFIGS_H
