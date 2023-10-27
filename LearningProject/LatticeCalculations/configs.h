//
// Created by andi on 25.08.23.
//

#ifndef LEARNINGPROJECT_CONFIGS_H
#define LEARNINGPROJECT_CONFIGS_H


enum class Calc {
    BinderCumulant,
    CorrLength,
    CorrLengthXY,
    SecondMomentCorr,
    StructFact,
    StructFactXY,
    QuenchProcess
};
// possible Calcs: BinderCumulant, CorrLength, SecondMomentCorr, StructFact
vector<Calc> calcs = {Calc::BinderCumulant, Calc::CorrLengthXY, Calc::StructFactXY};
// vector<Calc> calcs = {Calc::StructFact, Calc::QuenchProcess};
//fs::path root = "../../../Generated content/Trash/New/Overdamped Quenching 2/";
fs::path root = "../../../Generated content/XY/XY Peak";

map<string, double> StructFactConfig {
        {"cell_L", 128},
        {"cutup", 8}
};

map<string, double> CorrLengthHandlerConfig {
        {"starting_k", 10},
        {"nr_Ls", 20}
};

map<string, double> BinderHandlerConfig {
        {"starting_k", 8},
        {"nr_Ls", 14},
        {"min_L", 3},
        {"max_L", 50}
};
map<string, double> QuenchProcessHandlerConfig {
        {"cell_L", 128},
        {"cutup", 6}
};

#endif //LEARNINGPROJECT_CONFIGS_H
