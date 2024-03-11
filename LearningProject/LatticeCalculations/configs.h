//
// Created by andi on 25.08.23.
//

#ifndef LEARNINGPROJECT_CONFIGS_H
#define LEARNINGPROJECT_CONFIGS_H


enum class Calc {
    BinderCumulant,
    SurBinderCumulant,
    SurBinderCumulantSilicon,
    CorrLength,
    CorrLengthXY,
    SurCorrLengthXY,
    SecondMomentCorr,
    StructFact,
    StructFactXY,
    QuenchProcess
};
// possible Calcs: BinderCumulant, CorrLength, SecondMomentCorr, StructFact
//vector<Calc> calcs = {Calc::BinderCumulant, Calc::CorrLengthXY, Calc::StructFactXY};
//vector<Calc> calcs = {Calc::SurBinderCumulantSilicon, Calc::SurCorrLengthXY};
//vector<Calc> calcs = {Calc::SurCorrLengthXY};
vector<Calc> calcs = {Calc::StructFactXY};
//vector<Calc> calcs = {Calc::StructFactXY, Calc::QuenchProcess};
//fs::path root = "../../../Generated content/Trash/New/Overdamped Quenching 2/";
fs::path root = "../../../Generated content/Silicon/Subsystems/Suite/L_xi/scan-more-flips-more-vals/0.4161791450287818/Tc/120";

map<string, double> StructFactConfig {
        {"cell_L", 0},
        {"cutup", 1},
        {"subsystems", 1}
};

map<string, double> CorrLengthHandlerConfig {
        {"starting_k", 10},
        {"nr_Ls", 20},
        {"min_Lx", 5},
        {"max_Lx", 40},
        {"cut_zero_impuls", 0}
};

map<string, double> BinderHandlerConfig {
        {"starting_k", 8},
        {"nr_Ls", 14},
        {"min_L", 2},
        {"max_L", 40}
};
map<string, double> QuenchProcessHandlerConfig {
        {"cell_L", 128},
        {"cutup", 6}
};

#endif //LEARNINGPROJECT_CONFIGS_H
