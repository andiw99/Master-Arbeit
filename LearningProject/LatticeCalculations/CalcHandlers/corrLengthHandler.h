//
// Created by andi on 11.08.23.
//

#ifndef LEARNINGPROJECT_CORRLENGTHHANDLER_H
#define LEARNINGPROJECT_CORRLENGTHHANDLER_H

#include "../../Header/Helpfunctions and Classes.h"
#include "calcHandler.h"

class CorrLengthHandler : public calcHandler {
    using calcHandler::calcHandler;
public:

    void pre_routine() override {
        cout << "Calling CorrLength pre routine" << endl;
    }

    void realization_routine(vector<double> &lat_q, double T, double t) override {
        cout << "Calling realization routine of corr length" << endl;
    }

    void post_routine() override {
        cout << "Calling CorrLengthHandler post routine" << endl;
    }
};

#endif //LEARNINGPROJECT_CORRLENGTHHANDLER_H
