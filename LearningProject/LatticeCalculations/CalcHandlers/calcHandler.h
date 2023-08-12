//
// Created by andi on 11.08.23.
//

#ifndef LEARNINGPROJECT_CALCHANDLER_CPP
#define LEARNINGPROJECT_CALCHANDLER_CPP
#include "../../Header/Helpfunctions and Classes.h"


class calcHandler {
public:
    // Base class for the Handlers of the different calculations
    // Every Calculation will have a ofstream
    ofstream calcFile;
    fs::path root;      // probably very useful to have this as a class variable since some Handlers have to
    // write in different directories
    size_t lat_dim = 0;
    calcHandler(){
        root = "./";
    }
    explicit calcHandler(const fs::path& root): root(root){
    }
    calcHandler(const fs::path& root, size_t lat_dim): root(root), lat_dim(lat_dim){
    }
    virtual void pre_routine() {}           // some initialization to do before starting to iterate over the settings
    virtual void setting_pre_routine(fs::path setting_path) {
        // we adapt that before the calculation for the setting, the Handler gets the path to do all the
        // things that might be unusual, like reading tau
    }   // something to do for every setting before iterating over the realizations
    virtual void realization_routine(vector<double> &lat_q, double T, double t) {}   // something to do for every realization of every setting
    virtual void setting_post_routine(){}   // something to do after we iterated over all realizations
    virtual void post_routine() {}          // something to do after we are done with every setting

    void set_lat_dim(size_t lattice_dim) {
        lat_dim = lattice_dim;
    }
};



#endif //LEARNINGPROJECT_CALCHANDLER_CPP
