//
// Created by andi on 11.08.23.
//

#ifndef LEARNINGPROJECT_SIMULATION_H
#define LEARNINGPROJECT_SIMULATION_H

#include "CalcHandlers/BinderHandler.h"
#include "CalcHandlers/corrLengthHandler.h"
#include "CalcHandlers/secondCorrLenghtHandler.h"
#include "CalcHandlers/calcHandler.h"
#include "CalcHandlers/StructFactHandler.h"
#include "CalcHandlers/QuenchProcessHandler.h"




calcHandler* create(Calc calc, const fs::path& root){
    // i think we sadly need to return a pointer to the created object?
    // Never really understood that
    if (calc == Calc::BinderCumulant) {
        return new BinderHandler(root);
    } else if (calc == Calc::SurBinderCumulant) {
        return new SurBinderHandler(root);
    } else if (calc == Calc::SurBinderCumulantSilicon) {
        return new SurBinderHandlerSilicon(root);
    } else if (calc == Calc::CorrLength) {
        return new CorrLengthHandler(root);
    } else if (calc == Calc::CorrLengthXY) {
        return new CorrLengthHandlerXY(root);
    } else if (calc == Calc::SurCorrLengthXY) {
        return new SurCorrLengthHandlerXY(root);
    } else if (calc == Calc::SecondMomentCorr) {
        return new secondCorrLenghtHandler(root);
    } else if (calc == Calc::StructFact) {
        return new StructFactHandler(root);
    } else if (calc == Calc::StructFactXY) {
        return new StructFactHandlerXY(root);
    } else if (calc == Calc::QuenchProcess) {
        return new QuenchProcessHandler(root);
    }
    return nullptr;
}



struct unity_functor {
    template <class value_type>
    value_type operator()(value_type x) {
        return x;
    }
};

calcHandler* create(Calc calc){
    return create(calc, "./");
}

using namespace std;
namespace fs = std::filesystem;
using namespace fs;

template <class transform_functor>
class simulation {
    // We probably need many class members and we maybe should split the functions to different handlers
    fs::path root;
    size_t dim_size_x;
    size_t dim_size_y;
    bool chessTrafo = true;
    vector<fs::path> root_directories;
    vector<calcHandler*> handlers = {};
public:
    simulation(const fs::path& root): root(root) {
        // First thing we should get done is look for the lattice dim
        fs::path txt_file = findFirstTxtFile(root);
        cout << "found txt file: " << txt_file << endl;
        dim_size_x = (size_t)extractValueFromTxt(txt_file, "dim_size_x");
        dim_size_y = (size_t)extractValueFromTxt(txt_file, "dim_size_y");
        cout << "dim_size_x = " << dim_size_x << endl;
        cout << "dim_size_y = " << dim_size_y << endl;
    }

    simulation(const fs::path& root, bool chessTrafo): root(root), chessTrafo(chessTrafo) {
        // TODO duplicate code, maybe look into delegating constructor? Seems to be not to easy
    }

    void run(vector<Calc> calcs) {
        root_directories = list_dir_paths(root, false);

        // run is supposed to iterate over every csv file and do all operations that are requested
        // we need some system that tells run which values should be calculated.
        // vector of strings would work but thats sketchy i think
        // Okay it is really not to easy to solve, for now I think we will just use a
        // create function that consits of chained if else statements
        for(auto calc : calcs) {
            // We add the pointer to the new handler to our vector of pointers.
            // some of the handlers need parameters, how do i hand them over?
            handlers.push_back(create(calc, root));
        }
        cout << "Created Handlers" << endl;
        // now we have all handlers in the handler vector

        // I just thougt of: We could speed up the process even more if we
        // iterated over every lattice only once, but that might be difficult to implement
        // But i mean it should be possible to have a vector of objects, the objects itself
        // being handlers for the different calculations. Then in every stage of the run
        // we iterate through this vector and call the appropriate method.
        // I sense it could get complicated if we have the calculations for the subsystems etc.
        // but lets just start and see what will happen

        // Okay, some general things that we need to do when calling run?
        // cant think of any for now...

        // okay initialization calls for the handlers
        for (auto handler : handlers) {
            // Or get it into the pre routine?
            handler->pre_routine();
        }
        for(auto directory : root_directories) {
            // We atm do not now what kind of directory it is, could already be a setting directory (temperature or tau)
            // but it could also be one over
            directory_routine(directory);
        }

        for (auto handler : handlers) {
            handler->post_routine();
        }
    }

    void directory_routine(path& directory_path) {
        // We somehow need to find out whether it is a setting routine, lets just check whether there is a .csv file in
        // the directory

        if(isCSVFileInCurrentDirectory(directory_path)) {
            setting_routine(directory_path);
        } else {
            for (auto handler : handlers) {
                handler->directory_pre_routine(directory_path);
            }
            auto directories = list_dir_paths(directory_path);
            sort(directories.begin(), directories.end());
            for (auto path : directories) {
                cout << "Before inner directory routine" << endl;
                cout << path << endl;
                directory_routine(path);
            }
        }
    }

    void setting_routine(path &setting_path) {
        fs::path txt_file = findFirstTxtFile(setting_path);
        cout << "found txt file (setting routine): " << txt_file << endl;
        dim_size_x = (size_t)extractValueFromTxt(txt_file, "dim_size_x");
        dim_size_y = (size_t)extractValueFromTxt(txt_file, "dim_size_y");
        for (auto handler : handlers) {
            handler->set_dims(dim_size_x, dim_size_y);
            handler->setting_pre_routine(setting_path);
        }
        // ready to cycle over realizations?
        vector<path> csv_files = list_csv_files(setting_path);
        txt_file = findFirstTxtFile(setting_path);
        double J = extractValueFromTxt(txt_file, "J");
        if(J < 0) {
            chessTrafo = true;
        } else {
            chessTrafo = false;
        }
        for(auto csv_path : csv_files) {
            realization_routine(csv_path);

        }
        for (auto handler : handlers) {
            handler->setting_post_routine();
        }
    }

    void realization_routine(path &csv_path) {// find out if chess_trafo should be true or not

        ifstream file = safe_read(csv_path, false);
        double T, t;
        // TODO okay we only read in the last line here, so it doesnt work for the quench.process
// plot atm, but i guess we could fix that in the future?
        cout << "probably here?" << endl;
        auto lat_q = readDoubleValuesAt(file, -1,  T, t);
        //cout << "T = " << T << "   t = " << t << endl;
        if(chessTrafo) {
            chess_trafo_rectangular(lat_q, dim_size_x);
        }
        transform(lat_q.cbegin(), lat_q.cend(), lat_q.begin(), transform_functor());
        // calling the calcs... it would be even more efficient if I had to cycle over the
// lattice only once, but that would be really complicated to implement
        for (auto handler : handlers) {
            // for sure the realization routine needs the lattice
            handler->realization_routine(lat_q, T, t);
        }
    }
};

template <class transform_functor>
class Subsimulation : public simulation<transform_functor> {
    // We just start a simulation for every

};

#endif //LEARNINGPROJECT_SIMULATION_H
