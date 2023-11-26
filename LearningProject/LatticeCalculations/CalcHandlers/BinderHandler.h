//
// Created by andi on 11.08.23.
//

#ifndef LEARNINGPROJECT_BINDERHANDLER_CPP
#define LEARNINGPROJECT_BINDERHANDLER_CPP

#include "../../Header/Helpfunctions and Classes.h"
#include <map>
#include <boost/range/adaptor/reversed.hpp>
#include "calcHandler.h"


class BinderHandler : public calcHandler {
public:
    int starting_k = 10;                 // should find way to initialize them in create? Probably not
    int nr_Ls = 20; // to comfortable to always have to change them in the class
    int min_L, max_L;
    vector<int> L_vec;
    map<int, vector<pair<double, double>>> m_map = {};
    double T;
    double x_y_factor = 1;

    BinderHandler(fs::path root): calcHandler(root) {
        starting_k = (int) BinderHandlerConfig["starting_k"];
        nr_Ls = (int) BinderHandlerConfig["nr_Ls"];
        min_L = (int) BinderHandlerConfig["min_L"];
        max_L = (int) BinderHandlerConfig["max_L"];

    }

    void pre_routine() override {
        cout << "Calling BinderHandler pre routine" << endl;
        calcFile.open(root / "binder.cumulants");
        // Generating the vector of Ls
        // Okay i think we have to rewrite this entirely since we want to cut up into rectangles (or do we want to?)
        // TODO think about whether you want to calculate the binder cumulant in rectangular cells, for now i dont think so
        // L_vec = generate_L(starting_k, dim_size_x * dim_size_y, nr_Ls);
        nr_Ls = max_L - min_L + 1;
        vector<int> Ls(nr_Ls);
        iota(begin(Ls), end(Ls), min_L);
        L_vec = Ls;
        init_file_header(calcFile, L_vec);

    }

    void directory_pre_routine(path directory_path) override {
        nr_Ls = max_L - min_L + 1;
        vector<int> Ls(nr_Ls);
        iota(begin(Ls), end(Ls), min_L);
        L_vec = Ls;
        init_file_header(calcFile, L_vec);
        fs::path txt_file = findFirstTxtFile(directory_path);
        x_y_factor = extractValueFromTxt(txt_file, "x_y_factor");
    }

    void setting_pre_routine(fs::path setting_path) override {
        // we need a map that maps the system size to the vector containing all the m values
        // it has to be reset for every setting / temperature
        m_map = {};
    }

    void realization_routine(vector<double> &lat_q, double temp, double t) override {
        // for every subsystem size L we need to iterate over the lattice and extract all m's
        for(int Lx : L_vec) {
            int Ly = (int) (Lx * x_y_factor);
            int nr_cells = (int)(dim_size_x * dim_size_y / (Lx * Ly));
            for(int cell_nr = 0; cell_nr < nr_cells; cell_nr++) {
                // now we are in the ith cell and can start calculating the local magnetization powers
                // it might be easiest and maybe also the most efficient to extract the spins of the cell
                // out of the large lattice
                vector<double> q_cell = vector<double>(Lx * Ly, 0);
                extract_cell(lat_q, q_cell, cell_nr, Lx, Ly, dim_size_x);
                pair<double, double> m_L = calc_m(q_cell);
                m_map[Lx].push_back(m_L);
            }
        }
        T = temp;
    }

    virtual pair<double, double> calc_m(vector<double>& q_cell) {
        pair<double, double> m;
        m.first = transform_reduce(q_cell.begin(), q_cell.end(),  0.0, plus<double>(), cos_functor())  / (double)q_cell.size();
        m.second = transform_reduce(q_cell.begin(), q_cell.end(),  0.0, plus<double>(), sin_functor()) / (double)q_cell.size();
        return m;
    }

    void setting_post_routine() override {
        // setting is in this case the temperature, so we write the temp to file after iterating over all realizations
        // I need the tau instead of the t here...
        calcFile << T;
        // now we have the full m_map for the temperature, leaves to calculate the binder cumulant aswell as errors
        vector<double> cum, cum_error;
        calc_cum(cum, cum_error);
        write_cum(cum, cum_error);
        calcFile << endl;

    }

    void post_routine() override {
        cout << "Calling BinderHandler post routine" << endl;
    }



    vector<int> generate_L(int starting_k, int n, int nr_Ls) {
        vector<int> L = {(int)sqrt(n) / (2*starting_k)};
        for (int k = starting_k + 1; k < starting_k + nr_Ls; k++) {
            int L_val = (int)sqrt(n) / (2*k);
            cout << L.back() - L_val << endl;
            if (L.back() - L_val < 2) {
                L_val = L.back() - 2;
            }
            L.push_back(L_val);
        }
        return L;
    }

    void init_file_header(ofstream& file, const vector<int>& L_vec) {
        file << "T";
        for (int L : L_vec) {
            file << "," << L << "," << L << "_err";
        }
        file << endl;
    }



    void calc_cum (vector<double> &cum_vec, vector<double>& cum_error_vec) {
        for(auto L_m_vec_pair : m_map){
            double m_L2 = std::transform_reduce(L_m_vec_pair.second.begin(), L_m_vec_pair.second.end(),
                                                0.0, // initial value for the reduction (sum)
                                                std::plus<double>(), // transformation (square)
                                                [](::pair<double, double> m) -> double { return squared(m); });
            m_L2 /= L_m_vec_pair.second.size();
            double m_L2_err = std::transform_reduce(L_m_vec_pair.second.begin(), L_m_vec_pair.second.end(),
                                                    0.0, // initial value for the reduction (sum)
                                                    std::plus<double>(), // transformation (square)
                                                    [&m_L2](::pair<double, double> m) { return pow((squared(m) - m_L2), 2)  ; });
            m_L2_err /= (double)pow(L_m_vec_pair.second.size(), 2);
            m_L2_err = sqrt(m_L2_err);
            double m_L4 = std::transform_reduce(L_m_vec_pair.second.begin(), L_m_vec_pair.second.end(),
                                                0.0, // initial value for the reduction (sum)
                                                std::plus<>(), // transformation (square)
                                                [](::pair<double, double> m) { return (pow( squared(m), 2)); });
            m_L4 /= L_m_vec_pair.second.size();
            double m_L4_err = std::transform_reduce(L_m_vec_pair.second.begin(), L_m_vec_pair.second.end(),
                                                    0.0, // initial value for the reduction (sum)
                                                    std::plus<>(), // transformation (square)
                                                    [&m_L4](::pair<double, double> m) { return pow(pow(squared(m), 2) - m_L4, 2); });
            m_L4_err /= (double)pow(L_m_vec_pair.second.size(), 2);
            m_L4_err = sqrt(m_L4_err);
            double cum = m_L4 / (m_L2 * m_L2);
            double cum_error = sqrt(pow(1 / m_L2 / m_L2 * m_L4_err, 2) + pow(2 * m_L4 / pow(m_L2, 3) * m_L2_err, 2));
            cum_vec.push_back(cum);
            cum_error_vec.push_back(cum_error);
        }
    }

    void write_cum(vector<double> &cum_vec, vector<double>& cum_error_vec) {
        for(int i = 0; i < cum_vec.size(); i++) {
            calcFile << "," << cum_vec[i] << "," << cum_error_vec[i];
        }
    }
};

#define p_XY 2.57
class BinderHandlerSilicon : virtual public BinderHandler {
public:
    pair<double, double> calc_m(vector<double>& q_cell) override{
        pair<double, double> m;
        m.first = transform_reduce(q_cell.begin(), q_cell.end(),  0.0, plus<double>(), sin_functor(p_XY / 2.0))  / (double)q_cell.size();
        m.second = 0;
        return m;
    }
    using BinderHandler::BinderHandler;
};

class SurBinderHandler : virtual public BinderHandler {
    using BinderHandler::BinderHandler;
    map<size_t, vector<tuple<double, double, double>>> size_T_cum_map;
    size_t subsystem_Lx;
    void pre_routine() override {
        cout << "Calling SurBinderHandler pre routine" << endl;
        calcFile.open(root / "binder.cumulants");
    }

    void directory_pre_routine(path directory_path) override {
        fs::path txt_file = findFirstTxtFile(directory_path);
        x_y_factor = extractValueFromTxt(txt_file, "x_y_factor");
        subsystem_Lx = (size_t)extractValueFromTxt(txt_file, "subsystem_Lx");

        cout << "calling directory pre routine" << endl;
        cout << "subsystem_Lx = " << subsystem_Lx << endl;
        L_vec = {(int)subsystem_Lx};
        // if you now go through with the setting routines and stuff, you will be left with an m_map with one size
        // and the block m's
        size_T_cum_map[subsystem_Lx] = vector<tuple<double, double, double>>{};
    }


    void setting_post_routine() override {
        // setting is in this case the temperature, so we write the temp to file after iterating over all realizations
        // I need the tau instead of the t here...
        // now we have the full m_map for the temperature, leaves to calculate the binder cumulant aswell as errors
        vector<double> cum, cum_error;
        cout << "T = " << T << "   " << "subsystem_size = " << subsystem_Lx << endl;
        print_pair_vector(m_map[subsystem_Lx]);
        calc_cum(cum, cum_error);
        // we insert into the vector of size_T_cum_map for the current size, which is dim_size_x the triplet with the
        // current temperature, the cumulant and stuff.
        cout << "Setting post routine for T = " << T  << endl;
        size_T_cum_map[subsystem_Lx].push_back(tuple<double, double, double>(T, cum[0], cum_error[0]));
        cout << "subsystem size = " << subsystem_Lx << endl;
        for(auto tuples : size_T_cum_map[subsystem_Lx]) {
            cout << "(" << get<0>(tuples) << ", " << get<1>(tuples) << ", " << get<2>(tuples) << ")" << endl;
        }
    }

    void post_routine() override {
        // and now we just have to write this stuff?
        calcFile << "T";
        for(auto entry : size_T_cum_map) {
            calcFile << "," << entry.first << "," << entry.first << "_err";
        }
        calcFile << endl;
        int nr_temps = (int)size_T_cum_map.begin()->second.size();
        cout << "nr of temps: " << nr_temps << endl;
        for(int i = 0; i < nr_temps; i++) {
            calcFile << get<0>(size_T_cum_map.begin()->second[i]); // accesses the temp of the i-th tuple
            for(auto entry : size_T_cum_map) {
                calcFile << "," << get<1>(entry.second[i]) << "," << get<2>(entry.second[i]);
            }
            calcFile << endl;
        }
        calcFile.close();
    }
};

class SurBinderHandlerSilicon : public SurBinderHandler, public BinderHandlerSilicon{
public:
    SurBinderHandlerSilicon(fs::path root) : BinderHandlerSilicon(root), SurBinderHandler(root),
                                                             BinderHandler(root) {}
};

#endif //LEARNINGPROJECT_BINDERHANDLER_CPP
