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
    int starting_k = 10;                 // should find way to initialize them in create? Probably not
    int nr_Ls = 20;                     // to comfortable to always have to change them in the class
    vector<int> L_vec;
    map<int, vector<double>> m_map = {};
    double T;
public:

    BinderHandler(fs::path root): calcHandler(root) {
        starting_k = (int) BinderHandlerConfig["starting_k"];
        nr_Ls = (int) BinderHandlerConfig["nr_Ls"];
    }

    void pre_routine() override {
        cout << "Calling BinderHandler pre routine" << endl;
        calcFile.open(root / "binder.cumulants");
        // Generating the vector of Ls
        // Okay i think we have to rewrite this entirely since we want to cut up into rectangles (or do we want to?)
        // TODO think about whether you want to calculate the binder cumulant in rectangular cells, for now i dont think so
        L_vec = generate_L(starting_k, dim_size_x * dim_size_y, nr_Ls);

        init_file_header(calcFile, L_vec);

    }

    void setting_pre_routine(fs::path setting_path) override {
        // we need a map that maps the system size to the vector containing all the m values
        // it has to be reset for every setting / temperature
        m_map = {};
    }

    void realization_routine(vector<double> &lat_q, double temp, double t) override {
        // for every subsystem size L we need to iterate over the lattice and extract all m's
        for(int L : L_vec) {
            int nr_cells = (int)(dim_size_x * dim_size_y / (L * L));
            for(int cell_nr = 0; cell_nr < nr_cells; cell_nr++) {
                double m_L;
                // now we are in the ith cell and can start calculating the local magnetization powers
                // it might be easiest and maybe also the most efficient to extract the spins of the cell
                // out of the large lattice
                vector<double> q_cell = vector<double>(L * L, 0);
                extract_cell(lat_q, q_cell, cell_nr, L, L, dim_size_x);
                m_L = mean(q_cell);
                m_map[L].push_back(m_L);
            }
        }
        T = temp;
    }

    void setting_post_routine() override {
        // setting is in this case the temperature, so we write the temp to file after iterating over all realizations
        // I need the tau instead of the t here...
        calcFile << T;
        // now we have the full m_map for the temperature, leaves to calculate the binder cumulant aswell as errors
        calc_write_cum(m_map, calcFile);
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

    void calc_write_cum (map<int, vector<double>> m_map, ofstream& cumsList) {
        for(auto pair : boost::adaptors::reverse(m_map)){
            double m_L2 = std::transform_reduce(pair.second.begin(), pair.second.end(),
                                                0.0, // initial value for the reduction (sum)
                                                std::plus<double>(), // transformation (square)
                                                [](double b) { return b * b; });
            m_L2 /= pair.second.size();
            double m_L2_err = std::transform_reduce(pair.second.begin(), pair.second.end(),
                                                    0.0, // initial value for the reduction (sum)
                                                    std::plus<double>(), // transformation (square)
                                                    [&m_L2](double b) { return pow((b * b - m_L2), 2)  ; });
            m_L2_err /= (double)pow(pair.second.size(), 2);
            m_L2_err = sqrt(m_L2_err);
            double m_L4 = std::transform_reduce(pair.second.begin(), pair.second.end(),
                                                0.0, // initial value for the reduction (sum)
                                                std::plus<>(), // transformation (square)
                                                [](double b) { return (pow(b, 4)); });
            m_L4 /= pair.second.size();
            double m_L4_err = std::transform_reduce(pair.second.begin(), pair.second.end(),
                                                    0.0, // initial value for the reduction (sum)
                                                    std::plus<>(), // transformation (square)
                                                    [&m_L4](double b) { return pow(pow(b, 4) - m_L4, 2); });
            m_L4_err /= (double)pow(pair.second.size(), 2);
            m_L4_err = sqrt(m_L4_err);
            double cum = m_L4 / (m_L2 * m_L2);
            double cum_error = sqrt(pow(1 / m_L2 / m_L2 * m_L4_err, 2) + pow(2 * m_L4 / pow(m_L2, 3) * m_L2_err, 2));
            cumsList << "," << cum << "," << cum_error;
        }
    }
};

#endif //LEARNINGPROJECT_BINDERHANDLER_CPP
