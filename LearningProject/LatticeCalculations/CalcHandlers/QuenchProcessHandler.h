//
// Created by andi on 29.08.23.
//

#ifndef LEARNINGPROJECT_QUENCHPROCESSHANDLER_H
#define LEARNINGPROJECT_QUENCHPROCESSHANDLER_H

#include "calcHandler.h"

//
// Created by andi on 11.08.23.
//



#include <fftw3.h>
#include "../../Header/Helpfunctions and Classes.h"
#include "calcHandler.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <unsupported/Eigen/NonLinearOptimization>
#include "../../DomainSize/QuenchCorrLength.h"

class QuenchProcessHandler : public calcHandler {
    int cell_L = 0;
public:
    QuenchProcessHandler(const fs::path& root): calcHandler(root) {
        cell_L = (int) QuenchProcessHandlerConfig["cell_L"];
    }

    void pre_routine() override {
        if(cell_L == 0) {
            cell_L = lat_dim;
        }
    }

    void setting_pre_routine(fs::path setting_path) {
        // For every setting i and every timepoint in need an array where i can add up stuff
        vector<fs::path> csv_files = list_csv_files(setting_path);
        fs::path txt_file = findFirstTxtFile(setting_path);
        double J = extractValueFromTxt(txt_file, "J");
        int nr_values = (int)extractValueFromTxt(txt_file, "nr_save_values");
        bool chessTrafo;
        if(J < 0) {
            chessTrafo = true;
        } else {
            chessTrafo = false;
        }
        fftw_complex *in, *out;
        in = new fftw_complex[cell_L*cell_L];
        out = new fftw_complex[cell_L*cell_L];

        auto q = init_q(cell_L);
        auto p = vector<vector<array<double, 2>>>(
                cell_L, vector<array<double, 2>>(cell_L, array<double, 2>()));
        fill_p(q, p);
        auto k = p_to_vec(p);

        ofstream ofile;
        ofile.open(setting_path/"quench.process");
        ofile << "t,T,xi_x,xi_y,ampl_x,ampl_y" << endl;

        auto* ft_k = new double[cell_L];
        auto* ft_l = new double[cell_L];

        for(int time_nr = 0; time_nr < nr_values; time_nr++) {

            for(int l = 0; l < cell_L; l++) {
                ft_k[l] = 0;
                ft_l[l] = 0;
            }

            double T, t;
            for(auto csv_path : csv_files) {
                // find out if chess_trafo should be true or not
                cout << "Before reading" << endl;
                ifstream file = safe_read(csv_path, true);
                cout << "After reading" << endl;
                auto lat_q = readDoubleValuesAt(file, time_nr,  T, t);
                cout << "After extracting double values at" << endl;
                if(chessTrafo) {
                    chess_trafo(lat_q);
                }
                int nr_cells = (int)(lat_dim * lat_dim/ (cell_L * cell_L));
                cout << "nr of cells = " << nr_cells << endl;
                for(int cell_nr = 0; cell_nr < nr_cells; cell_nr++) {
                    vector<double> cell = vector<double>(cell_L * cell_L, 0);
                    extract_cell(lat_q, cell, cell_nr, cell_L);
                    // TODO defining the plan here?
                    fftw_plan plan;
                    plan = fftw_plan_dft_2d(cell_L, cell_L, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
                    for (int l = 0; l < cell_L * cell_L; l++) {
                        in[l][0] = cell[l];
                        in[l][1] = 0;
                    }
                    cout << "before executing plan" << endl;
                    fftw_execute(plan);
                    cout << "after executing plan" << endl;
                    sum_and_add(cell_L, out, ft_k, ft_l);
                    fftw_destroy_plan(plan);
                }
                file.close();
            }
            // normalizing
            for(int i = 0; i < cell_L; i++) {
                ft_k[i] /= pow(cell_L, 4);
                ft_l[i] /= pow(cell_L, 4);
            }

            // do the fitting
            vector<double> ft_vec_x = vector<double>(ft_k, ft_k + cell_L);
            Eigen::VectorXd paras_X = fit_lorentz_peak(k, ft_vec_x);
            // cout << paras_X << endl;


            vector<double> ft_vec_y = vector<double>(ft_l, ft_l + cell_L);
            Eigen::VectorXd paras_y = fit_lorentz_peak(k, ft_vec_y);
            // cout << paras_y << endl;
            cout << "writing t = " << t << endl;
            ofile << t << "," << T << "," << paras_X(1) << "," << paras_y(1) << "," << paras_X(0) << "," << paras_y(0) << endl;

        }
        delete[] ft_k;
        delete[] ft_l;
        ofile.close();
    }

    void realization_routine(vector<double> &lat_q, double T, double t) override {

    }

    void setting_post_routine() {
    }


    void sum_and_add(const int N, fftw_complex const (*out), double (*ft_squared_k), double (*ft_squared_l)) {
        for(int i = 0; i < N; i++) {
            // i counts in x dimension?
            for(int j = 0; j < N; j++) {
                // if i want to average over l i need to sum over the rows, so determine row i
                int k_ind = i * N + j;
                // I sum over k the squared absolute value
                // |sigma'_kl|^2 = (sqrt(sigma'_kl.real * sigma'_kl.real + sigma'_kl.imag * sigma'_kl.imag))^2
                ft_squared_k[i] += ((out[k_ind][0] * out[k_ind][0]) + (out[k_ind][1] * out[k_ind][1]));
                // if i want to average over k, i need to sum over the columns, so determine column by fixed +i, run over
                // rows with j*N
                int l_ind = j * N + i;
                ft_squared_l[i] += ((out[l_ind][0] * out[l_ind][0]) + (out[l_ind][1] * out[l_ind][1]));

            }
        }
    }


    vector<double> p_to_vec(vector<vector<array<double, 2>>>& p) {
        vector<double> k;
        for(int j = 0; j<p.size(); j++) {
            k.push_back(p[j][j][0]);
        }
        return k;
    }

};

#endif //LEARNINGPROJECT_QUENCHPROCESSHANDLER_H
