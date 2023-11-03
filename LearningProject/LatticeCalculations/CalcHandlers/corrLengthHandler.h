//
// Created by andi on 11.08.23.
//

#ifndef LEARNINGPROJECT_CORRLENGTHHANDLER_H
#define LEARNINGPROJECT_CORRLENGTHHANDLER_H

#include <fftw3.h>
#include "../../Header/Helpfunctions and Classes.h"
#include "calcHandler.h"
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Core>
#include <unsupported/Eigen/NonLinearOptimization>
#include "../../DomainSize/QuenchCorrLength.h"

class CorrLengthHandler : public calcHandler {
    // what does it do? I think this one calculates the correlation length for different cells for the L/xi plots
    int starting_k = 0;
    int nr_Ls = 10;
    vector<int> cutup_vec;
protected:
    vector<pair<int, int>> L_vec;
    ofstream corrList;
    ofstream corr2ndList;
    size_t nr_csv_files;
    map<pair<int, int>, vector<pair<double, double>>> m_map = {};
    map<pair<int, int>, double*> ft_k_map;
    map<pair<int, int>, double*> ft_l_map;
    double Temp = 0;
    int nr_of_used_cells = 0;
public:
    CorrLengthHandler(const fs::path& root): calcHandler(root) {
        corrList.open(root/"corr.lengths");
        corr2ndList.open(root/"corr2nd.lengths");
        corrList << "T";
        corr2ndList << "T";
        starting_k = (int) CorrLengthHandlerConfig["starting_k"];
        nr_Ls = (int)CorrLengthHandlerConfig["nr_Ls"];
    }

    void pre_routine() override {
        // now this one has to be completely rewritten because here I need rectangular cells
        /*      cutup_vec = vector<int>(nr_Ls);
          std::iota(std::begin(cutup_vec), std::end(cutup_vec), starting_k);
            for (int cut_factor : cutup_vec) {
                int Lx = dim_size_x / cut_factor;
                int Ly = dim_size_y / cut_factor;
                if (L_vec.size() > 0) {
                    if (L_vec.back().first - Lx < 2) {
                        // if the difference is not enough, we set it to be 2 less
                        Lx = L_vec.back().first - 2;
                    }
                    if (L_vec.back().second - Ly < 2) {
                        // if the difference is not enough, we set it to be 2 less
                        Ly = L_vec.back().second - 2;
                    }
                }
                L_vec.push_back(pair(Lx, Ly));      // TODO is this save or do I need to initialize the vector beforehand?
                corrList << "," << Lx << "," << Ly << "_y";
                corr2ndList << "," << Lx << "," << Ly << "_y";
                cout << "Lx = " << Lx << "    Ly = " << Ly << endl;
            }*/


    }

    void directory_pre_routine(path directory_path) override {
        double ratio_x_y = (double)dim_size_x / (double)dim_size_y;
        cout << "ratio_x_y = " << ratio_x_y << endl;
        int min_Lx = (int)CorrLengthHandlerConfig["min_Lx"];
        int max_Lx = (int)CorrLengthHandlerConfig["max_Lx"];
        vector<int> Lx_s = vector<int>(max_Lx - min_Lx);
        iota(begin(Lx_s), end(Lx_s), min_Lx);
        for (int Lx : Lx_s) {
            int Ly = (int)((double)Lx * ratio_x_y);
            L_vec.push_back(pair(Lx, Ly));      // TODO is this save or do I need to initialize the vector beforehand?
            corrList << "," << Lx << "," << Ly << "_y";
            corr2ndList << "," << Lx << "," << Ly << "_y";
            cout << "Lx = " << Lx << "    Ly = " << Ly << endl;
        }
    }

    void setting_pre_routine(fs::path setting_path) {
        // Here i need to refresh my maps?
        // I need maps that hold the FTs of the size for one temp
        ft_k_map = {};
        ft_l_map = {};
        vector<fs::path> csv_files = list_csv_files(setting_path);
        nr_csv_files = csv_files.size();
        for (pair<int, int> L_pair : L_vec) {
            // I have L values for the fourier transform i think
            int Lx = L_pair.first;
            int Ly = L_pair.second;
            double* ft_k = new double[Lx];
            double* ft_l = new double[Ly];
            // TODO This is probably not even necessary?
            for(int l = 0; l < Lx; l++) {
                // cout << ft_k[l] << endl;
                ft_k[l] = 0;
            }
            for(int l = 0; l < Ly; l++) {
                // cout << ft_k[l] << endl;
                ft_l[l] = 0;
            }

            // ft_k_map[L] = (double*) fftw_malloc(sizeof(double) * L * L);
            // ft_l_map[L] = (double*) fftw_malloc(sizeof(double) * L * L);
            ft_k_map[L_pair] = ft_k;
            ft_l_map[L_pair] = ft_l;
        }
        nr_of_used_cells = 0;
    }

    void realization_routine(vector<double> &lat_q, double T, double t) override {
        Temp = T;
        cout << "corrLength realization routine" << endl;


        for(pair<int, int> L_pair : L_vec) {
            // enumerate subsystems
            int Lx = L_pair.first;
            int Ly = L_pair.second;
            int nr_cells = (int)(dim_size_x * dim_size_y/ (Lx * Ly));
            fftw_complex *in, *out;
            in = new fftw_complex[Lx*Ly];
            out = new fftw_complex[Lx*Ly];
            fftw_plan plan;
            plan = fftw_plan_dft_2d(Ly, Lx, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
            // loop over cell and extract fourier trafo of the subsystem

            for(int cell_nr = 0; cell_nr < nr_cells; cell_nr++) {
                vector<double> cell = vector<double>(Lx * Ly, 0);
                extract_cell(lat_q, cell, cell_nr, Lx, Ly, dim_size_x);
                // TODO I think this is also not necessary?
                cell_routine(L_pair, in, out, plan, cell);

                // stuff for m and chi
                pair<double, double> m_L = calc_m(cell);
                m_map[L_pair].push_back(m_L);
                nr_of_used_cells += 1;
                }
            }

    }

    virtual void cell_routine(const pair<int, int> &L_pair, fftw_complex (*in), fftw_complex const (*out),
                      fftw_plan plan, const vector<double> &cell) {
        int Lx = L_pair.first;
        int Ly = L_pair.second;
        for (int l = 0; l < Lx * Ly; l++) {
            in[l][0] = cell[l];
            in[l][1] = 0;
        }
        // fourier trafo
        fftw_execute(plan);
                // out is basically    s~_kl^x and k and l are just the sum over which dimension it runs
        // for XY model we have to have one sinus transformation on cell and do the trafo with it and one cosine
        sum_and_add(Lx, Ly, out, ft_k_map[L_pair], ft_l_map[L_pair]);
    }

    void setting_post_routine() {
        average_setting();
        // we now need to fit and write for every L
        corrList << endl << Temp;
        corr2ndList << endl << Temp;
        cout << "Writing for temp = " << Temp << endl;

        for (pair<int, int> L_pair : L_vec) {
            double xix;
            double xiy;

            calc_corr_length(L_pair, xix, xiy);

            double xix_2nd;
            double xiy_2nd;

            calc_2nd_corr_length(L_pair, xix_2nd, xiy_2nd);

            corrList << "," << xix << "," << xiy;
            corr2ndList << "," << xix_2nd << "," << xiy_2nd;        // It is insanely wrong since even the fourier transform is wrong!
        }

        deallocate_ft();
    }

    void calc_2nd_corr_length(pair<int, int> L_pair, double& xix, double& xiy) {
        int Lx = L_pair.first;
        int Ly = L_pair.second;
        // index one is the correlation length
        // we now have the correlation length for one temperature for one L
        // We add it to a file that looks like
        // T    L1_x      L2      ...
        // 0.1  xix_11   xi_12   ...

        // we now want to calculate the susceptibility for a given L_pair
        auto m_vec = m_map[L_pair];
        double chi = transform_reduce(m_vec.begin(), m_vec.end(),
                                      0.0, plus<double>(),
                                      [](::pair<double, double> m) -> double { return squared(m);}) / (Lx * Ly);

        // so ft_k has Lx entries and ft_l Ly entries
        // I atm don't know where the maximum is, could be at the start...
        // lets just say it is at the second entry...
        double Fx = ft_k_map[L_pair][1];
        double Fy = ft_l_map[L_pair][1];
        // chi is the same for both i think
        xix = 1 / (2 * sin(M_PI/ Lx)) * sqrt(chi/Fx - 1);
        xiy = 1 / (2 * sin(M_PI/ Ly)) * sqrt(chi/Fy - 1);
    }

    void average_setting() {
        for(pair<int, int> L_pair : L_vec) {
            int Lx = L_pair.first;
            int Ly = L_pair.second;
            // The total number of cells, or systems, is nr_csv_files * (dim_size_x * dim_size_y) / (Lx * Ly)
            int nr_systems = nr_csv_files * (dim_size_x * dim_size_y) / (Lx * Ly);
            for(int l = 0; l < Lx; l++) {
                // TODO average by number of csv files
                // so average over all cells
                // and average over the the averaged dimension, so y in the case of ft_k
                ft_k_map[L_pair][l] /= (nr_systems * Ly);
            }
            for(int l = 0; l < Ly; l++) {
                ft_l_map[L_pair][l] /= (nr_systems * Lx);
            }
        }
    }

    void calc_corr_length(pair<int, int> L_pair, double& xix, double& xiy) {
        // We do the fit with vectors...
/*            vector<double> ft_vec_k = vector<double>(ft_k_map[L], ft_k_map[L] + sizeof(ft_k_map[L]) / (sizeof ft_k_map[L][0]));
            vector<double> ft_vec_l = vector<double>(ft_l_map[L], ft_l_map[L] + sizeof(ft_l_map[L]) / (sizeof ft_l_map[L][0]));*/
        int Lx = L_pair.first;
        int Ly = L_pair.second;
        // fitting
        // vectors for the parameters
        auto qx = init_q(Lx);
        auto px = vector<vector<array<double, 2>>>(
                Lx, vector<array<double, 2>>(Lx, array<double, 2>()));
        fill_p(qx, px);
        auto qy = init_q(Ly);
        auto py = vector<vector<array<double, 2>>>(
                Ly, vector<array<double, 2>>(Ly, array<double, 2>()));
        fill_p(qx, px);
        fill_p(qy, py);
        auto kx = p_to_vec(px);
        auto ky = p_to_vec(py);

        Eigen::VectorXd paras_x = fit_lorentz_peak(kx, ft_k_map[L_pair], Lx);
        Eigen::VectorXd paras_y = fit_lorentz_peak(ky, ft_l_map[L_pair], Ly);

        xix = paras_x(1);
        xiy = paras_y(1);
    }

    void deallocate_ft(){
        // deallocate memory of the maps?
        for (pair<int, int> L_pair : L_vec) {
            delete[] ft_k_map[L_pair];
            delete[] ft_l_map[L_pair];
        }
    }

    void post_routine() override {
        corrList.close();
        corr2ndList.close();
    }

    vector<int> generate_L(int starting_k, int n, int nr_Ls) {
        vector<int> L = {(int)sqrt(n) / (2*starting_k)};
        for (int k = starting_k + 1; k < starting_k + nr_Ls; k++) {
            int L_val = (int)sqrt(n) / (2*k);
            if (L.back() - L_val < 2) {
                L_val = L.back() - 2;
            }
            L.push_back(L_val);
        }
        return L;
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

    void sum_and_add(const int nx, const int ny, fftw_complex const (*out), double (*ft_squared_k), double (*ft_squared_l)) {
        for(int i = 0; i < nx; i++) {
            // i counts in x dimension?
            for(int j = 0; j < ny; j++) {
                // if i want to average over l i need to sum over the rows, so determine row i
                int k_ind = j * nx + i;
                // I sum over k the squared absolute value
                // |sigma'_kl|^2 = (sqrt(sigma'_kl.real * sigma'_kl.real + sigma'_kl.imag * sigma'_kl.imag))^2
                ft_squared_k[i] += ((out[k_ind][0] * out[k_ind][0]) + (out[k_ind][1] * out[k_ind][1]));
            }
        }
        // TODO it should also be possible to write this in one loop with ft_sqaured_l[j] and l_ind = i * nx + j
        for(int i = 0; i < ny; i++) {
            // i counts in x dimension?
            for(int j = 0; j < nx; j++) {
                int l_ind = i * nx + j;
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

class CorrLengthHandlerXY: public CorrLengthHandler {
    using CorrLengthHandler::CorrLengthHandler;

    void cell_routine(const pair<int, int> &L_pair, fftw_complex (*in), fftw_complex const (*out),
                              fftw_plan plan, const vector<double> &cell) override {
        int Lx = L_pair.first;
        int Ly = L_pair.second;
        for (int l = 0; l < Lx * Ly; l++) {
            // first sx
            in[l][0] = cos(cell[l]);
            in[l][1] = 0;
        }
        // fourier trafo
        fftw_execute(plan);
        // out is basically    s~_kl^x and k and l are just the sum over which dimension it runs
        // for XY model we have to have one sinus transformation on cell and do the trafo with it and one cosine
        sum_and_add(Lx, Ly, out, ft_k_map[L_pair], ft_l_map[L_pair]);
        for (int l = 0; l < Lx * Ly; l++) {
            // then sy
            in[l][0] = sin(cell[l]);
            in[l][1] = 0;
        }
        // fourier trafo
        fftw_execute(plan);
        // out is basically    s~_kl^x and k and l are just the sum over which dimension it runs
        // for XY model we have to have one sinus transformation on cell and do the trafo with it and one cosine
        sum_and_add(Lx, Ly, out, ft_k_map[L_pair], ft_l_map[L_pair]);
    }
};

class SurCorrLengthHandlerXY: public CorrLengthHandlerXY{
    using CorrLengthHandlerXY::CorrLengthHandlerXY;
    map<size_t, vector<tuple<double, double, double>>> size_T_corr_length_map{};
    size_t subsystem_Lx;
    double x_y_factor;

    void pre_routine() override {
        cout << "Calling SurCorrLenghtHandlerXY pre routine" << endl;
    }

    void directory_pre_routine(path directory_path) override {
        fs::path txt_file = findFirstTxtFile(directory_path);
        subsystem_Lx = (size_t)extractValueFromTxt(txt_file, "subsystem_Lx");
        x_y_factor = extractValueFromTxt(txt_file, "x_y_factor");
        cout << "calling directory pre routine of SurCorrlengthHandlerXY" << endl;
        cout << "subsystem_Lx = " << subsystem_Lx << endl;
        L_vec = {pair<int, int>((int)subsystem_Lx, (int) x_y_factor * subsystem_Lx)};
        // if you now go through with the setting routines and stuff, you will be left with an m_map with one size
        // and the block m's
        size_T_corr_length_map[subsystem_Lx] = vector<tuple<double, double, double>>{};
    }

    void setting_post_routine() override {
        // setting is in this case the temperature, so we write the temp to file after iterating over all realizations
        // I need the tau instead of the t here...
        // now we have the full m_map for the temperature, leaves to calculate the binder cumulant aswell as errors
        average_setting();
        cout << "T = " << Temp << "   " << "subsystem_size = " << subsystem_Lx << endl;
        cout << "total cells used for this setting: " << nr_of_used_cells << endl;
        // We want to fit xi now for the only L_pair we got (that is as large as the system)
        double xix, xiy;
        pair<int, int> L_pair = L_vec[0];       // L_vec should have only one entry
        calc_corr_length(L_pair, xix, xiy);
        size_T_corr_length_map[subsystem_Lx].push_back(tuple<double, double, double>(Temp, xix, xiy));
        for(auto tuples : size_T_corr_length_map[subsystem_Lx]) {
            cout << "(" << get<0>(tuples) << ", " << get<1>(tuples) << ", " << get<2>(tuples) << ")" << endl;
        }
        deallocate_ft();
    }


    void post_routine() override {
        // and now we just have to write this stuff?
        cout << "SurCorrLength post routine called" << endl;
        for(auto entry : size_T_corr_length_map) {
            cout << entry.first << endl;
            corrList << "," << entry.first << "," << (int)((double)entry.first * x_y_factor) << "_y";
        }
        corrList << endl;
        int nr_temps = (int)size_T_corr_length_map.begin()->second.size();
        cout << "nr of temps: " << nr_temps << endl;
        for(int i = 0; i < nr_temps; i++) {
            corrList << get<0>(size_T_corr_length_map.begin()->second[i]); // accesses the temp of the i-th tuple
            for(auto entry : size_T_corr_length_map) {
                corrList << "," << get<1>(entry.second[i]) << "," << get<2>(entry.second[i]);
            }
            corrList << endl;
        }
        corrList.close();
    }

};
#endif //LEARNINGPROJECT_CORRLENGTHHANDLER_H
