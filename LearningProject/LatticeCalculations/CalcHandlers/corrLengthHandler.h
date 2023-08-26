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
    int starting_k = 0;
    int nr_Ls = 10;
    vector<int> L_vec;
    ofstream corrList;
    map<int, double*> ft_k_map;
    map<int, double*> ft_l_map;
    double Temp = 0;
public:
    CorrLengthHandler(const fs::path& root): calcHandler(root) {
        corrList.open(root/"corr.lengths");
        corrList << "T";
        starting_k = (int) CorrLengthHandlerConfig["starting_k"];
        nr_Ls = (int)CorrLengthHandlerConfig["nr_Ls"];
    }

    void pre_routine() override {
        L_vec = generate_L(starting_k, lat_dim * lat_dim, nr_Ls);
        for (int L : L_vec) {
            corrList << "," << L << "," << L << "_y";
        }
    }

    void setting_pre_routine(fs::path setting_path) {
        // Here i need to refresh my maps?
        // I need maps that hold the FTs of the size for one temp
        ft_k_map = {};
        ft_l_map = {};

    }

    void realization_routine(vector<double> &lat_q, double T, double t) override {
        Temp = T;
        for (int L : L_vec) {
            // I have L values for the fourier transform i think
            double* ft_k = new double[L];
            double* ft_l = new double[L];
            for(int l = 0; l < L; l++) {
                ft_k[l] = 0;
                ft_l[l] = 0;
            }

            // ft_k_map[L] = (double*) fftw_malloc(sizeof(double) * L * L);
            // ft_l_map[L] = (double*) fftw_malloc(sizeof(double) * L * L);
            ft_k_map[L] = ft_k;
            ft_l_map[L] = ft_l;
        }

        for(int L : L_vec) {
            // enumerate subsystems
            int nr_cells = (int)(lat_dim * lat_dim/ (L * L));
            fftw_complex *in, *out;
            in = new fftw_complex[L*L];
            out = new fftw_complex[L*L];
            fftw_plan plan;
            plan = fftw_plan_dft_2d(L, L, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
            // loop over cell and extract fourier trafo of the subsystem

            for(int cell_nr = 0; cell_nr < nr_cells; cell_nr++) {
                vector<double> cell = vector<double>(L * L, 0);
                extract_cell(lat_q, cell, cell_nr, L);
                for (int l = 0; l < L * L; l++) {
                    in[l][0] = cell[l];
                    in[l][1] = 0;
                }
                // fourier trafo
                fftw_execute(plan);
                sum_and_add(L, out, ft_k_map[L], ft_l_map[L]);
            }
            }

    }

    void setting_post_routine() {
        for(int L : L_vec) {
            for(int l = 0; l < L; l++) {
                // TODO average by number of csv files
                ft_l_map[L][l] /= (pow(L, 4));
                ft_k_map[L][l] /= (pow(L, 4));
            }
        }
        // we now need to fit and write for every L
        corrList << endl << Temp;
        cout << "Writing for temp = " << Temp << endl;

        for (int L : L_vec) {
            // We do the fit with vectors...
/*            vector<double> ft_vec_k = vector<double>(ft_k_map[L], ft_k_map[L] + sizeof(ft_k_map[L]) / (sizeof ft_k_map[L][0]));
            vector<double> ft_vec_l = vector<double>(ft_l_map[L], ft_l_map[L] + sizeof(ft_l_map[L]) / (sizeof ft_l_map[L][0]));*/

            // fitting
            // vectors for the parameters

            auto q = init_q(L);
            auto p = vector<vector<array<double, 2>>>(
                    L, vector<array<double, 2>>(L, array<double, 2>()));
            fill_p(q, p);
            auto k = p_to_vec(p);

            Eigen::VectorXd paras_x = fit_lorentz_peak(k, ft_k_map[L], L);
            Eigen::VectorXd paras_y = fit_lorentz_peak(k, ft_l_map[L], L);
            // index one is the correlation length
            // we now have the correlation length for one temperature for one L
            // We add it to a file that looks like
            // T    L1_x      L2      ...
            // 0.1  xix_11   xi_12   ...
            corrList << "," << paras_x(1) << "," << paras_y(1);
        }

        // deallocate memory of the maps?
        for (int L : L_vec) {
            delete[] ft_k_map[L];
            delete[] ft_l_map[L];
        }
    }

    void post_routine() override {
        corrList.close();
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


    vector<double> p_to_vec(vector<vector<array<double, 2>>>& p) {
        vector<double> k;
        for(int j = 0; j<p.size(); j++) {
            k.push_back(p[j][j][0]);
        }
        return k;
    }

};

#endif //LEARNINGPROJECT_CORRLENGTHHANDLER_H
