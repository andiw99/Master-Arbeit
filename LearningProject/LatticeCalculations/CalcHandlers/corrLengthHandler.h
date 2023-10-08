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
    vector<pair<int, int>> L_vec;
    vector<int> cutup_vec;
    ofstream corrList;
    map<pair<int, int>, double*> ft_k_map;
    map<pair<int, int>, double*> ft_l_map;
    double Temp = 0;
public:
    CorrLengthHandler(const fs::path& root): calcHandler(root) {
        corrList.open(root/"corr.lengths");
        corrList << "T";
        starting_k = (int) CorrLengthHandlerConfig["starting_k"];
        nr_Ls = (int)CorrLengthHandlerConfig["nr_Ls"];
    }

    void pre_routine() override {
        // now this one has to be completely rewritten because here I need rectangular cells
        cutup_vec = vector<int>(nr_Ls);
        std::iota(std::begin(cutup_vec), std::end(cutup_vec), starting_k);
        for (int cut_factor : cutup_vec) {
            int Lx = dim_size_x / cut_factor;
            int Ly = dim_size_y / cut_factor;
            L_vec.push_back(pair(Lx, Ly));      // TODO is this save or do I need to initialize the vector beforehand?
            corrList << "," << Lx << "," << Ly << "_y";
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
        cout << "corrLength realization routine" << endl;
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
                for (int l = 0; l < Lx * Ly; l++) {
                    in[l][0] = cell[l];
                    in[l][1] = 0;
                }
                // fourier trafo
                fftw_execute(plan);
                sum_and_add(Lx, Ly, out, ft_k_map[L_pair], ft_l_map[L_pair]);
            }
            }

    }

    void setting_post_routine() {
        for(pair<int, int> L_pair : L_vec) {
            int Lx = L_pair.first;
            int Ly = L_pair.second;
            for(int l = 0; l < Lx; l++) {
                // TODO average by number of csv files
                ft_k_map[L_pair][l] /= (pow(Ly, 4));
            }
            for(int l = 0; l < Ly; l++) {
                ft_l_map[L_pair][l] /= (pow(Lx, 4));
            }
        }
        // we now need to fit and write for every L
        corrList << endl << Temp;
        cout << "Writing for temp = " << Temp << endl;

        for (pair<int, int> L_pair : L_vec) {
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

            print_vector(kx);
            cout << endl;
            print_array(ft_k_map[L_pair], Lx);

            print_vector(ky);
            cout << endl;
            print_array(ft_l_map[L_pair], Ly);


            Eigen::VectorXd paras_x = fit_lorentz_peak(kx, ft_k_map[L_pair], Lx);
            Eigen::VectorXd paras_y = fit_lorentz_peak(ky, ft_l_map[L_pair], Ly);
            // index one is the correlation length
            // we now have the correlation length for one temperature for one L
            // We add it to a file that looks like
            // T    L1_x      L2      ...
            // 0.1  xix_11   xi_12   ...
            cout << "xi_x = " << paras_x(1) << endl;
            cout << "xi_y = " << paras_y(1) << endl;
            corrList << "," << paras_x(1) << "," << paras_y(1);
        }

        // deallocate memory of the maps?
        for (pair<int, int> L_pair : L_vec) {
            delete[] ft_k_map[L_pair];
            delete[] ft_l_map[L_pair];
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

#endif //LEARNINGPROJECT_CORRLENGTHHANDLER_H
