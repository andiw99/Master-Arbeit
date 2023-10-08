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
    int cutup = 1;
public:
    QuenchProcessHandler(const fs::path& root): calcHandler(root) {
        cutup = (int) QuenchProcessHandlerConfig["cutup"];
    }

    void pre_routine() override {

    }

    void setting_pre_routine(fs::path setting_path) {
        // get the cell_Ls
        int Lx = dim_size_x / cutup;
        int Ly = dim_size_y / cutup;

        // For every setting i and every timepoint in need an array where i can add up stuff
        vector<fs::path> csv_files = list_csv_files(setting_path);
        fs::path txt_file = findFirstTxtFile(setting_path);
        double J = extractValueFromTxt(txt_file, "J");
        int nr_values = (int)extractValueFromTxt(txt_file, "nr_saves");
        bool chessTrafo;
        if(J < 0) {
            chessTrafo = true;
        } else {
            chessTrafo = false;
        }
        fftw_complex *in, *out;
        in = new fftw_complex[Lx*Ly];
        out = new fftw_complex[Lx*Ly];

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

        ofstream ofile;
        ofile.open(setting_path/"quench.process");
        ofile << "t,T,xi_x,xi_y,ampl_x,ampl_y" << endl;

        auto* ft_k = new double[Lx];
        auto* ft_l = new double[Ly];

        cout << "nr_values = " << nr_values << endl;

        for(int time_nr = 0; time_nr < nr_values; time_nr++) {

            // TODO probably not necessary
            for(int l = 0; l < Lx; l++) {
                ft_k[l] = 0;
            }

            for(int l = 0; l < Ly; l++) {
                ft_l[l] = 0;
            }

            double T, t;
            for(auto csv_path : csv_files) {
                // find out if chess_trafo should be true or not
                ifstream file = safe_read(csv_path, true);
                auto lat_q = readDoubleValuesAt(file, time_nr,  T, t);
                if(chessTrafo) {
                    chess_trafo(lat_q);
                }
                int nr_cells = (int)(dim_size_x * dim_size_y/ (Lx * Ly));       // actually to complicated, its cust cutup²
                cout << "nr of cells = " << nr_cells << endl;
                for(int cell_nr = 0; cell_nr < nr_cells; cell_nr++) {
                    vector<double> cell = vector<double>(Lx * Ly, 0);
                    extract_cell(lat_q, cell, cell_nr, Lx, Ly, dim_size_x);
                    // TODO defining the plan here?
                    fftw_plan plan;
                    plan = fftw_plan_dft_2d(Ly, Lx, in, out, FFTW_FORWARD, FFTW_ESTIMATE);
                    for (int l = 0; l < Lx * Ly; l++) {
                        in[l][0] = cell[l];
                        in[l][1] = 0;
                    }
                    fftw_execute(plan);
                    sum_and_add(Lx, Ly, out, ft_k, ft_l);
                    fftw_destroy_plan(plan);
                }
                file.close();
            }
            // normalizing
            for(int i = 0; i < Lx; i++) {
                ft_k[i] /= pow(Ly, 4);
            }
            for(int i = 0; i < Ly; i++) {
                ft_l[i] /= pow(Lx, 4);
            }

            // do the fitting
            vector<double> ft_vec_x = vector<double>(ft_k, ft_k + Lx);
            Eigen::VectorXd paras_X = fit_lorentz_peak(kx, ft_vec_x);
            // cout << paras_X << endl;


            vector<double> ft_vec_y = vector<double>(ft_l, ft_l + Ly);
            Eigen::VectorXd paras_y = fit_lorentz_peak(ky, ft_vec_y);
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

#endif //LEARNINGPROJECT_QUENCHPROCESSHANDLER_H
