//
// Created by andi on 24.08.23.
//

#ifndef LEARNINGPROJECT_STRUCTFACTHANDLER_H
#define LEARNINGPROJECT_STRUCTFACTHANDLER_H

#include <fftw3.h>
#include "../../Header/Helpfunctions and Classes.h"
#include "calcHandler.h"

class StructFactHandler : public calcHandler {
    using calcHandler::calcHandler;
    fftw_complex *in, *out;
    fftw_plan plan;
    int n;
    int L = 64;         // TODO it should not just stand here but whatever

    vector<vector<double>>  ft_squared_k;
    vector<vector<double>>  ft_squared_l;
    vector<double>          mean_ft_squared_k = {};
    vector<double>          mean_ft_squared_l = {};
    vector<double>          stddev_mean_k = {};
    vector<double>          stddev_mean_l = {};

    fs::path writepath;
public:

    void pre_routine() override {
        // need to find out the lattice dimension
        // TODO Assuming we have all the same size?
        // ah no lets not do that
    }

    void setting_pre_routine(fs::path setting_path) override {
        // we need a map that maps the system size to the vector containing all the m values
        // it has to be reset for every setting / temperature
        fs::path txt_path = findFirstTxtFile(setting_path);
        writepath = setting_path / "struct.fact";
        n = lat_dim * lat_dim;
        cout << "n = " << n << endl;
        cout << n << endl;
        // Allocating the arrays for the fourier trafos
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
        plan = fftw_plan_dft_2d(L, L, in, out, FFTW_FORWARD, FFTW_MEASURE);
         // reset everything that kept track of stuff
        ft_squared_k = {};
        ft_squared_l = {};
        mean_ft_squared_k = vector<double>(L, 0);
        mean_ft_squared_l = vector<double>(L, 0);
        stddev_mean_k = vector<double>(L, 0);
        stddev_mean_l = vector<double>(L, 0);
    }

    void realization_routine(vector<double> &lat_q, double T, double t) override {
        vector<double> ft_k_squared_of_this_csv_file = vector<double>(L, 0);
        vector<double> ft_l_squared_of_this_csv_file = vector<double>(L, 0);
        // do stuff
        // okay we want to split in cells again
        int nr_of_cells = n / (L * L);
        cout << "nr of cells = " << nr_of_cells << endl;
        vector<double> cell = vector<double>(L * L, 0);
        for(int cell_nr = 0; cell_nr < nr_of_cells; cell_nr++) {
            extract_cell(lat_q, cell, cell_nr, L);
            for(int i = 0; i < n; i++) {
                in[i][0] = cell[i];
                in[i][1] = 0;
            }
            fftw_execute(plan);
            // get them values into the vectors
            for(int i = 0; i < L; i++) {
                // i counts in x dimension?
                for(int j = 0; j < L; j++) {
                    // if i want to average over l i need to sum over the rows, so determine row i
                    int k_ind = i * L + j;
                    // I sum over k the squared absolute value
                    // |sigma'_kl|^2 = (sqrt(sigma'_kl.real * sigma'_kl.real + sigma'_kl.imag * sigma'_kl.imag))^2
                    ft_k_squared_of_this_csv_file[i] += ((out[k_ind][0] * out[k_ind][0]) + (out[k_ind][1] * out[k_ind][1]));
                    // if i want to average over k, i need to sum over the columns, so determine column by fixed +i, run over
                    // rows with j*N
                    int l_ind = j * L + i;
                    ft_l_squared_of_this_csv_file[i] += ((out[l_ind][0] * out[l_ind][0]) + (out[l_ind][1] * out[l_ind][1]));

                }
            }
            // add them to the vectors that keep track of everything
            ft_squared_k.push_back(ft_k_squared_of_this_csv_file);
            ft_squared_l.push_back(ft_l_squared_of_this_csv_file);
        }
    }

    void setting_post_routine() {
        // I need to average over the csv files and
        // thermal average
        int nr_realizaitons = ft_squared_l.size();
        for(int p = 0; p < L; p++) {
            for(int realization = 0; realization < nr_realizaitons; realization++) {
                mean_ft_squared_k[p] += ft_squared_k[realization][p];
                mean_ft_squared_l[p] += ft_squared_l[realization][p];
            }
        }
        // averaging
        for(int p = 0; p < L; p++) {
            mean_ft_squared_k[p] /= (double)nr_realizaitons * pow(L, 4);
            mean_ft_squared_l[p] /= (double)nr_realizaitons * pow(L, 4);
        }
        write_to_file();
    }

    void post_routine() override {
        cout << "Calling CorrLengthHandler post routine (UNIMPLEMENTED)" << endl;
    }

    void write_to_file() {
        auto q = init_q(L);
        auto p = vector<vector<array<double, 2>>>(
                L, vector<array<double, 2>>(L, array<double, 2>()));
        fill_p(q, p);
        ofstream ofile;
        ofile.open(writepath);
        ofile << "px," << "ft_avg_y,stddev_y," << "py," << "ft_avg_x,stddev_x\n";
        for(int j = 0; j<L; j++) {
            int K = L/2;
            int i = (j + K < L) ? j + K : j - K;


            ofile <<  p[i][i][0] << ", " << mean_ft_squared_k[j] << ", " << 0 << "," << p[i][i][1] << ", "
                  << mean_ft_squared_l[j] <<"," << 0;
/*            cout <<  p[i][i][0] << ", " << ft_squared_y[i] << ", " << p[i][i][1] << ", " << ft_squared_x[i] << endl;
        cout << i << "  " << p[0][i][0] << "  " << ft_squared_y[i] << "  " << p[i][0][1] << "   " << ft_squared_x[i] << endl;
        */
            if(j < L - 1) {
                ofile << endl;
            }
        }
    }

    vector<vector<array<double, 2>>> init_q(size_t lat_dim) {
        // later we can either overload or adjust this function to work with lattice spacings
        vector<vector<array<double, 2>>> q = vector<vector<array<double, 2>>>(
                lat_dim, vector<array<double, 2>>(lat_dim, array<double, 2>()));
        // iterate over the lattice
        for(int m = 0; m < lat_dim; m++) {
            for(int n = 0; n < lat_dim; n++) {
                q[m][n][0]= 1.0 * n;
                q[m][n][1] = 1.0 * m;
            }
        }
        return q;
    }
    void fill_p(const vector<vector<array<double, 2>>> &q, vector<vector<array<double, 2>>> &p) {
        int lat_dim = q.size();         // f has size of lattice dim
        int N = lat_dim;
        // To center the impulses around the center, we shift by K = N/2
        int K = N / 2;
        // qi corresponds to xn, so i is the col
        // find out the lattice spacings a_x = x_1 - x_0
        double ax = q[0][1][0] - q[0][0][0];
        double ay = q[1][0][1] - q[0][0][1];
        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                int i_ft = i < K ? i : i - N;
                int j_ft = j < K ? j : j - N;
                double p_i = 2 * M_PI * (double)i_ft / N / ax;
                double p_j = 2 * M_PI * (double)j_ft / N / ay;
                p[i][j] = array<double, 2>{p_i, p_j};
            }
        }
    }
};
#endif //LEARNINGPROJECT_STRUCTFACTHANDLER_H
