//
// Created by andi on 24.08.23.
//

#ifndef LEARNINGPROJECT_STRUCTFACTHANDLER_H
#define LEARNINGPROJECT_STRUCTFACTHANDLER_H

#include <fftw3.h>
#include "../../Header/Helpfunctions and Classes.h"
#include "calcHandler.h"
#include "../configs.h"

class StructFactHandler : public calcHandler {
    fftw_complex *in, *out;
    fftw_plan plan;
    int n;
    size_t Lx;         // TODO it should not just stand here but whatever
    size_t Ly;
    int cutup;
    vector<vector<double>>  ft_squared_k;
    vector<vector<double>>  ft_squared_l;
    vector<double>          mean_ft_squared_k = {};
    vector<double>          mean_ft_squared_l = {};
    vector<double>          stddev_mean_k = {};
    vector<double>          stddev_mean_l = {};

    fs::path writepath;
public:

    StructFactHandler(fs::path root) : calcHandler(root) {
        // Thing is now that i have different sizes in x and y dimension to be able to capture larger corr lengths
        // in one dimension, it does not make sense to cut it up into square pieces
        // -> use a number of cells that i want to get? should be 2², 3² etc
        // -> define cells per row/col to stay true to the ratio?
        Lx = (int) StructFactConfig["cell_L"];
        cutup = (int)StructFactConfig["cutup"];

    }

    void pre_routine() override {
        // need to find out the lattice dimension
        // TODO Assuming we have all the same size?
        // ah no lets not do that
    }

    void setting_pre_routine(fs::path setting_path) override {
        if (cutup == 0) {
            Lx = dim_size_x;
            Ly = dim_size_y;
            cout << "cell Lx equals system size Lx = " << Lx << endl;
            cout << "cell Ly equals system size Ly = " << Ly << endl;
        } else {
            // so cutup is something like 2, meaning 2 cells per row
            Lx = dim_size_x / cutup;
            Ly = dim_size_y / cutup;
        }
        // we need a map that maps the system size to the vector containing all the m values
        // it has to be reset for every setting / temperature
        fs::path txt_path = findFirstTxtFile(setting_path);
        writepath = setting_path / "struct.fact";
        n = dim_size_x * dim_size_y;
        cout << "n = " << n << endl;
        cout << n << endl;
        // Allocating the arrays for the fourier trafos
        in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
        out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * n);
        // TODO is the order right here?
        plan = fftw_plan_dft_2d(Ly, Lx, in, out, FFTW_FORWARD, FFTW_MEASURE);
         // reset everything that kept track of stuff
        ft_squared_k = {};
        ft_squared_l = {};
        mean_ft_squared_k = vector<double>(Lx, 0);       // lets say k is "x-direction" ? k is
        mean_ft_squared_l = vector<double>(Ly, 0);
        stddev_mean_k = vector<double>(Lx, 0);
        stddev_mean_l = vector<double>(Ly, 0);
    }

    void realization_routine(vector<double> &lat_q, double T, double t) override {
        vector<double> ft_k_squared_of_this_csv_file = vector<double>(Lx, 0);
        vector<double> ft_l_squared_of_this_csv_file = vector<double>(Ly, 0);
        // do stuff
        // okay we want to split in cells again
        int nr_of_cells = n / (Lx * Ly);
        cout << "nr of cells = " << nr_of_cells << endl;
        vector<double> cell = vector<double>(Lx * Ly, 0);
        for(int cell_nr = 0; cell_nr < nr_of_cells; cell_nr++) {
            extract_cell(lat_q, cell, cell_nr, Lx, Ly, dim_size_x);
            // printAsMatrix(cell, Ly, Lx);
            for(int i = 0; i < Lx * Ly; i++) {
                in[i][0] = cell[i];
                in[i][1] = 0;
            }
/*            cout << endl << "fftw in array:" << endl;
            printComplexMatrix(in, Ly, Lx);*/
            fftw_execute(plan);
/*            cout << endl << "fftw out array:" << endl;
            printComplexMatrix(out, Ly, Lx);
            cout << "ft_k_squared before:" << endl;
            for(int i = 0; i < Lx; i++) {
                cout << ft_k_squared_of_this_csv_file[i] << "  ";
            }
            cout << endl;*/
            // get them values into the vectors
            for(int i = 0; i < Lx; i++) {
                // i counts in x dimension?
                for(int j = 0; j < Ly; j++) {
                    // we determine the col by i and sum over the rows by incrementing j 'average over y-direction'
                    int k_ind = j * Lx + i;
                    // I sum over k the squared absolute value
                    // |sigma'_kl|^2 = (sqrt(sigma'_kl.real * sigma'_kl.real + sigma'_kl.imag * sigma'_kl.imag))^2
                    ft_k_squared_of_this_csv_file[i] += ((out[k_ind][0] * out[k_ind][0]) + (out[k_ind][1] * out[k_ind][1]));
                }
            }
/*            cout << "ft_k_squared after:" << endl;
            for(int i = 0; i < Lx; i++) {
                cout << ft_k_squared_of_this_csv_file[i] << "  ";
            }
            cout << endl;*/
            for(int i = 0; i < Ly; i++) {
                for(int j = 0; j < Lx; j++) {
                    // I want to average over the row? so i fix the row by constant i and increment j
                    int l_ind = i * Lx + j;
                    ft_l_squared_of_this_csv_file[i] += ((out[l_ind][0] * out[l_ind][0]) + (out[l_ind][1] * out[l_ind][1]));
                }
            }
/*            cout << "ft_l_squared:" << endl;
            for(int i = 0; i < Ly; i++) {
                cout << ft_l_squared_of_this_csv_file[i] << "  ";
            }
            cout << endl;*/

            // add them to the vectors that keep track of everything
            ft_squared_k.push_back(ft_k_squared_of_this_csv_file);
            ft_squared_l.push_back(ft_l_squared_of_this_csv_file);
        }
    }

    void setting_post_routine() {
        // I need to average over the csv files and
        // thermal average
        int nr_realizaitons = ft_squared_l.size();
        for(int p = 0; p < Lx; p++) {
            for(int realization = 0; realization < nr_realizaitons; realization++) {
                mean_ft_squared_k[p] += ft_squared_k[realization][p];
                if (p < Ly) {
                    mean_ft_squared_l[p] += ft_squared_l[realization][p];
                }
            }
        }
        // averaging
        for(int p = 0; p < Lx; p++) {
            // when averaging it is actually swapped? because i think i sum over Ly entries to avg the ft_k
            mean_ft_squared_k[p] /= (double)nr_realizaitons * pow(Ly, 4);
            if (p < Ly) {
                mean_ft_squared_l[p] /= (double)nr_realizaitons * pow(Lx, 4);
            }
        }
        write_to_file();
    }

    void post_routine() override {
        cout << "Calling CorrLengthHandler post routine (UNIMPLEMENTED)" << endl;
    }

    void write_to_file() {
        auto py = get_frequencies(Ly);
        auto px = get_frequencies(Lx);
        ofstream ofile;
        ofile.open(writepath);
        ofile << "px," << "ft_avg_y,stddev_y," << "py," << "ft_avg_x,stddev_x\n";

        for(int j = 0; j<Lx; j++) {
            int Kx = Lx/2;
            int ix = (j + Kx < Lx) ? j + Kx : j - Kx;
            int ft_ind_x = ix;

            int Ky = Ly/2;
            int iy = (j + Ky < Ly) ? j + Ky : j - Ky;
            int ft_ind_y = iy;

            ofile <<  px[j] << ", " << mean_ft_squared_k[ft_ind_x] << ", " << 0;
            // TODO we just assume here that Lx is greater than Ly, but we can always do that right? only good to rmember
            if (j < Ly) {
                ofile << "," << py[j] << ", "
                        << mean_ft_squared_l[ft_ind_y] <<"," << 0;

            }
/*            cout <<  p[i][i][0] << ", " << ft_squared_y[i] << ", " << p[i][i][1] << ", " << ft_squared_x[i] << endl;
        cout << i << "  " << p[0][i][0] << "  " << ft_squared_y[i] << "  " << p[i][0][1] << "   " << ft_squared_x[i] << endl;
        */
            if(j < Lx - 1) {
                ofile << endl;
            }
        }
    }

};
#endif //LEARNINGPROJECT_STRUCTFACTHANDLER_H
