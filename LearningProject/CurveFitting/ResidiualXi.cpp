//
// Created by andi on 31.07.23.
//

#include "ResidiualXi.h"
//
// Created by andi on 20.07.23.
//

#include "../DomainSize/QuenchCorrLength.h"
#include "../DomainSize/FourierLatticeTrafo.h"
#include "../Header/Helpfunctions and Classes.h"
#include "../../CudaProject/parameters.cuh"
#include <filesystem>

void write_to_file(const int lat_dim, const vector<vector<array<double, 2>>> &p, const double *ft_squared_y,
                   const double *ft_squared_x, const fs::path &writepath) {
    ofstream ofile;
    ofile.open(writepath);
    ofile << "px, " << "ft_avg_y, " << "py, " << "ft_avg_x \n";
    for(int j = 0; j<lat_dim; j++) {
        // so px(i) is just the p_x value of every entry of p of the i-th col
        // p[0] is first row, p[0][i] is i-th entry of the first row, and p[0][i][0] is px value of the entry
        // px = p[0][i][0]
        // p[i] is the ith-row, p[i][0] is the first entry of the i-th row, which has all the same py values
        // py = p[i][0][1]
        int K = lat_dim/2;
        int i = (j + K < lat_dim) ? j + K : j - K;


        ofile <<  p[i][i][0] << ", " << ft_squared_y[i] << ", " << p[i][i][1] << ", " << ft_squared_x[i];
/*            cout <<  p[i][i][0] << ", " << ft_squared_y[i] << ", " << p[i][i][1] << ", " << ft_squared_x[i] << endl;
        cout << i << "  " << p[0][i][0] << "  " << ft_squared_y[i] << "  " << p[i][0][1] << "   " << ft_squared_x[i] << endl;
        */
        if(j < lat_dim - 1) {
            ofile << endl;
        }
    }
}

using namespace std;
namespace fs = std::filesystem;
using namespace fs;



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


void trafo_routine(const int N, fftw_complex (*in), fftw_complex const (*out), fftw_plan plan, double *ft_squared_k,
                   double *ft_squared_l, fs::path csv_path) {
    ifstream file = safe_read(csv_path);
    double T = 0;
    double t = 0;
    auto data = readDoubleValuesAt(file, -1, T, t);
    // copying the data to the in array
// fftw_complex is just double[2].


    for(int i = 0; i < N * N; i++) {
        in[i][0] = data[i];
        in[i][1] = 0;
    }
    fftw_execute(plan);

    // printComplexMatrix(out, N-1);
// sum and add
    sum_and_add(N, out, ft_squared_k, ft_squared_l);
}


vector<double> p_to_vec(vector<vector<array<double, 2>>>& p) {
    vector<double> k;
    for(int j = 0; j<p.size(); j++) {
        k.push_back(p[j][j][0]);
    }
    return k;
}

vector<int> generate_L(int starting_k, int n, int nr_Ls) {
    vector<int> L = {(int)sqrt(n) / (starting_k)};
    for (int k = starting_k + 1; k < starting_k + nr_Ls; k++) {
        int L_val = (int)sqrt(n) / (k);
        if (L.back() - L_val < 2) {
            L_val = L.back() - 2;
        }
        L.push_back(L_val);
    }
    return L;
}

int main(int argc, char* argv[]) {
    // okay we need a system that calculates the correlation func for every file or at least for every directory
    // path of the root where we have multiple directories with multiple csv for every Temp
    fs::path root;

    if(argc >= 2) {
        root = argv[1];
    } else {
        cout << "PLEASE MAKE SURE TO ADJUST LATTICE DIM" << endl;
        root = "../../../Generated content/New/Coulomb/Critical Exponent";
    }
    // lattice dim
    const int lat_dim = lattice_dim;
    const int N = 250;
    const int starting_k = 2;
    const int nr_Ls = 14;
    cout << "Lattice dim = " << N << endl;
    const vector<int> L_vec = generate_L(starting_k, N*N, nr_Ls);

    print_vector(L_vec);

    vector<fs::path> temp_directories = list_dir_paths(root);
    print_vector(temp_directories);

    auto q = init_q(N);
    auto p = vector<vector<array<double, 2>>>(
            N, vector<array<double, 2>>(N, array<double, 2>()));
    fill_p(q, p);
    auto k = p_to_vec(p);


    // We need plans and arrays for every size
    map<int, fftw_complex*> in_map = {};
    map<int, fftw_complex*> out_map = {};
    map<int, fftw_plan> plan_map = {};

    ofstream corrList;
    corrList.open(root/"corr.lengths");
    corrList << "T";

    // initializing
    for(int L : L_vec) {
        // in_map[L] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * L * L);
        // out_map[L] = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * L * L);
        fftw_complex *in, *out;
        in = new fftw_complex[L*L];
        out = new fftw_complex[L*L];
        in_map[L] = in;
        out_map[L] = out;
        plan_map[L] = fftw_plan_dft_2d(L, L, in, out, FFTW_FORWARD, FFTW_MEASURE);
        corrList << "," << L << "," << L << "_y";
        for(int i = 0; i < L * L; i++){
            out_map[L][i][0] = 0;
            out_map[L][i][1] = 0;
        }

    }


    for(auto path : temp_directories) {
        // so what do we want to do here? path is temp folder
        // i want to calculate the correlation length of a lattice of size L < n
        // by fourier transforming every sublattice, averaging at the end and then fitting
        // we could have problems for to small systems that are completely homogenous.
        // But maybe we can try using only larger sublattices.
        // After we have averaged the fourier transform of every sublattice, we do the fit
        // and get one value for the size and temp

        // for every temp i need a map that maps the system size to its running lattice fourier trafo
        map<int, double*> ft_k_map = {};
        map<int, double*> ft_l_map = {};
        // allocate memory for the maps
        for (int L : L_vec) {
            // I have L values for the fourier transform i think
            // ft_k_map[L] = new double[L];
            // ft_l_map[L] = new double[L];
            double* ft_k = new double[L];
            double* ft_l = new double[L];
            for(int l = 0; l < L; l++) {
                ft_k[l] = 0;
                ft_l[l] = 0;
            }
            ft_k_map[L] = ft_k;
            ft_l_map[L] = ft_l;
        }


        vector<fs::path> csv_files = list_csv_files(path);
        const int nr_csv_files = csv_files.size();
        // first looping over requested sizes or first looping over csv files?
        // if we looped over csv files, we would have to read the file once
        // but keep track of all the lattice ffts
        // looping over the size first we would have to read the csv files
        // for every Length
        // -> I dont know, probably both no problem but lets just read the csv file once

        // Temperature stays the same
        double T = 0, t;
        // we loop over every csv file to
        for(int i = 0; i < nr_csv_files; i++) {
            // read file
            // cout << csv_files[i] << endl;
            ifstream file = safe_read(csv_files[i], false);
            auto lat_q = readDoubleValuesAt(file, -1,  T, t);

            // fourier trafo for every size
            for(int L : L_vec) {
                // enumerate subsystems
                int nr_cells = (int)(N * N/ (L * L));
                // loop over cell and extract fourier trafo of the subsystem
                for(int cell_nr = 0; cell_nr < nr_cells; cell_nr++) {
                    // extract cell
                    vector<double> cell = vector<double>(L * L, 0);
                    extract_cell(lat_q, cell, cell_nr, L);
                    // copy data into array
                    for(int l = 0; l < L * L; l++) {
                        in_map[L][l][0] = cell[l];
                        in_map[L][l][1] = 0;
                    }

                    // fourier trafo
                    fftw_execute(plan_map[L]);
                    // add to running fourier trafo
                    sum_and_add(L, out_map[L], ft_k_map[L], ft_l_map[L]);
                }

            }
            // okay so for this file i added everything up and now have the averaged (but not normalized!)
            // fourier transformations
            // but only for one file, i also want to average over the number of csv files?
        }
        // okay here we should have add up everything and want to average now
        for(int L : L_vec) {
            for(int l = 0; l < L; l++) {
                ft_l_map[L][l] /= nr_csv_files * (pow(L, 4));
                ft_k_map[L][l] /= nr_csv_files * (pow(L, 4));
            }
        }


        // we now need to fit and write for every L
        corrList << endl << T;
        for (int L : L_vec) {
            auto q = init_q(L);
            auto p = vector<vector<array<double, 2>>>(
                    L, vector<array<double, 2>>(L, array<double, 2>()));
            fill_p(q, p);
            auto k = p_to_vec(p);

            // print_array(ft_k_map[L], L);
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
    for (int L : L_vec) {
        delete[] in_map[L];
        delete[] out_map[L];
        fftw_destroy_plan(plan_map[L]);
    }
}


