//
// Created by andi on 20.07.23.
//

#include "QuenchCorrLength.h"
#include "FourierLatticeTrafo.h"
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


template <size_t N>
void average_and_add(fftw_complex ft[N][N], double (&ft_squared_y)[N], double (&ft_squared_x)[N]) {
    // supposed to calculate ft_avg_y for one lattice and add it to ft_avg_y

    array<array<double, 2>, N> ft_y = {};
    array<array<double, 2>, N> ft_x = {};

    array<double, N> ft_avg_x = {};
    array<double, N> ft_avg_y = {};

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            // cout << ft[j][i][0] << ", " << ft[i][j][1];
/*            if(i == 0 & j == 0) {
                cout << *ft[i][j] <<  ", " << ft[i][j][0]<<  ", " << ft[i][j][1] << endl;
            }*/
            ft_y[i][0] += ft[j][i][0];
            ft_y[i][1] += ft[j][i][1];
            ft_x[i][0] += ft[i][j][0];
            ft_x[i][1] += ft[i][j][1];
        }
        // actually N^8 i dont know really why but i think i probably knew what I was doing
/*        if(i == 1) {
            cout << "ft_y:" << endl;
            cout << ft_y[i][0] << ", " << ft_y[i][1]  << endl;
        }*/
        ft_avg_y[i] = (ft_y[i][0] * ft_y[i][0]) / pow((double) N, 4); // +  (ft_y[i][1] * ft_y[i][1]) / pow((double) N, 4);
        ft_avg_x[i] = (ft_x[i][0] * ft_x[i][0]) / pow((double) N, 4); // +  (ft_x[i][1] * ft_x[i][1]) / pow((double) N, 4);
    }
    // now another loop that adds the now calculated squared absolute values to the running squared abs
    for(int i = 0; i < N; i++) {
        ft_squared_y[i] += ft_avg_y[i];
        ft_squared_x[i] += ft_avg_x[i];
    }
    // that should be it
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

void fft_trafo_routine(const int N, int ind, fftw_complex (*in), fftw_complex const (*out), fftw_plan plan, double *ft_squared_k,
                       double *ft_squared_l, fs::path csv_path, double &T, double& t, bool chessTrafo=true) {
    // needs to take in the references to the ft_avg_y and ft_avg_x and add them up for several files
    // and now it needs to add up everytime for every file, so it needs to take in the filepath
/*    cout << "reading: " << filepath << endl;
    ifstream file(filepath);
    // check if file is opened
    if(file.is_open()) {
        // cout << "File successfully opened" << endl;
    } else {
        cout << "Failed to open file" << endl;
        // abort if file is not opened
        exit(0);
    }
    // and now we do the stuff we did all the time?
    // read last line and put it into a 2D vector
    vector<complex<double>> values = readValuesAt(file, ind, T, t);
    vector<vector<complex<double>>> f = oneD_to_twoD(values);


    // init the empty arrays for the fftw
    fftw_complex f_fftw[N][N] = {};
    fftw_complex ft_fftw[N][N] = {};
    // copy the values from f to f_fftw
    copy_values2D<N>(f, f_fftw);
    // print2DContainer(f_fftw);
    // do the trafo
    fftw_plan plan2d;
    plan2d = fftw_plan_dft_2d(N, N, &f_fftw[0][0], &ft_fftw[0][0], FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan2d);

    // now we got to calculate the avg abs squared for every row and col and add them to ft_avg_y and ft_avg_x
    average_and_add(ft_fftw, ft_squared_y, ft_squared_x);*/
    ifstream file = safe_read(csv_path, false);
    auto data = readDoubleValuesAt(file, ind, T, t);
    if(chessTrafo) {
        chess_trafo(data);
    }
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



/*
Eigen::VectorXd fit_lorentz_peak(vector<double>& k_values, vector<double>& ft_values) {
    // constrtuct the matrix that holds the k and ft values
    Eigen::MatrixXd X_Y_vals(k_values.size(), 2);
    X_Y_vals = construct_matrix(k_values, ft_values);
    printMatrixXd(X_Y_vals);

    Eigen::VectorXd params(2);
    // the params of the fit are the scaling and the correlation length? params(0) = amplitude params(1) = xi
    params << 1.0, 1.0;
    LorentzianPeakFunctor functor(X_Y_vals);
    Eigen::NumericalDiff<LorentzianPeakFunctor> numericalDiff(functor);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<LorentzianPeakFunctor>> lm(numericalDiff);

    Eigen::LevenbergMarquardtSpace::Status status = lm.minimize(params);
    std::cout << "status: " << status << std::endl;
    return params;
}
*/

vector<double> p_to_vec(vector<vector<array<double, 2>>>& p) {
    vector<double> k;
    for(int j = 0; j<p.size(); j++) {
        k.push_back(p[j][j][0]);
    }
    return k;
}

int main(int argc, char* argv[]) {
    // okay we need a system that calculates the correlation func for every file or at least for every directory
    // path of the root where we have multiple directories with multiple csv for every Temp
    fs::path root;

    if(argc >= 2) {
        root = argv[1];
    } else {
        cout << "PLEASE MAKE SURE TO ADJUST LATTICE DIM AND TURN CHESS TRAFO ON/OFF" << endl;
        root = "../../../Generated content/AA/AA Quench";
    }
    // lattice dim
    const int lat_dim = lattice_dim;
    const int N = 160;
    bool chessTrafo = true;
    cout << "Lattice dim = " << N << endl;

    vector<fs::path> temp_directories = list_dir_paths(root);
    print_vector(temp_directories);

    auto q = init_q(N);
    auto p = vector<vector<array<double, 2>>>(
            N, vector<array<double, 2>>(N, array<double, 2>()));
    fill_p(q, p);
    auto k = p_to_vec(p);

    fftw_complex *in, *out, *real_out;
    double *real_in;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
    real_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
    real_in = (double*) fftw_malloc(sizeof(double) * N * N);

    fftw_plan plan, p2;
    plan = fftw_plan_dft_2d(N, N, in, out, FFTW_FORWARD, FFTW_MEASURE);

    for(auto path : temp_directories) {
        // so path corresponds to a directory containing multiple realizations for one tau.
        // We now want to calculate for every Temperature of those Quenches the correlation length
        // To do that, we need to fourier transform every realization at a specific temperature and add up the averaged
        // fourier transform
        // then we average the fourier transforms
        // then we fit and extract xi_x and xi_y
        // save xi for every temperature (and time)

        // cout << endl << endl << "Dealing with tau path " << path << endl;
        // those are the different realizations of the Quenches.
        vector<fs::path> csv_files = list_csv_files(path);

        // we need a vector / array / matrix that contains the time, temperature and calculated correlation length
        vector<double> corr_lengths_x = {};
        vector<double> corr_lengths_y = {};
        vector<double> times = {};
        vector<double> temps = {};

        // we have to loop over the different times now. We should find out how many rows our csv files have
        // do we need that or do we just keep going until we are done?
        int nr_rows = getNrRows(csv_files[0]);
        fs::path writepath = path / "quench.process";
        ofstream file;
        file.open(writepath);
        // Create Header
        file << "t,T,xi_x,xi_y,ampl_x,ampl_y" << endl;

        for(int i = 0; i < nr_rows; i++) {
            // Now we do everything we always did but for every line
            // we need running arrays for the averages over the lattices
            // those have to be reset for every temperature
            double ft_squared_y[N] = {};
            double ft_squared_x[N] = {};
            double T;
            double t;

            for(auto csv_path :csv_files) {
                // here i am going through line i on every realization and add up ft_squared_y
                fft_trafo_routine(N, i, in, out, plan, ft_squared_y, ft_squared_x, csv_path, T, t, chessTrafo);
            }

            // and now just average over the run size
            for(int i = 0; i < N; i++) {
                ft_squared_x[i] /= (double)csv_files.size();
                ft_squared_y[i] /= (double)csv_files.size();
            }

            // Okay now we need to do the fit as we have the lorentz peak.
            // the dimension of ft_squared should be the same as the dimonsion of the p vector that was constructed one
            // time at the beginning
            vector<double> ft_vec_x = vector<double>(ft_squared_x, ft_squared_x + sizeof(ft_squared_x) / (sizeof ft_squared_x[0]));
            Eigen::VectorXd paras_X = fit_lorentz_peak(k, ft_vec_x);
            // cout << paras_X << endl;

            corr_lengths_x.push_back(paras_X(1));

            vector<double> ft_vec_y = vector<double>(ft_squared_y, ft_squared_y + sizeof(ft_squared_y) / (sizeof ft_squared_y[0]));
            Eigen::VectorXd paras_y = fit_lorentz_peak(k, ft_vec_y);
            // cout << paras_y << endl;

            corr_lengths_y.push_back(paras_y(1));
            times.push_back(t);
            temps.push_back(T);
            file << t << "," << T << "," << paras_X(1) << "," << paras_y(1) << "," << paras_X(0) << "," << paras_y(0) << endl;
            // Find out temp and time...
        }
    }
}


