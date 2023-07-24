//
// Created by andi on 22.07.23.
//

#include "LatticeFFT.h"
#include <fftw3.h>
#include "../Header/Helpfunctions and Classes.h"
#include "../../CudaProject/parameters.cuh"


using namespace std;

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
void write_to_file(const int lat_dim, const vector<vector<array<double, 2>>> &p, const vector<double> ft_squared_y,
                   const vector<double> ft_squared_x, const fs::path &writepath) {
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

int main(int argc, char* argv[]) {

    // We just read from the parameter file and compile this when we start the run, but we run when we have finished
    // the run
    const int N = 50;
    fs::path root;
    if(argc >= 2) {
        // if we give some argument, doesnt even matter what argument, we take the parameter file values
        root = "../" + adaptive_tempscan_root;
    } else {
        root = "../../../Generated content/Coulomb/system size test/Detailed-50";
    }

    vector<fs::path> temp_directories = list_dir_paths(root);
    print_vector(temp_directories);

    // need to transform
    auto q = init_q(N);
    auto p = vector<vector<array<double, 2>>>(
            N, vector<array<double, 2>>(N, array<double, 2>()));
    fill_p(q, p);

    fftw_complex *in, *out, *real_out;
    double *real_in;


    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
    real_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
    real_in = (double*) fftw_malloc(sizeof(double) * N * N);

    fftw_plan plan, p2;
    plan = fftw_plan_dft_2d(N, N, in, out, FFTW_FORWARD, FFTW_MEASURE);

    for(auto path : temp_directories) {
        cout << endl << endl << "Dealing with temp path " << path << endl;
        vector<fs::path> csv_files = list_csv_files(path);
        double ft_squared_k[N] = {};    // summed over l so still a function of k
        double ft_squared_l[N] = {};

        for(auto csv_path :csv_files) {
            // so trafo routine opens the file, fills 'in' in, does the FT of the lattice
            // adds to ft_squared_k
            trafo_routine(N, in, out, plan, ft_squared_k, ft_squared_l, csv_path);
        }
        // maybe average
        for(int i = 0; i < N; i++) {
            ft_squared_k[i] /= (double)csv_files.size() * pow(N, 4);
            ft_squared_l[i] /= (double)csv_files.size() * pow(N, 4);
        }
        // and write it
        fs::path writepath = path / "struct.fact";
        write_to_file(N, p, ft_squared_k, ft_squared_l, writepath);
    }


    //Now i just do this for every csv file?
    return 0;
}




