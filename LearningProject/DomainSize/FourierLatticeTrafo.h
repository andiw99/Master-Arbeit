//
// Created by andi on 28.04.23.
//

#ifndef LEARNINGPROJECT_FOURIERLATTICETRAFO_H
#define LEARNINGPROJECT_FOURIERLATTICETRAFO_H
#include <fstream>
#include <string>
#include <vector>
#include <iostream>
#include <chrono>
#include <array>
#include <complex>
#include <fftw3.h>
#include <filesystem>
#include "../Header/Helpfunctions and Classes.h"

using namespace std;
namespace fs = std::filesystem;


vector<fs::path> list_dir_paths(const fs::path& root)
{
    vector<fs::path> dir_paths;

    if (fs::is_directory(root))
    {
        for (const auto& entry : fs::directory_iterator(root))
        {
            // if regular file, do nothing?
            // if directory, add
            if (fs::is_directory(entry.path()) && entry.path().filename() != "plots")
            {
                dir_paths.push_back(entry.path());
                vector<fs::path> sub_dir_paths = list_dir_paths(entry.path());
                dir_paths.insert(dir_paths.end(), sub_dir_paths.begin(), sub_dir_paths.end());
            }
        }
    }
    return dir_paths;
}


vector<fs::path> list_csv_files(const fs::path& root)
{
    vector<fs::path> csv_files;

    if (fs::is_directory(root))
    {
        for (const auto& entry : fs::directory_iterator(root))
        {
            if (fs::is_regular_file(entry.path()) && entry.path().extension() == ".csv")
            {
                csv_files.push_back(entry.path());
            }
            else if (fs::is_directory(entry.path()))
            {
                vector<fs::path> sub_csv_files = list_csv_files(entry.path());
                csv_files.insert(csv_files.end(), sub_csv_files.begin(), sub_csv_files.end());
            }
        }
    }

    return csv_files;
}


template<typename T>
std::vector<double> linspace(T start_in, T end_in, int num_in)
{

    std::vector<double> linspaced;

    double start = static_cast<double>(start_in);
    double end = static_cast<double>(end_in);
    double num = static_cast<double>(num_in);

    if (num == 0) { return linspaced; }
    if (num == 1)
    {
        linspaced.push_back(start);
        return linspaced;
    }

    double delta = (end - start) / (num - 1);

    for(int i=0; i < num-1; ++i)
    {
        linspaced.push_back(start + delta * i);
    }
    linspaced.push_back(end); // I want to ensure that start and end
    // are exactly the same as the input
    return linspaced;
}

//
// Created by andi on 28.04.23.
//



void exponential_decay(vector<double> &f, vector<double> r, double xi) {
    size_t n = r.size();
    cout << "n = " << n << endl;
    for(int i = 0; i < n; i++){
        f[i] = exp(- r[i] / xi);
    }
}
// that would be to easy, we need one that can handle arrays of type fftw_complex
template <size_t n>
void exponential_decay(fftw_complex f[n], vector<double> r, double xi = 1) {
    for(int i = 0; i < n; i++){
        // TODO this shouldnt work? f[i] should be an array
        f[i][0] = exp(- r[i] / xi);

    }
}

template <int N>
void nanscan(fftw_complex (&f)[N][N]) {
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            if(isnan(f[i][j][0])) {
                cout << "nan found" << endl;
            }
            if(isnan(f[i][j][1])) {
                cout << "nan found" << endl;
            }
        }
    }
    cout << "nan scan done" << endl;
}

void nanscan(const vector<vector<complex<double>>> (&f)) {
    int N = f.size();
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            if(isnan(f[i][j].real())) {
                cout << "nan found" << endl;
            }
            if(isnan(f[i][j].imag())) {
                cout << "nan found" << endl;
            }
        }
    }
    cout << "nan scan done" << endl;
}

string readLastLine(std::ifstream& file) {
    std::string line;
    std::string lastLine;
    while (std::getline(file, line))
    {
        lastLine = line;
    }
    return lastLine;
}


vector<complex<double>> readLastValues(ifstream& file, double& T) {
    string lastline = readLastLine(file);
    vector<complex<double>> values;

    stringstream llss(lastline);
    string token;
    while (std::getline(llss, token, ',')) {
        if(token != "t") {
            complex<double> value = stod(token);
            values.push_back(value);
        }
    }
    // we delete the last value as it is the temperature
    // or do we need the temperature?
    // we might need it but we cannot really give it back
    T = values.back().real();
    values.pop_back();

    // TODO we erase t for now
    values.erase(values.begin());

    return values;
}


vector<vector<complex<double>>> oneD_to_twoD(vector<complex<double>> &q) {
    int N = q.size();
    int lat_dim = (int)sqrt(q.size());


    vector<vector<complex<double>>> q_2D = vector<vector<complex<double>>>( lat_dim, vector<complex<double>>(lat_dim, 0));
    // now i fill it with values
    // actually there was probably a function for this in some package...
    int n;
    int m;

    for(int i = 0; i < N; i++) {
        n = (i % lat_dim);
        m = i / lat_dim;
        // Zeile kreuz spalte?
        q_2D[m][n] = q[i];
    }

    return q_2D;
}

complex<double> two_d_fourier(const vector<vector<complex<double>>> &f, const vector<vector<array<double, 2>>> &q,
                              double pi, double pj) {
    complex<double> ft_ij = 0;
    int lat_dim = f.size();         // f has size of lattice dim

    for(int n = 0; n < lat_dim; n++) {
        for(int m = 0; m < lat_dim; m++) {
            // only one array access
            // performance: we don't need these three variables but this won't make it much faster
            // we could probably do the summation on the gpu for all 250000 lattice sites on the gpu and add them up
            // to ft_ij afterwards
            // this would probably make it very fast
            // but we can also just use this self written stuff to check wheter i am using fftw right
            // but if fftw is still to slow, gpu would save me again i guess.
            double xn = q[m][n][0];
            double ym = q[m][n][1];
            complex<double> f_mn = f[m][n];
            // TODO isnt this lacking a 2pi factor?
            ft_ij += f_mn * exp(1i * (xn * pi + ym * pj));
        }
    }
    return ft_ij;
}


void latticeTrafo(const vector<vector<complex<double>>> &f, const vector<vector<array<double, 2>>> &q,
                  vector<vector<complex<double>>> &ft, vector<vector<array<double, 2>>> &p) {
    // We have x and f, even though with a = 1 we wouldnt really need q but it is more general that way
    // but if we give q, q would be a tuple (x, y) for every lattice site, do we want to do it that way?
    // i mean, we should, right?
    // in our case N=M, do we want to make it even more genral or will we only ever use quadratic lattices?
    // i guess in the moment i don't see the use of a rectangular lattice
    // N is the number of x values i have, so the lattice dim
    // The thing is, we also need the impulses, so we probably have to hand over empty arrays and fill them up here

    int lat_dim = f.size();         // f has size of lattice dim
    int N = lat_dim;
    // To center the impulses around the center, we shift by K = N/2
    int K = N / 2;
    // qi corresponds to xn, so i is the col
    // find out the lattice spacings a_x = x_1 - x_0
    double ax = q[0][1][0] - q[0][0][0];
    double ay = q[1][0][1] - q[0][0][1];

    // now i think we can start to loop?
    // we have to iterate over every lattice site, this will also take very long
    for(int i = 1 - K; i <= N-K; i++) {
        cout << i << endl;
        for(int j = 1 - K; j <= N-K; j++) {
            // 250000 iterations for 500 x 500 lattice
            // i,j are the indexes to calculate the p_i, p_j
            // we need one more pair of indeces to fill up the vectors
            int i_ind = i + K - 1;
            int j_ind = j + K - 1;
            // meaning for j = 0 -> p = 0 we are in the center of our vector
            // calculate the p_i, p_j
            double p_i = 2 * M_PI * i / N / ax;
            double p_j = 2 * M_PI * j / N / ay;
            // calculate the Trafo
            // again 250000, makes again about 10^12 operations which is fucking expensive
            complex<double> ft_ij = two_d_fourier(f, q, p_i, p_j);
            // fill the vectors
            // first the impuls vector, is a pair of  p_i and p_j
            p[j_ind][i_ind] = array<double, 2>{p_i, p_j};
            // now the fouriertransform
            ft[j_ind][i_ind] = ft_ij;
        }
    }
    // okay I think thats it

}

template <int N>
void latticeTrafo1D(const complex<double> (&f)[N][N], const vector<vector<array<double, 2>>> &q,
                    complex<double> (&fx_k)[N], complex<double> (&fy_k)[N], vector<double> &px, vector<double> &py) {
    // To center the impulses around the center, we shift by K = N/2
    int K = N / 2;
    // qi corresponds to xn, so i is the col
    // find out the lattice spacings a_x = x_1 - x_0
    double ax = q[0][1][0] - q[0][0][0];
    double ay = q[1][0][1] - q[0][0][1];

    // so, we need to make clear again what we want to do here
    // What we want, is the 1D Fourier Transform of every row j and save that into an array
    // after that, we can average this over all j and get the Structure Function for k
    // so we have rows:     j=0 : [f_01, f_02, ..., f_0N]
    //                      j=1 : [f_11, f_12, ..., f_1N]
    // and for every row, we do a 1D Fourier Trafo, meaning we get an array:
    //       [f_j1, f_j2, ..., f_jN] -> [px_j(k=-N/2), px_j(k=-N/2 + 1), ..., px_j(k=N/2)]
    // for every j.
    // then we calculate the absolute squared value of this fourier trafo
    // [px_j(k=-N/2) * px_j(k=-N/2)^* , ..., px_j(k=N/2) * px_j(k=N/2)^*] = [fx_j(k=-N/2), ..., fx_j(k=N/2)]
    // afterwards we average/sum
    //                                  [fx_0(k=-N/2), fx_0(k=-N/2 + 1), ..., fx_0(k=N/2)]
    //                                  [fx_1(k=-N/2), fx_1(k=-N/2 + 1), ..., fx_1(k=N/2)]
    //                                  ...
    //  to
    //                                  [fx(k=-N/2), fx(k=-N/2 + 1), ..., fx(k=N/2)]  = S_x(k)
    // We will do the same thing for every column and calculate fy
    for(int j = 0; j < N; j++){
        // This will be the loop that sums um the fourier trafos of the different rows
        complex<double> px_jk[N];
        // we do the same thing for every col
        complex<double> py_jk[N];
        for(int i = 0; i < N; i++) {
            // This will be the loop that does the fourier trafo or row j, so sums over all cols i
            // We need to calculate for all possible k
            for(int k = -N/2; k < N/2; k++) {
                // f[i][j] is in general a complex number, so will be p_ij
                px_jk[k + N/2] += f[i][j] * exp(2 * M_PI * 1i * (double)(i * k));
                // i think the only difference is that we need to switch i and j...
                py_jk[k + N/2] += f[j][i] * exp(2 * M_PI * 1i * (double)(i * k));

                px[k + N/2] = 2 * M_PI * k / N / ax;
                py[k + N/2] = 2 * M_PI * k / N / ay;
            }
        }
        // sum over all j and squared
        fx_k[j] += px_jk[j] * conj(px_jk[j]);
        fy_k[j] += py_jk[j] * conj(py_jk[j]);
        // TODO averaging?
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
            if (i == 0) {
                cout << i << endl;
                cout << "in loop" << p[j][i][0] << endl;
            }
        }
        cout << p[0][i][0] << endl;
    }
}

void fill_p_1d(const vector<double> &q, vector<double> &p) {
    int lat_dim = q.size();         // f has size of lattice dim
    int N = lat_dim;
    // To center the impulses around the center, we shift by K = N/2
    int K = N / 2;
    // qi corresponds to xn, so i is the col
    // find out the lattice spacings a_x = x_1 - x_0
    double a = 1;
    for(int i = 0; i < N; i++) {
        int i_ft = i < K ? i : i - N;
        double p_i = 2 * M_PI * i_ft / N / a;
        p[i] = p_i;
    }
}

template<size_t n>
void write_1d(double ft[n][2], vector<double> &p, double f[n][2], vector<double> &x, string root, string name = "corr-func-trafo.csv") {
    ofstream file;
    create_dir(root);
    file.open(root + name);

    file << "x," << "f," << "p," << "ft\n";

    for(int i = 0; i < n; i++) {
        file << x[i] << ", " << sqrt(f[i][0] * f[i][0] + f[i][1] * f[i][1]) << ", " << p[i] << ", " <<
             sqrt(ft[i][0] * ft[i][0] + ft[i][1] * ft[i][1]) << "\n";
    }

}


template<size_t n>
void write_1d_real(double ft[n][2], vector<double> &p, double f[n][2], vector<double> &x, string root, string name = "corr-func-trafo.csv") {
    ofstream file;
    create_dir(root);
    file.open(root + name);

    file << "x," << "f," << "p," << "ft\n";

    for(int i = 0; i < n; i++) {
        file << x[i] << ", " << f[i][0] << ", " << p[i] << ", " <<
             ft[i][0] << "\n";
    }

}

void print_coords(vector<vector<array<double, 2>>> &q, int val_per_row = 0) {
    size_t lat_dim = q.size();
    if(val_per_row == 0) {
        val_per_row = lat_dim;
    }
    if (val_per_row > lat_dim) {
        val_per_row = lat_dim;
    }
    cout << "[";
    // iterate over the lattice
    for(int m = 0; m < val_per_row; m++) {
        cout << "[";
        for(int n = 0; n < val_per_row; n++) {
            cout << "(" << q[m][n][0] << ", "<< q[m][n][1] << ")   " ;
        }
        cout << "\n ";
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


/*
 * we need another function that plots the Fouriertransform, or at least the squared abs value
 * The structure function is the averaged squared abs, but averaged over what?
 * probably averaged over the all p with same abs value, meaning px + py = const
 * But i could also plot over px and average over py, I think this is what gernot told me
 * Okay we don't plot for the moment since it is somehow insanely hard to integrate matplotlib into cpp
 */
void average_and_write(vector<vector<complex<double>>>& ft, vector<vector<array<double, 2>>> &p, string filename = "./struct.func") {
    // lat dim
    size_t lat_dim = ft.size();

    // so i think we first average
    vector<complex<double>> ft_avg_y(lat_dim, 0.0);
    vector<complex<double>> ft_avg_x(lat_dim, 0.0);

    for(int i = 0; i < lat_dim; i++) {
        // I think we can average simulatenously
        for(int j = 0; j < lat_dim; j++) {
            // avg_y is the average over the y dimension meaning ft[j][i] and j has to run over the lat dim
            ft_avg_y[i] += ft[j][i];
            // avg_x is the average over the x dimension meaning ft[i][j] and j has to run over the lat dim
            ft_avg_x[i] += ft[i][j];
        }
        // finally we also have to get the abs and square
        // and average!!
        ft_avg_y[i] = abs(ft_avg_y[i] * ft_avg_y[i]) / (double)lat_dim;
        ft_avg_x[i] = abs(ft_avg_x[i] * ft_avg_x[i]) / (double)lat_dim;
    }

    // and this is the stuff we plot?
    // no we export
    // open file
    ofstream file;
    file.open(filename);

    // write the stuff in the form:             px      ft_avg_y(px)
    //                                          px_1        avg_1
    //                                          ...         ...
    //                                          py      ft_avg_x(py)
    //                                          ...         ...
    // but we don't actually write the stirings 'px' and 'py' since we will get stress with pandas
    // or maybe we could use another format?    px      ft_avg_y(px)        py          ft_avg_x(py)
    // then we would have a normal heade r
    file << "px," << "ft_avg_y," << "py," << "ft_avg_x\n";
    for(int i = 0; i<lat_dim; i++) {
        // so px(i) is just the p_x value of every entry of p of the i-th col
        // p[0] is first row, p[0][i] is i-th entry of the first row, and p[0][i][0] is px value of the entry
        // px = p[0][i][0]
        // p[i] is the ith-row, p[i][0] is the first entry of the i-th row, which has all the same py values
        // py = p[i][0][1]
        file <<  p[0][i][0] << ", " << abs(ft_avg_y[i]) << ", " << p[i][0][1] << ", " << abs(ft_avg_x[i]);
        if(i < lat_dim - 1) {
            file << "\n";
        }
    }
}

template <size_t N>
void average_and_write(fftw_complex ft[N][N], vector<vector<array<double, 2>>> &p,
                       string filename = "./../../../Generated content/structfunc.fftw") {
    // lat dim

    // so i think we first average

    array<array<double, 2>, N> ft_y{};
    array<array<double, 2>, N> ft_x{};
    array<double, N> ft_avg_x{};
    array<double, N> ft_avg_y{};



    for(int i = 0; i < N; i++) {
        // I think we can average simulatenously
        for(int j = 0; j < N; j++) {
            // avg_y is the average over the y dimension meaning ft[j][i] and j has to run over the lat dim
            ft_y[i][0] += ft[j][i][0];
            ft_y[i][1] += ft[j][i][1];
            // avg_x is the average over the x dimension meaning ft[i][j] and j has to run over the lat dim
            ft_x[i][0] += ft[i][j][0];
            ft_x[i][1] += ft[i][j][1];
/*            if(isnan(ft[j][i][0])) {
                cout << "nan found in ft " << i << " " << j << endl;
            }
            if(isnan(ft[j][i][1])) {
                cout << "nan found in img ft " << i << " " << j << endl;
            }
            if(isnan(ft_x[i][0])) {
                cout << "nan found in ftx real" << endl;
            }
            if(isnan(ft_x[i][1])) {
                cout<< "nan found in ftx img" << endl;
            }*/
        }
        // finally we also have to get the abs and square
        // and average!!
        // okay the average is a bit tricky, we have to average the real and the imaginary part and we actually would
        // have to do it before and we need to average by N * N since that is how many values we hav
        ft_avg_y[i] = (ft_y[i][0] * ft_y[i][0] + ft_y[i][1] * ft_y[i][1]) / pow((double)N, 8);
        ft_avg_x[i] = (ft_x[i][0] * ft_x[i][0] + ft_x[i][1] * ft_x[i][1]) / pow((double)N, 8);
        // cout << "after  " <<  ft_avg_y[i] << endl;
    }

    // and this is the stuff we plot?
    // no we export
    // open file
    ofstream file;
    file.open(filename);

    // write the stuff in the form:             px      ft_avg_y(px)
    //                                          px_1        avg_1
    //                                          ...         ...
    //                                          py      ft_avg_x(py)
    //                                          ...         ...
    // but we don't actually write the stirings 'px' and 'py' since we will get stress with pandas
    // or maybe we could use another format?    px      ft_avg_y(px)        py          ft_avg_x(py)
    // then we would have a normal heade r
    file << "px, " << "ft_avg_y, " << "py, " << "ft_avg_x \n";
    for(int i = 0; i<N; i++) {
        // so px(i) is just the p_x value of every entry of p of the i-th col
        // p[0] is first row, p[0][i] is i-th entry of the first row, and p[0][i][0] is px value of the entry
        // px = p[0][i][0]
        // p[i] is the ith-row, p[i][0] is the first entry of the i-th row, which has all the same py values
        // py = p[i][0][1]
        if(isnan(ft_avg_y[i])) {
            cout << "y nan gefunden " << i << endl;
        }
        if(isnan(ft_avg_x[i])) {
            cout << "x nan gefunden " << i << endl;
        }
        file <<  p[0][i][0] << ", " << ft_avg_y[i] << ", " << p[i][0][1] << ", " << ft_avg_x[i];
        if(i < N - 1) {
            file << "\n";
        }
    }
}

template <size_t N>
void write_lattice_trafo(complex<double> (&ft_x)[N], complex<double> ft_y[N], vector<double> px, vector<double> py,
                         string filename="./../../../Generated content/struct1d.func") {
    ofstream file;
    file.open(filename);

    file << "px, " << "ft_avg_y, " << "py, " << "ft_avg_x \n";
    for(int i = 0; i<N; i++) {
        // write only real part? no i dont think so
        file <<  px[i] << ", " << abs(ft_x[i]) << ", " << py[i] << ", " << abs(ft_y[i]);
        if(i < N - 1) {
            file << "\n";
        }
    }
}


void copy_values(const vector<complex<double>> &values, fftw_complex* &f_fftw) {
    // values and f_fftw have to have the same size
    size_t n = values.size();

    for(int i = 0; i < n; i++) {
        f_fftw[i][0] = values[i].real();
        f_fftw[i][1] = values[i].imag();
    }

}

template<size_t lat_dim>
void copy_values2D(const vector<vector<complex<double>>> &f, fftw_complex (&f_fftw)[lat_dim][lat_dim]) {
    for(int i = 0; i < lat_dim; i++) {
        for(int j = 0; j < lat_dim; j++) {
            f_fftw[i][j][0] = f[i][j].real();
            f_fftw[i][j][1] = f[i][j].imag();
            // cout << f_fftw[i][j][0] << "  ";
        }
    }
}

template<size_t N>
void copy_values2D(const vector<vector<complex<double>>> &f, complex<double> (&f_arr)[N][N]) {
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            f_arr[i][j] = f[i][j];
            // cout << f_fftw[i][j][0] << "  ";
        }
    }
}


class DomainSize {

};



#endif //LEARNINGPROJECT_FOURIERLATTICETRAFO_H
