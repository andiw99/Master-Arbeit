#include "FourierLatticeTrafo.h"
#include "../Header/Helpfunctions and Classes.h"

using namespace std;

namespace fs = std::filesystem;


template <size_t N>
void average_and_add(fftw_complex ft[N][N], double (&ft_squared_y)[N], double (&ft_squared_x)[N]) {
    // supposed to calculate ft_avg_y for one lattice and add it to ft_avg_y

    array<array<double, 2>, N> ft_y{};
    array<array<double, 2>, N> ft_x{};
    array<double, N> ft_avg_x{};
    array<double, N> ft_avg_y{};


    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            ft_y[i][0] += ft[j][i][0];
            ft_y[i][1] += ft[j][i][1];
            ft_x[i][0] += ft[i][j][0];
            ft_x[i][1] += ft[i][j][1];
        }
        ft_avg_y[i] = (ft_y[i][0] * ft_y[i][0] + ft_y[i][1] * ft_y[i][1]) / pow((double) N, 8);
        ft_avg_x[i] = (ft_x[i][0] * ft_x[i][0] + ft_x[i][1] * ft_x[i][1]) / pow((double) N, 8);
    }
    // now another loop that adds the now calculated squared absolute values to the running squared abs
    for(int i = 0; i < N; i++) {
        ft_squared_y[i] += ft_avg_y[i];
        ft_squared_x[i] += ft_avg_x[i];
    }
    // that should be it
}

template <size_t N>
void trafo_routine(string filepath, double (&ft_squared_y)[N], double (&ft_squared_x)[N]) {
    // needs to take in the references to the ft_avg_y and ft_avg_x and add them up for several files
    // and now it needs to add up everytime for every file, so it needs to take in the filepath
    cout << "reading: " << filepath << endl;
    ifstream file(filepath);
     // check if file is opened
    if(file.is_open()) {
        cout << "File successfully opened" << endl;
    } else {
        cout << "Failed to open file" << endl;
        // abort if file is not opened
        exit(0);
    }
    // and now we do the stuff we did all the time?
    // read last line and put it into a 2D vector
    double T;
    vector<complex<double>> values = readLastValues(file, T);
    vector<vector<complex<double>>> f = oneD_to_twoD(values);

    // init the empty arrays for the fftw
    fftw_complex f_fftw[N][N], ft_fftw[N][N];
    // copy the values from f to f_fftw
    copy_values2D<N>(f, f_fftw);
    // do the trafo
    fftw_plan plan2d;
    plan2d = fftw_plan_dft_2d(N, N, &f_fftw[0][0], &ft_fftw[0][0], FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan2d);

    // now we got to calculate the avg abs squared for every row and col and add them to ft_avg_y and ft_avg_x
    average_and_add(ft_fftw, ft_squared_y, ft_squared_x);
}

int main() {
    // I tried to implement matplotlib with the following flags
    // -I/usr/include/python3.10 -lpython3.10 -Xlinker -export-dynamic -lpthread -lutil -ldl
    // okay so first things firs, we need to read in the csv
    // /home/andi/Documents/Master-Arbeit Code/Generated content/DomainSize/eta=5.00/T=70.00/dt=0.0050/n=62500/alpha=5.00/beta=10.00/J=50.00

    // lattice dim
    const int lat_dim = 100;

    fs::path root = "../../../Generated content/Adaptive Stepsize 2/";
    vector<fs::path> temp_directories = list_dir_paths(root);
    print_vector(temp_directories);


    for(auto path : temp_directories) {
        vector<fs::path> csv_files = list_csv_files(path);
        print_vector(csv_files);

        cout << filesystem::current_path() << endl;
        // we need running arrays for the averages over the lattices
        double ft_squared_y[lat_dim];
        double ft_squared_x[lat_dim];


        for(auto csv_path :csv_files) {
            trafo_routine(csv_path, ft_squared_y, ft_squared_x);
        }
        // and now just average over the run size
        for(int i = 0; i < lat_dim; i++) {
            ft_squared_x[i] /= (double)csv_files.size();
            ft_squared_y[i] /= (double)csv_files.size();
        }

        auto q = init_q(lat_dim);
        auto p = vector<vector<array<double, 2>>>(
                lat_dim, vector<array<double, 2>>(lat_dim, array<double, 2>()));
        fill_p(q, p);


        // write it

        fs::path writepath = path / "struct.fact";
        ofstream ofile;
        ofile.open(writepath);
        ofile << "px, " << "ft_avg_y, " << "py, " << "ft_avg_x \n";
        cout << lat_dim << endl;
        for(int j = 0; j<lat_dim; j++) {
            // so px(i) is just the p_x value of every entry of p of the i-th col
            // p[0] is first row, p[0][i] is i-th entry of the first row, and p[0][i][0] is px value of the entry
            // px = p[0][i][0]
            // p[i] is the ith-row, p[i][0] is the first entry of the i-th row, which has all the same py values
            // py = p[i][0][1]
            int K = lat_dim/2;
            int i = (j + K < lat_dim) ? j + K : j - K;


            ofile <<  p[i][i][0] << ", " << ft_squared_y[i] << ", " << p[i][i][1] << ", " << ft_squared_x[i];
            cout <<  p[i][i][0] << ", " << ft_squared_y[i] << ", " << p[i][i][1] << ", " << ft_squared_x[i] << endl;
            cout << i << "  " << p[0][i][0] << "  " << ft_squared_y[i] << "  " << p[i][0][1] << "   " << ft_squared_x[i] << endl;
            if(j < lat_dim - 1) {
                ofile << endl;
            }
        }
    }


    exit(0);

    const int N = lat_dim;

    auto q = init_q(N);
    auto p = vector<vector<array<double, 2>>>(
            N, vector<array<double, 2>>(N, array<double, 2>()));

    /*
     *

    const int nn = 1000;
    double xi = 1.2;
    string root = "../../Generated content/LatticeTrafo/1D/";
    string name = "1000+20+1.2";
    create_dir(root);
    double end = 20;
    vector<double> r = linspace(0, 20, nn);
    vector<double> pr(nn, 0);
    fill_p_1d(r, pr);
    fftw_complex fu[nn];
    fftw_complex ftu[nn];
    exponential_decay<nn>(fu, r, xi);

    fftw_plan plan;
    plan = fftw_plan_dft_1d(nn, fu, ftu, FFTW_FORWARD, FFTW_ESTIMATE);
    fftw_execute(plan);
    write_1d_real<nn>(ftu, pr, fu, r, root, name);
    exit(0);
*/
    ifstream file("/home/weitze73/Documents/Master-Arbeit/Code/Generated content/SCP/0.csv");
    double T;
    if(file.is_open()) {
        cout << "File successfully opened" << endl;
    } else {
        cout << "Failed to open file" << endl;
        exit(0);
    }
    vector<complex<double>> values = readLastValues(file, T);
    vector<vector<complex<double>>> f = oneD_to_twoD(values);

    // now we can try wether the fourier transform will run through
    // but i guess we should build some markers in so tht we can see where we are
    // we allocate a vector that has the dimensions of f
    cout << f.size() << endl;
    vector<vector<complex<double>>> ft = vector<vector<complex<double>>>(f.size(),
                                                                         vector<complex<double>>(f.size(), 0.0));

    fftw_complex f_fftw[N][N], ft_fftw[N][N];
    complex<double> f_arr[N][N];
    // got to copy the values from 'values' to f_fftw
    copy_values2D<N>(f, f_fftw);
    // check for nans?
    nanscan(f);
    cout << "nan scan done" << endl;
    copy_values2D<N>(f, f_arr);

    auto q1 = init_q(f.size());

    // okay lets try our 1D fourier transform
    // we need arrays for the p values
    vector<double> px(f.size(), 0.0);
    vector<double> py(f.size(), 0.0);
    // We have the f values in fftw format
    // We still need the arrays for the ft_x, ft_y
    complex<double> ft_x[N], ft_y[N];
    // now we can apply?
    latticeTrafo1D(f_arr, q, ft_x, ft_y, px, py);
    // now we write?
    write_lattice_trafo<N>(ft_x, ft_y, px, py);
    auto p1 = vector<vector<array<double, 2>>>(
            f.size(), vector<array<double, 2>>(f.size(), array<double, 2>()));

    print_coords(q, 10);

    // I think my implementation works but it is slow, I will try now to use fftw
    // latticeTrafo(f, q, ft, p);
    // I already wrote a function that gets me the corresponding p
    fill_p(q, p);
    // now i need to do this plan stuff
    // first we need our data in an fftw_complex array
    // i guess for this the size has to be clear at compile time sadly

    /*
     *
    fftw_complex *f_fftw, *ft_fftw;
    // best way is allocating memory with fftw_malloc but i don't now how I would use this to get an array
    // can only be one dimensional and i have to retransform it afterwars, but that is actually not a huge deal?
    f_fftw = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
    ft_fftw = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N * N);
     */

    // initialize this plan, just a 1D trafo now, doesnt make sense. Get back static arrays
    fftw_plan plan2d;
    plan2d = fftw_plan_dft_2d(N, N, &f_fftw[0][0], &ft_fftw[0][0], FFTW_FORWARD, FFTW_ESTIMATE);
    // so now we already have ft_fftw now an could print or write to file?
    fftw_execute(plan2d);
    nanscan(ft_fftw);
    cout << "plan executed " << ft_fftw[0][0] << endl;
    cout << ft_fftw[0][0][0] << endl;
    average_and_write(ft_fftw, p);
    // average_and_write(ft, p);


/*
    for (int i = 0; i < 300; i++) {
        if(i == 250) {
            cout << endl;
        }
        cout << values[i] << "  ";
    }
    cout << endl << "values size: " << values.size() << endl;
    // reorder the vector, this should work as expected?
    cout << endl;
    for (int i = 0; i < 10; i++) {
        for (int j = 0; j < 10; j++) {
            cout << q[i][j] << " ";
        }
        cout << endl;
    }
*/
}