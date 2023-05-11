//
// Created by weitze73 on 11.05.23.
//
#include "FourierLatticeTrafo.h"


void calc_corr(vector<vector<complex<double>>> &f, vector<double> &C_x, vector<double> &C_y) {
    // i guess we could also use just normal doubles, would also give some performance
    // but i actually really don't think that it will be very computationally difficult
    int lat_dim = f.size();
    // I don't think that we have to pay attention to the lattice spacings. If we set them one we are just measuring
    // the distance in multiples of the lattice spacing
    // because of PBC, the maximum distance is half the lat_dim
    int d_max = lat_dim / 2;

    // loop over the rows(cols)
    for(int i = 0; i < lat_dim; i++) {
        // for every row we loop over all possible distances, which are 0, 1, ..., d_max
        for(int d = 0; d <= d_max; d++) {
            // for every distance we loop over every lattice site, so over every col (row)
            for(int j = 0; j < lat_dim; j++) {
                // we only have to add up the +d terms since we have pbc?
                C_x[d] += f[i][j].real() * f[i][(j + d) % lat_dim].real();
                C_y[d] += f[j][i].real() * f[(j + d) % lat_dim][i].real();
            }
        }
    }
    // meaning now for C_x[d] we have had lat_dim(i) * lat_dim(j) additions, so normalization
    for(int d = 0; d <= d_max; d++) {
        C_x[d] /= (lat_dim * lat_dim);
        C_y[d] /= (lat_dim * lat_dim);
    }

}


void write_corr_func(vector<double> &C_x, vector<double> C_y, string writepath) {
    ofstream file;
    file.open(writepath + ".csv");

    file << "C_x,C_y" << endl;
    for(int i = 0; i < C_x.size(); i++) {
        file << C_x[i] << "," << C_y[i];
        if(i < C_x.size() - 1) {
            file << endl;
        }
    }
}


int main() {
    // We still need to read in the values
    string readpath = "../../Generated content/Quadratic/0.csv";
    string writepath = "../../Generated content/Quadratic/corr.func";

    cout << "reading: " << readpath << endl;
    ifstream file(readpath);
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

    // wanna work with arrays or vectors? I think since we only calculate C_x and C_y, we can use vectors
    // initialize vectors for those
    int nr_dists = (int)f.size() / 2 + 1;
    cout << "number of possible distances: " << nr_dists << endl;
    vector<double> C_x(nr_dists, 0.0);
    vector<double> C_y(nr_dists, 0.0);
    // calculate the corr func and save them into those vectors
    calc_corr(f, C_x, C_y);
    // write them to a file
    write_corr_func(C_x, C_y, writepath);
    return 0;
}
