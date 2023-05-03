//
// Created by weitze73 on 27.04.23.
//

#include "CalcDomainSize.h"
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <iostream>
#include <unordered_map>
#include <cmath>
#include <chrono>
#include <cstdlib>
#include <complex>


using namespace std;

string readLastLine(std::ifstream& file) {
    std::string line;
    std::string lastLine;
    while (std::getline(file, line))
    {
        lastLine = line;
    }
    return lastLine;
}

vector<double> readLastValues(ifstream& file, double& T) {
    string lastline = readLastLine(file);
    vector<double> values;

    stringstream llss(lastline);
    string token;
    while (std::getline(llss, token, ',')) {
        if(token != "t") {
            double value = stod(token);
            values.push_back(value);
        }
    }
    // we delete the last value as it is the temperature
    // or do we need the temperature?
    // we might need it but we cannot really give it back
    T = values.back();
    values.pop_back();

    return values;
}

/*
 * Calculates the correlation function for the vector of q-values 'values'
 * It is preliminary that the values vector only consists of q values that are made to an 1D array
 */
unordered_map<int, double> calcCorrFuncManhatten(vector<double> values) {
    // the lattice size is the sqrt of the size of values
    int n = values.size();
    int lat_dim = (int)sqrt(values.size());
    // TODO do we need to calculate the average value to calculate the correlation function?
    // for now not...
    // Okay, since we are looking at the manhatten metrik the largest possible distance is (lat_dim - 1) / 2 + (lat_dim - 1) / 2
    //  *   *   *   *   *
    //                  |
    //  *   *   *   *   *
    //                  |
    //  *   *   * ->* ->*

    //  *   *   *   *   *

    //  *   *   *   *   *
    int d_max = lat_dim - 1;
    // my first approach was to just go over every lattic site, calculate the distance to every lattice site and the
    // correlation and save it in the map, but that would make for 1 million lattice sites 10^6 * 10^6 operations
    // which would be pretty large?
    // i mean we can try an potentially even let it run on gpu later but isn't there a better way?
    // we could probably cut the operations in half if we would not calculate the interaction between i an j twice
    // i guess if we calculated all values for i=1 with j = 2,...,10^6, in the next step in only have to calculate the
    // Interaction between i = 2 and j = 3, ..., 10^6 since the Interaction i = 1, j=2 is already covered

    // we need a map that maps the distance to the value of the correlation function
    unordered_map<int, double> corrFunc;
    // we need a second map to count the calculations that we did for one specific distance to average later
    unordered_map<int, int> countCalcs;

    // okay we need to optimize this code, it takes forever
    // most obvious thing would be to initialized the needed varialbes outside to loop, but i don't know if the
    // compiler doesn't do that already.
    int i_hor;
    int j_hor;
    int d_horizontal;
    int i_ver;
    int j_ver;
    int d_vertical;
    int d;
    int d_hor_naive;
    int d_ver_naive;

    // if we want every interaction to be calculated, we can't do less loops
    // but we can try to minimize the operations inside the loop
    // instead of doing two calculations, we could use lambdas or if statements
    auto start = chrono::high_resolution_clock::now();
    for(int i = 0; i < n; i++) {
        // first index, runs over every lattice sit
        i_hor = (i % lat_dim);
        i_ver = i / lat_dim;
        long duration_pre_calcs = 0;
        long duration_map_stuff = 0;


        for(int j = i; j < n; j++) {
            // second index, only runs over the following lattice sites after i
            // Now we calculate the distance between i and j
            // the horizonal difference is just the difference of the index values modulo the lattice sites
            // but because of PBC, we have to take the minimum of that and the addition of the two
            auto loop_start = chrono::high_resolution_clock::now();
            j_hor = (j % lat_dim);
            // possibility 1:
            // d_horizontal = min(abs(i_hor - j_hor), (i_hor + j_hor) % lat_dim);
            // possibility 2:
            d_hor_naive = abs(i_hor - j_hor);
            d_horizontal = (d_hor_naive >= d_max) ? (i_hor + j_hor) % lat_dim : d_hor_naive;
            // the vertical distance is the difference between the indexes of the rows
            // to get the index of the row, we devide the index by n and round down. Integer division?
            j_ver = j / lat_dim;
            // we know that j is always greater than i, so we can save us the abs calculation?
            // possibility 1
            // d_vertical = min(j_ver - i_ver, (j_ver + i_ver) % lat_dim);
            // pssibility 2
            d_ver_naive = abs(i_ver - j_ver);
            d_vertical = (d_ver_naive >= d_max) ? (i_ver + j_ver) % lat_dim : d_ver_naive;
            // in manhatten metrik, the whole distance is just the addition of both
            d = d_horizontal + d_vertical;
            auto end_pre_calcs = chrono::high_resolution_clock::now();
            auto duration_pre_calc = chrono::duration_cast<std::chrono::nanoseconds>(end_pre_calcs - loop_start);
            duration_pre_calcs += duration_pre_calc.count();
            // does it work like that or do we have to initialize the value to be zero?
            auto start_map_stuff = chrono::high_resolution_clock::now();
            corrFunc[d] += values[i] * values[j];
            countCalcs[d] += 1;
            auto end_map_stuff = chrono::high_resolution_clock::now();
            auto duration_map = chrono::duration_cast<std::chrono::nanoseconds>(end_map_stuff - start_map_stuff);
            duration_map_stuff += duration_map.count();
        }
        if ( i % 1000 == 0) {
            cout << "i = " << i << endl;
            auto end = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<std::chrono::seconds>(end - start);
            cout << "execution took for last 1000 sites " << duration.count() << "s" << endl;
            cout << "execution took for pr calcs " << duration_pre_calcs / 1000 << "ms" << endl;
            cout << "execution took for map stuff " << duration_map_stuff / 1000 << "ms" << endl;
            start = chrono::high_resolution_clock::now();
        }
    }
    // averaging
    // could be that we want an ordered map later, but lets see for now if we can do it with the unordered ones
    for (auto kv : corrFunc) {
        corrFunc[kv.first] /= countCalcs [kv.first];
    }

    return corrFunc;
}


std::vector<double> calcCorrFuncManhattenVector(const vector<double> &values) {
    // the lattice size is the sqrt of the size of values
    const int n = values.size();
    const int lat_dim = (int)sqrt(values.size());
    // TODO do we need to calculate the average value to calculate the correlation function?
    // for now not...
    // Okay, since we are looking at the manhatten metrik the largest possible distance is (lat_dim - 1) / 2 + (lat_dim - 1) / 2
    //  *   *   *   *   *
    //                  |
    //  *   *   *   *   *
    //                  |
    //  *   *   * ->* ->*

    //  *   *   *   *   *

    //  *   *   *   *   *
    const int d_max = lat_dim - 1;

    vector<double> corrFunc(lat_dim, 0);
    vector<double> countCalcs(lat_dim, 0);


    int i_hor;
    int j_hor;
    int d_horizontal;
    int i_ver;
    int j_ver;
    int d_vertical;
    int d;
    int d_hor_naive;
    int d_ver_naive;


    // if we want every interaction to be calculated, we can't do less loops
    // but we can try to minimize the operations inside the loop
    // instead of doing two calculations, we could use lambdas or if statements
    for(int i = 0; i < n; i++) {
        // first index, runs over every lattice sit
        i_hor = (i % lat_dim);
        i_ver = i / lat_dim;
        for(int j = i; j < n; j++) {
            j_hor = (j % lat_dim);
            d_hor_naive = abs(i_hor - j_hor);
            d_horizontal = (d_hor_naive >= d_max) ? (i_hor + j_hor) % lat_dim : d_hor_naive;
            j_ver = j / lat_dim;
            d_ver_naive = abs(i_ver - j_ver);
            d_vertical = (d_ver_naive >= d_max) ? (i_ver + j_ver) % lat_dim : d_ver_naive;
            d = d_horizontal + d_vertical;
            corrFunc[d] += values[i] * values[j];
            countCalcs[d] += 1;

        }
        if ( i % 1000 == 0) {
            cout << "i = " << i << endl;
        }
    }
    // averaging
    // could be that we want an ordered map later, but lets see for now if we can do it with the unordered ones
    cout << "d      C(d)" << endl;
    for (int i = 0; i < lat_dim; i++) {
        corrFunc[i] /= countCalcs[i];
        cout << i << "   " << corrFunc[i] << endl;
    }
    cout << "How can that be?" << endl;
    cout << corrFunc.size() << endl;

    return corrFunc;
}

std::vector<double> getVector(int n) {
    std::vector<double> result(n, 0.0);
    for (int i = 0; i < n; i++) {
        result[i] = 1.1 * i;
    }
    return result;
}


template <size_t lat_dim>
array<double, lat_dim> calcCorrFuncManhattenArray(const vector<double> &values) {
    // the lattice size is the sqrt of the size of values
    const int n = lat_dim*lat_dim;
/*    cout << "values size " << values.size() << endl;
    cout << values[13000] << endl;
    cout << values[13001] << endl;*/
    // TODO do we need to calculate the average value to calculate the correlation function?
    // for now not...
    // Okay, since we are looking at the manhatten metrik the largest possible distance is (lat_dim - 1) / 2 + (lat_dim - 1) / 2
    //  *   *   *   *   *
    //                  |
    //  *   *   *   *   *
    //                  |
    //  *   *   * ->* ->*

    //  *   *   *   *   *

    //  *   *   *   *   *
    const int d_max = lat_dim - 1;
    // my first approach was to just go over every lattic site, calculate the distance to every lattice site and the
    // correlation and save it in the map, but that would make for 1 million lattice sites 10^6 * 10^6 operations
    // which would be pretty large?
    // i mean we can try an potentially even let it run on gpu later but isn't there a better way?
    // we could probably cut the operations in half if we would not calculate the interaction between i an j twice
    // i guess if we calculated all values for i=1 with j = 2,...,10^6, in the next step in only have to calculate the
    // Interaction between i = 2 and j = 3, ..., 10^6 since the Interaction i = 1, j=2 is already covered

    // we need a map that maps the distance to the value of the correlation function
    // we need a second map to count the calculations that we did for one specific distance to average later

    // honestly the first one took forever, mainly because there are so many calculations to do but secondly also
    // because maps are fucking slow. lets try to use std arrays
    // the good thing is, since in manhatten metrik we have all distances from 0 to d_max present, we just initalize
    // an array of that size
    // okay the thing is the size of this array has to be clear at compile time. This would be uncomfortable since
    // I would have to know th size and hardcode it all the time (since i am reading a file)
    // we try vectors first but i want to compare the two speedwise
    array<double, lat_dim> corrFunc{};
    array<double, lat_dim> countCalcs{};

    // okay we need to optimize this code, it takes forever
    // most obvious thing would be to initialized the needed varialbes outside to loop, but i don't know if the
    // compiler doesn't do that already.
    int i_hor;
    int j_hor;
    int d_horizontal;
    int i_ver;
    int j_ver;
    int d_vertical;
    int d;
    int d_hor_naive;
    int d_ver_naive;


    // if we want every interaction to be calculated, we can't do less loops
    // but we can try to minimize the operations inside the loop
    // instead of doing two calculations, we could use lambdas or if statements
    auto start = chrono::high_resolution_clock::now();
    for(int i = 0; i < n; i++) {
        // first index, runs over every lattice sit
        i_hor = (i % lat_dim);
        i_ver = i / lat_dim;
        long duration_pre_calcs = 0;
        long duration_map_stuff = 0;

        for(int j = i; j < n; j++) {
            cout << j << endl;
            // second index, only runs over the following lattice sites after i
            // Now we calculate the distance between i and j
            // the horizonal difference is just the difference of the index values modulo the lattice sites
            // but because of PBC, we have to take the minimum of that and the addition of the two
            j_hor = (j % lat_dim);
            // possibility 1:
            // d_horizontal = min(abs(i_hor - j_hor), (i_hor + j_hor) % lat_dim);
            // possibility 2:
            d_hor_naive = abs(i_hor - j_hor);
            d_horizontal = (d_hor_naive >= d_max) ? (i_hor + j_hor) % lat_dim : d_hor_naive;
            // the vertical distance is the difference between the indexes of the rows
            // to get the index of the row, we devide the index by n and round down. Integer division?
            j_ver = j / lat_dim;
            // we know that j is always greater than i, so we can save us the abs calculation?
            // possibility 1
            // d_vertical = min(j_ver - i_ver, (j_ver + i_ver) % lat_dim);
            // pssibility 2
            d_ver_naive = abs(i_ver - j_ver);
            d_vertical = (d_ver_naive >= d_max) ? (i_ver + j_ver) % lat_dim : d_ver_naive;
            // in manhatten metrik, the whole distance is just the addition of both
            d = d_horizontal + d_vertical;

            // Hier passiert genau bei j=13000 ein SIGSEGV, der im run vorher noch nicht passiert
            cout << "before array stuff" << endl;
            cout << "values[13000] = " << values[13000] << endl;
            cout << "values[j+1] = " << values[j+1] << endl;
            cout << "values[j] = " << values[j] << endl;
            corrFunc[d] += values[i] * values[j];
            countCalcs[d] += 1;

        }
        if ( i % 1000 == 0) {
            cout << "i = " << i << endl;
            auto end = chrono::high_resolution_clock::now();
            auto duration = chrono::duration_cast<std::chrono::seconds>(end - start);
            cout << "execution took for last 1000 sites " << duration.count() << "s" << endl;
            cout << "execution took for pr calcs " << duration_pre_calcs / 1000 << "ms" << endl;
            cout << "execution took for vector stuff " << duration_map_stuff / 1000 << "ms" << endl;
            start = chrono::high_resolution_clock::now();
        }
    }
    // averaging
    // could be that we want an ordered map later, but lets see for now if we can do it with the unordered ones
    for (int i = 0; i < lat_dim; i++) {
        corrFunc[i] /= countCalcs[i];
    }

    return corrFunc;
}

/*
 * Okay we write a 2D DFT on our own to understand what happens, afterwards we can still use the fftw
 * Even though i don't think the DFT will be the part of the program that will be time consuming but lets see
 * Okay lets think about what we have, what needs to come in and what we want to come out.
 * We have the 2D lattice in a 1D array. For modularity, I think i might should write a function that first puts the
 * 1D Vector in a 2D array and then a 2D DFT Function that transforms this lattice
 * For now the lattice spacing in every direction is 1, so just the 2D array has all information over the lattice
 */

/*
 * Okay so what are we working with? vector<double, n> to vector<vector<double, lat_dim>, lat_dim>?
 * We can again not use arrays since their memory has to be known at compile time.
 * I guess vector worked fine the last times...
 */

vector<vector<double>> oneD_to_twoD(vector<double> &q) {
    int N = q.size();
    int lat_dim = (int)sqrt(q.size());

    vector<vector<double>> q_2D = vector<vector<double>>( lat_dim, vector<double>(0, lat_dim));
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


complex<double> two_d_fourier(const vector<vector<complex<double>>> &f, const vector<vector<pair<double, double>>> &q,
                     double pi, double pj) {
    complex<double> ft_ij = 0;
    int lat_dim = f.size();         // f has size of lattice dim

    for(int n = 0; n < lat_dim; n++) {
        for(int m = 0; m < lat_dim; m++) {
            // only one array access
            pair<double, double> xy = q[m][n];
            // performance: we don't need these three variables but this won't make it much faster
            // we could probably do the summation on the gpu for all 250000 lattice sites on the gpu and add them up
            // to ft_ij afterwards
            // this would probably make it very fast
            // but we can also just use this self written stuff to check wheter i am using fftw right
            // but if fftw is still to slow, gpu would save me again i guess.
            double xn = xy.first;
            double ym = xy.second;
            complex<double> f_mn = f[m][n];
            ft_ij += f_mn * exp(1i * (xn * pi + ym * pj));
        }
    }
    return ft_ij;
}


void latticeTrafo(const vector<vector<complex<double>>> &f, const vector<vector<pair<double, double>>> &q,
                                        vector<vector<complex<double>>> &ft, vector<vector<pair<double, double>>> &p) {
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
    double ax = q[0][1].first - q[0][0].first;
    double ay = q[1][0].second - q[0][0].second;

    // now i think we can start to loop?
    // we have to iterate over every lattice site, this will also take very long
    for(int i = 1 - K; i <= N-K; i++) {
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
            p[j_ind][i_ind] = pair<double, double>(p_i, p_j);
            // now the fouriertransform
            ft[j_ind][i_ind] = ft_ij;
        }
    }
    // okay I think thats it
}

int main() {
    // okay so first things firs, we need to read in the csv
    ifstream file("DomainSize/1.csv");
    vector<double> values;
    double T;

    // reorder the vector, this should work as expected?
    vector<vector<double>> q = oneD_to_twoD(values);
    for(int i = 0; i < 300; i ++ ) {
        cout << values[i] << "  ";
    }
    cout << endl;
    for(int i = 0; i < 250; i ++ ) {
        for(int j = 0; j < 250; j++) {
            cout << q[i][j] << " ";
        }
        cout << endl;
    }


    array<double, 500> corrfuncarray{};
    values = readLastValues(file, T);
    vector<double> corrfuncvector = calcCorrFuncManhattenVector(values);
    cout << "I am not exiting somewhere am I?" << endl;
    for(int i = 0; i < corrfuncvector.size(); i++) {
        cout << i << ", " << corrfuncvector[i] << endl;
    }

/*
    unordered_map<int, double> corrfunc;
    auto start = chrono::high_resolution_clock::now();

    corrfunc = calcCorrFuncManhatten(values);

    auto end = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<std::chrono::seconds>(end - start);
    cout << "execution took " << duration.count() << "s" << endl;

    for (auto kv : corrfunc) {
        cout << kv.first << ": " << kv.second << endl;
    }
*/


    // okay so we got all values


    return 0;
}