//
// Created by andi on 09.06.23.
//

#include "fftwMinimal.h"
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <iostream>
#include <array>
#include <cmath>
#include <filesystem>

using namespace std;

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

int main()
{
    int N = 10;

    fftw_complex *in, *out;
    fftw_plan p;


    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);


    fftw_execute(p); /* repeat as needed */

    fftw_complex a[3];

    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
    cout << "Testing minimal 4" << endl;


    cout << time(NULL) << endl;
}
