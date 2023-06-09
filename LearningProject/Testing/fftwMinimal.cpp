//
// Created by andi on 09.06.23.
//

#include "fftwMinimal.h"
#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>


int main()
{
    int N = 10;

    fftw_complex *in, *out;
    fftw_plan p;


    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    p = fftw_plan_dft_1d(N, in, out, FFTW_FORWARD, FFTW_ESTIMATE);


    fftw_execute(p); /* repeat as needed */


    fftw_destroy_plan(p);
    fftw_free(in); fftw_free(out);
}
