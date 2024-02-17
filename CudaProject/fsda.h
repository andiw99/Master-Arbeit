#include <iostream>
#include <cufft.h>
#include <cufftXt.h>
#include <cublas_v2.h>
#include <thrust/device_vector.h>

// Function to calculate the autocorrelation using Fourier transform
void autocorrelation(const thrust::device_vector<double>& signal, thrust::device_vector<double>& result)
{
    // Size of the signal
    int N = signal.size();

    // Allocate device memory for the Fourier transform
    thrust::device_vector<cufftDoubleComplex> fft_signal(N);

    // Plan for cuFFT
    cufftHandle plan;
    cufftPlan1d(&plan, N, CUFFT_Z2Z, 1);

    // Perform forward FFT
    cufftExecZ2Z(plan, reinterpret_cast<cufftDoubleComplex*>(thrust::raw_pointer_cast(signal.data())),
                 reinterpret_cast<cufftDoubleComplex*>(thrust::raw_pointer_cast(fft_signal.data())),
                 CUFFT_FORWARD);

    // Compute the power spectral density by multiplying the Fourier transform by its complex conjugate
    thrust::transform(fft_signal.begin(), fft_signal.end(), fft_signal.begin(), fft_signal.begin(), thrust::multiplies<cufftDoubleComplex>());

    // Plan for cuFFT
    cufftHandle inv_plan;
    cufftPlan1d(&inv_plan, N, CUFFT_Z2Z, 1);

    // Perform inverse FFT
    cufftExecZ2Z(inv_plan, reinterpret_cast<cufftDoubleComplex*>(thrust::raw_pointer_cast(fft_signal.data())), reinterpret_cast<cufftDoubleComplex*>(thrust::raw_pointer_cast(result.data())), CUFFT_INVERSE);

    // Scale the result by 1/N
    thrust::transform(result.begin(), result.end(), result.begin(), thrust::bind1st(thrust::multiplies<double>(), 1.0 / N));

    // Destroy plans
    cufftDestroy(plan);
    cufftDestroy(inv_plan);
}

int main()
{
    // Example signal (replace with your own data)
    thrust::device_vector<double> signal = {1.0, 2.0, 3.0, 4.0, 5.0};

    // Autocorrelation result
    thrust::device_vector<double> autocorr(signal.size());

    // Calculate autocorrelation
    autocorrelation(signal, autocorr);

    // Print the autocorrelation result
    std::cout << "Autocorrelation function: ";
    for (auto value : autocorr)
    {
        std::cout << value << " ";
    }
    std::cout << std::endl;

    return 0;
}
