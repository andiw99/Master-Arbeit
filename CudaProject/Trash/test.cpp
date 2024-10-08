/*
 * Copyright 2020 NVIDIA Corporation.  All rights reserved.
 *
 * NOTICE TO LICENSEE:
 *
 * This source code and/or documentation ("Licensed Deliverables") are
 * subject to NVIDIA intellectual property rights under U.S. and
 * international Copyright laws.
 *
 * These Licensed Deliverables contained herein is PROPRIETARY and
 * CONFIDENTIAL to NVIDIA and is being provided under the terms and
 * conditions of a form of NVIDIA software license agreement by and
 * between NVIDIA and Licensee ("License Agreement") or electronically
 * accepted by Licensee.  Notwithstanding any terms or conditions to
 * the contrary in the License Agreement, reproduction or disclosure
 * of the Licensed Deliverables to any third party without the express
 * written consent of NVIDIA is prohibited.
 *
 * NOTWITHSTANDING ANY TERMS OR CONDITIONS TO THE CONTRARY IN THE
 * LICENSE AGREEMENT, NVIDIA MAKES NO REPRESENTATION ABOUT THE
 * SUITABILITY OF THESE LICENSED DELIVERABLES FOR ANY PURPOSE.  IT IS
 * PROVIDED "AS IS" WITHOUT EXPRESS OR IMPLIED WARRANTY OF ANY KIND.
 * NVIDIA DISCLAIMS ALL WARRANTIES WITH REGARD TO THESE LICENSED
 * DELIVERABLES, INCLUDING ALL IMPLIED WARRANTIES OF MERCHANTABILITY,
 * NONINFRINGEMENT, AND FITNESS FOR A PARTICULAR PURPOSE.
 * NOTWITHSTANDING ANY TERMS OR CONDITIONS TO THE CONTRARY IN THE
 * LICENSE AGREEMENT, IN NO EVENT SHALL NVIDIA BE LIABLE FOR ANY
 * SPECIAL, INDIRECT, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, OR ANY
 * DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
 * WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS
 * ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE
 * OF THESE LICENSED DELIVERABLES.
 *
 * U.S. Government End Users.  These Licensed Deliverables are a
 * "commercial item" as that term is defined at 48 C.F.R. 2.101 (OCT
 * 1995), consisting of "commercial computer software" and "commercial
 * computer software documentation" as such terms are used in 48
 * C.F.R. 12.212 (SEPT 1995) and is provided to the U.S. Government
 * only as a commercial end item.  Consistent with 48 C.F.R.12.212 and
 * 48 C.F.R. 227.7202-1 through 227.7202-4 (JUNE 1995), all
 * U.S. Government End Users acquire the Licensed Deliverables with
 * only those rights set forth herein.
 *
 * Any use of the Licensed Deliverables in individual and commercial
 * software must include, in the user documentation and internal
 * comments to the code, the above Disclaimer and U.S. Government End
 * Users Notice.
 */

#include <array>
#include <complex>
#include <iostream>
#include <random>
#include <vector>
#include <hip/hip_runtime.h>
#include <hipfft/hipfft.h>


using dim_t = std::array<int, 3>;

int main(int argc, char *argv[]) {
    hipfftHandle plan;
    hipStream_t stream = NULL;

    int n = 2;
    dim_t fft = {n, n, n};
    int batch_size = 2;
    int fft_size = fft[0] * fft[1] * fft[2];

    using scalar_type = float;
    using data_type = std::complex<scalar_type>;

    std::vector<data_type> data(fft_size * batch_size);

    for (int i = 0; i < data.size(); i++) {
        data[i] = data_type(i, -i);
    }

    std::printf("Input array:\n");
    for (auto &i : data) {
        std::printf("%f + %fj\n", i.real(), i.imag());
    }
    std::printf("=====\n");

    hipfftComplex *d_data = nullptr;

    hipfftCreate(&plan);
    // inembed/onembed being nullptr indicates contiguous data for each batch, then the stride and dist settings are ignored
    hipfftPlanMany(&plan, fft.size(), fft.data(),
                             nullptr, 1, 0, // *inembed, istride, idist
                             nullptr, 1, 0, // *onembed, ostride, odist
                             HIPFFT_C2C, batch_size);

    hipStreamCreateWithFlags(&stream, hipStreamNonBlocking);
    hipfftSetStream(plan, stream);

    // Create device data arrays
    hipMalloc(reinterpret_cast<void **>(&d_data), sizeof(data_type) * data.size());
    hipMemcpyAsync(d_data, data.data(), sizeof(data_type) * data.size(),
                                 hipMemcpyHostToDevice, stream);

    /*
     * Note:
     *  Identical pointers to data and output arrays implies in-place transformation
     */
    hipfftExecC2C(plan, d_data, d_data, HIPFFT_FORWARD);
    hipMemcpyAsync(data.data(), d_data, sizeof(data_type) * data.size(),
                                 hipMemcpyDeviceToHost, stream);
    hipStreamSynchronize(stream);

    std::printf("Output array after Forward transform:\n");
    for (auto &i : data) {
        std::printf("%f + %fj\n", i.real(), i.imag());
    }
    std::printf("=====\n");


    /* free resources */
    hipFree(d_data);
    hipfftDestroy(plan);
    hipStreamDestroy(stream);
    hipDeviceReset();

    return EXIT_SUCCESS;
}
