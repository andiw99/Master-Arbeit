#include <iostream>
#include <chrono>
using namespace std;

int right_neighbor_branchless(size_t i, size_t dim_size){
    return (i + 1) % dim_size + (i / dim_size) * dim_size;
}

int right_neighbor_lambda(size_t i, size_t dim_size) {
    return (i % dim_size == dim_size - 1) ? i - (dim_size - 1) : i + 1;
}

int main() {
    // Write C++ code here
    std::cout << "Hello world!" << endl;
    size_t a = 9;
    size_t b = 5;
    cout << a / b << endl;

    for(int i = 0; i < 15; i++) {
        cout << "right neighbor of i = " << i << " is " << right_neighbor_branchless(i, 5) << endl;
        cout << "right neighbor of i = " << i << " is " << right_neighbor_lambda(i, 5) << endl;
    }

    size_t reps = 100000;

    auto startpoint = chrono::high_resolution_clock::now();
    for(int i = 0; i < reps; i++) {
        right_neighbor_branchless(i, 5);
    }
    auto endpoint = chrono::high_resolution_clock::now();
    auto duration = chrono::duration_cast<chrono::microseconds>(endpoint - startpoint);
    cout << "Branchless took " << duration.count() << " us" << endl;

    startpoint = chrono::high_resolution_clock::now();
    for(int i = 0; i < reps; i++) {
        right_neighbor_lambda(i, 5);
    }
    endpoint = chrono::high_resolution_clock::now();
    duration = chrono::duration_cast<chrono::microseconds>(endpoint - startpoint);
    cout << "Lambda took " << duration.count() << " us" << endl;

    return 0;
}