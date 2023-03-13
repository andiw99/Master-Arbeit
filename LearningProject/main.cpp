#include <iostream>
#include <chrono>

long fibonacci(long number) {
    if (number < 3) {
        return 1;
    } else {
        return fibonacci(number - 1) + fibonacci(number - 2);
    }
}

/*
int main() {
    std::cout << "Welche Stelle der Fibonacci Reihe möchtest du wissen? ";
    long number;
    std::cin >> number;
    auto time_before_calc = std::chrono::high_resolution_clock::now();
    std::cout << "Die " << number << ". Fibonacci Zahl ist " << fibonacci(number) << "." << std::endl;
    auto time_after_calc = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            time_after_calc - time_before_calc);
    std::cout << "Die Berechnung dauerte " << duration.count() << "ms.";

    return 0;
}
*/
