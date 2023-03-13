#include <iostream>
#include <chrono>

long fibonacci(long number) {
    if (number < 3) {
        return 1;
    } else {
        return fibonacci(number - 1) + fibonacci(number - 2);
    }
}

size_t py_mod2(int a, int b) {
    // Implementation einer Funktion die sich verhält wie die python modulo funktion
    std::cout << b << std::endl;
    std::cout << a << std::endl;
    std::cout << (a % b) << std::endl;
    std::cout << (b + (a % b)) << std::endl;
    std::cout << ((b + (a % b)) % b) << std::endl;
    return ((b + (a % b)) % b);
}

int main() {
    std::cout << (-1) % 10 << std::endl;
    std::cout << (10 + ((-1) % 10)) % 10 << std::endl;
    std::cout << py_mod2((-1), 20) << std::endl;
    std::size_t a = 5;
    int b = 4;
    std::cout << (b < a);


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

