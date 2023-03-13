#=
main:
- Julia version: 
- Author: weitze73
- Date: 2023-03-07
=#
using BenchmarkTools
using Printf

BenchmarkTools.DEFAULT_PARAMETERS.samples = 1

fibonnaci(n) =
    if n < 3
        return 1
    else
        return fibonnaci(n - 2) + fibonnaci(n - 1)
    end

print("Welche Fibonacci Zahl möchtest du wissen? ")
n = readline()
n = parse(Int64, n)
b = @benchmarkable @printf("Die %d. Fibonacci Zahl ist %d \n", n, fibonnaci(n)) samples=1
print(run(b))