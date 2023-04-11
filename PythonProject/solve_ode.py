import sympy as sp
from sympy import *


# Define the variables and constants
b, d, C = sp.symbols('b d C')
a, D, m, eta, k = sp.symbols('a D m \eta k', positive=True)
y, x2 = sp.symbols('y x^2', cls=sp.Function)
x, t = sp.symbols('x t')

# Define the differential equation
eq = sp.Eq(y(x).diff(x, 3) + a*y(x).diff(x, 2) + b*y(x).diff(x) + d*y(x), C)
eq = sp.Eq(y(x).diff(x, 2) + a*y(x), 0)
eq = sp.Eq(x2(t).diff(t), (2 * D) - 2 * k / eta * x2(t))
# Find the general solution

x0 = 0
x20 = 0
t0 = 0
y0 = sp.symbols('C2')
dydx = 0

ics = {
    #y(x0) : y0,
    #y(x).diff(x).subs(x, x0): dydx
    x2(t0): x0,
}

sol = sp.dsolve(eq, x2(t), ics=ics)

# Print the solution
print("The general solution is:")
print(sol)
print(latex(sol))