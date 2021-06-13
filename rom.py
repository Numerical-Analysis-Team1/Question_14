import math
import sympy as sp
from sympy.utilities.lambdify import lambdify


def trapezoidal_method(f, n, rng):
    a, b = rng
    h = (b - a) / n
    x = a
    s = 0

    for i in range(n):
        s += (f(x) + f(x + h)) / 2
        x += h

    return s * h


def Romberg_method(f, n, rng):
    r = [[0 for j in range(i)] for i in range(1, n + 1)]

    for i in range(0, n):
        r[i][0] = trapezoidal_method(f, i + 1, rng)

    for i in range(1, n):
        for j in range(1, i + 1):
            r[i][j] = r[i][j - 1] + 1/(4 ** j - 1) * (r[i][j - 1] - r[i - 1][j - 1])

    print(r[n - 1][n - 1])


x = sp.symbols('x')
f = x*math.e**(-x**2 + 5*x)*(2*x**2 - 3*x - 5)
f = lambdify(x, f)

Romberg_method(f, 10, [0.5, 1])
