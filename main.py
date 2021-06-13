import math
import sympy as sp
from sympy.utilities.lambdify import lambdify


def Integration_Simpson_Method(f, n, rng):
    f = lambdify(x, f)
    a, b = rng[0], rng[1]
    h = (b - a) / n
    s = f(a) + f(b)
    X = a
    for i in range(0, n-1):
        X += h
        if i % 2 == 0:
            s += 4*f(X)
        else:
            s += 2*f(X)
        print(str((h/3) * s))
    print("Integration Value = " + str((h/3) * s))


def Bisection_Method(f, little_range, epsilon):
    f = lambdify(x, f)
    a, b = little_range
    k = math.ceil(- math.log(epsilon/(b - a), math.e) / math.log(2, math.e))
    counter = 0

    while abs(b - a) > epsilon:
        c = (a + b) / 2

        if f(a) * f(c) > 0:
            a = c
        else:
            b = c
    counter += 1
    print(c)

    if counter <= k:
        return c, counter


def Newton_Raphson(pol, little_range, epsilon):
    f = lambdify(x, pol)
    df = lambdify(x, sp.diff(pol, x))
    x1 = sum(little_range) / 2
    x2 = x1 - f(x1) / df(x1)

    print(x2)
    counter = 1

    while abs(x2 - x1) > epsilon:
        x1 = x2
        x2 = (x1 - f(x1) / df(x1))

        counter += 1
        print(x2)

    return x2, counter


def roots_Solver(pol, big_range, epsilon, step, method):
    roots_Function_Solver(pol, big_range, epsilon, step, method)
    roots_Derivative_solver(pol, big_range, epsilon, step, method)


def roots_Function_Solver(pol, big_range, epsilon, step, method):
    f = lambdify(x, pol)
    left_bound, right_bound = big_range
    a, b = left_bound, left_bound + step

    while b <= right_bound:

        if f(a) * f(b) < 0:
            solution = method(pol, (a, b), epsilon)
            sol, iterations = solution
            if solution is not None:
                if abs(sol) < epsilon:
                    sol = 0
                print("x = " + str(sol) + ", number of iteration: " + str(iterations))
        a += step
        b += step


def roots_Derivative_solver(pol, big_range, epsilon, step, method):
    f = lambdify(x, pol)
    df = lambdify(x, sp.diff(pol, x))
    left_bound, right_bound = big_range
    a, b = left_bound, left_bound + step

    while b <= right_bound:

        if df(a) * df(b) < 0:

            if abs(df(b)) < epsilon or abs(df(a) < epsilon):
                solution = method(sp.diff(pol, x), (a - step, b + step), epsilon)
                sol, iterations = solution

            else:
                solution = method(sp.diff(pol, x), (a, b), epsilon)
                sol, iterations = solution

            if solution is not None and abs(f(sol)) < epsilon:
                if abs(sol) < epsilon:
                    sol = 0
                print("x = " + str(sol) + ", number of iterations : " + str(iterations))
                solution = None
            else:
                print("Not Converge")
        a += step
        b += step


x = sp.symbols('x')
f = x*math.e**(-x**2 + 5*x)*(2*x**2 - 3*x - 5)
root_rng = (0, 3)
integration_rng = [0.5, 1]
epsilon = 0.0001
step = 0.1
n = 10

print("f(x) = " + "x*math.e**(-x**2 + 5*x)*(2*x**2 - 3*x - 5)")
print("Roots: [0,3]")
print("Bisection Method:")
roots_Solver(f, root_rng, epsilon, step, Bisection_Method)
print("Newton Raphson Method:")
roots_Solver(f, root_rng, epsilon, step, Newton_Raphson)
print("\nIntegration: [0.5,1]")
print("Simpson Method:")
Integration_Simpson_Method(f, 10, integration_rng)



