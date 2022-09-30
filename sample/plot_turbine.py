import compflow as cf
import numpy as np
import timeit
import matplotlib.pyplot as plt

# import sympy as sp

# m, g, AR, Yp = sp.symbols('Ma, ga, AR, Yp')
# flowfun = g/sp.sqrt(g-1.)*m*(1.+((g-1.)*m*m)/2.)**((g+1.)/(g-1.)/-2.)
# PR = 1.+Yp*(1.-(1.+((g-1.)*m*m)/2.)**(-g/(g-1.)))


Ma1_in = 0.3
Ma2_in = np.linspace(0.0, 1.0, 1000)
AR_in = 0.39
Yp_in = 0.05
ga_in = 1.4
Ma_tol = 0.000001

f, a = plt.subplots()
Q2b = (
    cf.mcpTo_APo_from_Ma(Ma1_in, ga_in)
    / AR_in
    * (1.0 + Yp_in * (1.0 - 1.0 / cf.Po_P_from_Ma(Ma2_in, ga_in)))
)
Q2a = cf.mcpTo_APo_from_Ma(Ma2_in, ga_in)
a.plot(Ma2_in, Q2a)
a.plot(Ma2_in, Q2b)


def solve_newton(Ma1, ga, AR, Yp):
    Niter = 10
    Ma_guess = np.ones((Niter,)) * np.nan
    K = cf.mcpTo_APo_from_Ma(Ma1, ga) / AR
    Ma_guess[0] = 0.5

    if cf.mcpTo_APo_from_Ma(1.0, ga) < K * (
        1.0 + Yp * (1.0 - 1.0 / cf.Po_P_from_Ma(0.0, ga))
    ):
        return np.nan
    print("beans")

    def f(x):
        return K * (
            1.0 + Yp * (1.0 - 1.0 / cf.Po_P_from_Ma(x, ga))
        ) - cf.mcpTo_APo_from_Ma(x, ga)

    def fp(x):
        return K * Yp * cf.der_Po_P_from_Ma(x, ga) / cf.Po_P_from_Ma(
            x, ga
        ) ** 2.0 - cf.der_mcpTo_APo_from_Ma(x, ga)

    for n in range(Niter - 1):
        Ma_guess[n + 1] = Ma_guess[n] - f(Ma_guess[n]) / fp(Ma_guess[n])
        if np.abs(Ma_guess[n + 1] - Ma_guess[n]) < Ma_tol:
            # Q_guess = cf.mcpTo_APo_from_Ma(Ma_guess,ga)
            # if plot_stuff:
            #     a.plot(Ma_guess,Q_guess,'ko-')
            return Ma_guess[n + 1]
            break
    return np.nan


def solve_fixed_point(Ma1, ga, AR, Yp):
    Niter = 100
    Ma_guess = np.ones((Niter,)) * np.nan
    # Q_guess = np.ones((Niter,))*np.nan
    K = cf.mcpTo_APo_from_Ma(Ma1, ga) / AR
    Ma_guess[0] = 0.5
    for n in range(Niter - 1):
        Ma_guess[n + 1] = (
            K
            * (1.0 + Yp * (1.0 - 1.0 / cf.Po_P_from_Ma(Ma_guess[n], ga)))
            / cf.mcpTo_APo_from_Ma(Ma_guess[n], ga)
            * Ma_guess[n]
        )
        if np.abs(Ma_guess[n + 1] - Ma_guess[n]) < Ma_tol:
            # print('Fixed point broke at %d' % n)
            # if plot_stuff:
            #     Q_guess = cf.mcpTo_APo_from_Ma(Ma_guess,ga)
            #     a.plot(Ma_guess,Q_guess,'kx-')
            return Ma_guess[n + 1]
            break
    return np.nan


T_fixed = timeit.Timer(
    "solve_fixed_point(Ma1_in, ga_in, AR_in, Yp_in)",
    "from __main__ import solve_fixed_point,Ma1_in, ga_in, AR_in, Yp_in",
)
T_newton = timeit.Timer(
    "solve_newton(Ma1_in, ga_in, AR_in, Yp_in)",
    "from __main__ import solve_newton,Ma1_in, ga_in, AR_in, Yp_in",
)
# print('fixed', T_fixed.timeit(number=10000))
print("netwon", T_newton.timeit(number=10000))


print(solve_fixed_point(Ma1_in, ga_in, AR_in, Yp_in))
print(solve_newton(Ma1_in, ga_in, AR_in, Yp_in))
plt.show()
