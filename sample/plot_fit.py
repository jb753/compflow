"""Fit a polynomial to flow function."""
import compflow as cf
import numpy as np
import matplotlib.pyplot as plt
import timeit

# # Lay out the true function
ga = 1.365
N = 100
Ma = np.linspace(0.0, 1.0, 10001)
X = cf.mcpTo_APo_from_Ma(Ma, ga)
poly_sub = np.polyfit(X[Ma < 1.0], Ma[Ma < 1.0], 4)
print(poly_sub)

# X_sub = np.sort(X[Ma<1.])
# fig, ax = plt.subplots()
# ax.plot(np.polyval(poly_sub,X_sub), X_sub, '-')
# ax.plot(Ma, X, 'kx')
# plt.show()


def evaluate_speedup(ga):

    T_datum = timeit.Timer(
        "cf.fortran._Ma_from_mcpTo_APo_old(X,ga_now)",
        "from __main__ import Ma,cf; ga_now = %f; X = cf.mcpTo_APo_from_Ma(Ma, ga_now)"
        % ga,
    )
    calls, time = T_datum.autorange()
    time_per_call_datum = time / calls

    T_datum = timeit.Timer(
        "cf.Ma_from_mcpTo_APo(X,ga_now)",
        "from __main__ import Ma,cf; ga_now = %f; X = cf.mcpTo_APo_from_Ma(Ma, ga_now)"
        % ga,
    )
    calls, time = T_datum.autorange()
    time_per_call_new = time / calls

    return time_per_call_datum / time_per_call_new


gas = np.linspace(1.33, 1.4, 10)
speedup = [evaluate_speedup(gi) for gi in gas]

fig, ax = plt.subplots()
ax.plot(gas, speedup, "k-x")
ax.set_ylim((0.0, 1.5))
ax.set_ylabel("Speedup Factor")
ax.set_xlabel("Specific Heat Ratio, $\gamma$")
plt.show()
