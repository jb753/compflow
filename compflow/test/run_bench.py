"""Produce benchmarks for the Fortran and native implementations."""
import numpy as np
import timeit
import matplotlib.pyplot as plt
from matplotlib import rcParams
from context import compflow as cf

if __name__ == "__main__":

    # Set up plotting
    rcParams["text.usetex"] = True
    rcParams["font.family"] = "serif"
    rcParams["font.serif"] = "cm"
    rcParams["axes.titlesize"] = "medium"
    rcParams["font.serif"] = "cm"
    tick_marks = 10 ** np.array(np.arange(6))

    # Define input data
    ga = 1.4
    N = np.logspace(0, 5, 14).astype(int)
    reps = 5
    Ma_max = 2.0
    Ma_min = 0.0
    Xmax = 0.1
    Xmin = 1.1

    # FORWARD EVALUATION
    print("Benchmarking forward evaluation...")

    # Loop over array sizes
    dt_native = []
    dt_fort = []
    for Ni in N:

        # Randomise input data
        Ma = np.random.rand(Ni) * (Ma_max - Ma_min) + Ma_min

        # Timers
        T_fort = timeit.Timer(
            "cf.mcpTo_APo_from_Ma(Ma,ga)", "from __main__ import Ma,ga,cf"
        )
        T_native = timeit.Timer(
            "cf.native.mcpTo_APo_from_Ma(Ma,ga)", "from __main__ import Ma,ga,cf"
        )

        # Repeat a number of runs, calculate time per call
        time_per_call = np.empty((2, reps))
        for i in range(reps):
            res_fort = T_fort.autorange()
            res_native = T_native.autorange()
            time_per_call[0, i] = res_fort[1] / res_fort[0]
            time_per_call[1, i] = res_native[1] / res_native[0]

        # Add fastest time to list
        dt_fort.append(time_per_call[0, :].min())
        dt_native.append(time_per_call[1, :].min())

    # Plot
    f, a = plt.subplots()
    f.set_size_inches((4.0, 3.0))
    a.loglog(N, dt_native, label="Native")
    a.loglog(N, dt_fort, label="Fortran")
    a.set_xlabel(r"Array Length, $n$")
    a.set_ylabel(r"Time per call, $\Delta t$/s")
    a.grid(True)
    a.set_title("Benchmark evaluation of $\dot{m}\sqrt{c_pT_0}/Ap_0$")
    a.legend()
    a.set_xlim((1, 1e5))
    a.set_ylim((1e-7, 1e-2))
    a.set_xticks(tick_marks)
    f.tight_layout(pad=0.1)
    plt.savefig("bench_forward.png", dpi=250)

    speedup = np.array(dt_fort) / np.array(dt_native)
    print(
        "Fortran speedup: ",
        1.0
        / speedup[
            (0, -1),
        ],
    )

    # FORWARD EVALUATION
    print("Benchmarking inversion...")

    # Initialise lookup table
    cf.lookup_mcpTo_APo(0.4, ga)

    # Loop over array sizes
    dt_native = []
    dt_fort = []
    dt_lookup = []
    for Ni in N:

        X = np.random.rand(Ni) * (Xmax - Xmin) + Xmin

        # Set up timers
        T_fort = timeit.Timer(
            "cf.Ma_from_mcpTo_APo(X,ga)", "from __main__ import X,ga,cf"
        )
        T_native = timeit.Timer(
            "cf.native.Ma_from_mcpTo_APo(X,ga)", "from __main__ import X,ga,cf"
        )
        T_lookup = timeit.Timer(
            "cf.lookup_mcpTo_APo(X,ga)", "from __main__ import X,ga,cf"
        )

        # Repeat a number of runs, calculate time per call
        time_per_call = np.empty((3, reps))
        for i in range(reps):
            res_fort = T_fort.autorange()
            res_native = T_native.autorange()
            res_lookup = T_lookup.autorange()
            time_per_call[0, i] = res_fort[1] / res_fort[0]
            time_per_call[1, i] = res_native[1] / res_native[0]
            time_per_call[2, i] = res_lookup[1] / res_lookup[0]

        dt_fort.append(time_per_call[0, :].min())
        dt_native.append(time_per_call[1, :].min())
        dt_lookup.append(time_per_call[2, :].min())

    # Make plot
    f, a = plt.subplots()
    f.set_size_inches((4.0, 3.0))
    a.loglog(N, dt_native, label="Native")
    a.loglog(N, dt_fort, label="Fortran")
    a.loglog(N, dt_lookup, label="Lookup")
    a.set_xlabel(r"Array Length, $n$")
    a.set_ylabel(r"Time per call, $\Delta t$/s")
    a.set_title("Benchmark inversion of $\dot{m}\sqrt{c_pT_0}/Ap_0$")
    a.grid(True)
    a.legend()
    a.set_xlim((1, 1e5))
    a.set_xticks(tick_marks)
    a.set_ylim((1e-7, 1e-1))
    f.tight_layout(pad=0.1)
    plt.savefig("bench_inverse.png", dpi=250)

    speedup = np.array(dt_fort) / np.array(dt_native)
    print(
        "Fortran speedup: ",
        1.0
        / speedup[
            (0, -1),
        ],
    )
    speedup = np.array(dt_lookup) / np.array(dt_native)
    print(
        "Lookup speedup: ",
        1.0
        / speedup[
            (0, -1),
        ],
    )

    plt.show()
