"""Produce benchmarks for a function."""
import numpy as np
import timeit
import matplotlib.pyplot as plt
from matplotlib import rcParams
from context import compflow as cf

if __name__ == '__main__':

    # Set up plotting
    rcParams['text.usetex'] = True
    rcParams['font.family'] = 'serif'
    rcParams['font.serif'] = 'cm'
    rcParams['axes.titlesize'] = 'medium'
    rcParams['font.serif'] = 'cm'

    ga = 1.4
    N = np.logspace(0,5,14).astype(int)
    xt = 10**np.array(np.arange(6))
    print(xt)

    Xmax = 1.4
    Xmin = 1.04

    dt_native = []
    dt_fort = []
    for N1 in N:
        N2 = 1000
        X = np.random.rand(N1)*(Xmax-Xmin) + Xmin
        T_fort = timeit.Timer('cf.mcpTo_APo_from_Ma(X,ga)','from __main__ import X,ga,cf')
        T_native = timeit.Timer('cf.native.mcpTo_APo_from_Ma(X,ga)','from __main__ import X,ga,cf')
        res_fort = T_fort.autorange()
        res_native = T_native.autorange()
        dt_fort.append(res_fort[1]/res_fort[0])
        dt_native.append(res_native[1]/res_native[0])
        print(N1)

    f,a  = plt.subplots()
    f.set_size_inches((4., 3.))
    a.loglog(N,dt_native,label='Native')
    a.loglog(N,dt_fort,label='Fortran')
    a.set_xlabel(r'Array Length, $n$')
    a.set_ylabel(r'Time per call, $\Delta t$/s')
    a.grid(True)
    a.set_title('Benchmark for evaluation of $\dot{m}\sqrt{c_pT_0}/Ap_0$')
    a.legend()
    a.set_xlim((1,1e5))
    a.set_ylim((1e-7,1e-2))
    a.set_xticks(xt)
    f.tight_layout(pad=0.1)
    plt.savefig('bench_forward.svg')

    speedup = np.array(dt_fort)/np.array(dt_native)
    print(speedup)

    # Initialise lookup table
    cf.to_Ma('mcpTo_APo',0.4,ga,use_lookup=True)

    Xmax = 0.1
    Xmin = 1.1

    dt_native = []
    dt_fort = []
    dt_lookup = []
    for N1 in N:
        X = np.random.rand(N1)*(Xmax-Xmin) + Xmin
        T_fort = timeit.Timer('cf.Ma_from_mcpTo_APo(X,ga)','from __main__ import X,ga,cf')
        T_native = timeit.Timer('cf.native.Ma_from_mcpTo_APo(X,ga)','from __main__ import X,ga,cf')
        T_lookup = timeit.Timer('cf.to_Ma("mcpTo_APo",X,ga,use_lookup=True)','from __main__ import X,ga,cf')
        res_fort = T_fort.autorange()
        res_native = T_native.autorange()
        res_lookup = T_lookup.autorange()
        dt_fort.append(res_fort[1]/res_fort[0])
        dt_native.append(res_native[1]/res_native[0])
        dt_lookup.append(res_lookup[1]/res_lookup[0])
        print(N1)

    f,a  = plt.subplots()
    f.set_size_inches((4., 3.))
    a.loglog(N,dt_native,label='Native')
    a.loglog(N,dt_fort,label='Fortran')
    a.loglog(N,dt_lookup,label='Lookup')
    a.set_xlabel(r'Array Length, $n$')
    a.set_ylabel(r'Time per call, $\Delta t$/s')
    a.set_title('Benchmark for inversion of $\dot{m}\sqrt{c_pT_0}/Ap_0$')
    a.grid(True)
    a.legend()
    a.set_xlim((1,1e5))
    a.set_xticks(xt)
    a.set_ylim((1e-7,1e-1))
    f.tight_layout(pad=0.1)
    plt.show()
    plt.savefig('bench_inverse.svg')

    speedup = np.array(dt_fort)/np.array(dt_native)
    print(speedup)
    speedup = np.array(dt_lookup)/np.array(dt_native)
    print(speedup)


