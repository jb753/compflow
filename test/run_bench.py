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

    dt_native = []
    dt_fort = []
    for N1 in N:
        N2 = 1000
        X = np.linspace(1.04,1.4,N1)
        T_fort = timeit.Timer('cf.V_cpTo_from_Ma(X,ga)','from __main__ import X,ga,cf')
        T_native = timeit.Timer('cf.native.V_cpTo_from_Ma(X,ga)','from __main__ import X,ga,cf')
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
    a.set_title('Benchmark for $V/\sqrt{c_pT_0}$')
    a.legend()
    f.tight_layout()
    plt.savefig('bench_V_cpTo.svg')

    # Initialise lookup table
    cf.to_Ma('mcpTo_APo',0.4,ga,use_lookup=True)

    dt_native = []
    dt_fort = []
    dt_lookup = []
    for N1 in N:
        N2 = 1000
        X = np.linspace(0.1,1.1,N1)
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
    a.loglog(N,dt_lookup,label='Lookup')
    a.loglog(N,dt_fort,label='Fortran')
    a.set_xlabel(r'Array Length, $n$')
    a.set_ylabel(r'Time per call, $\Delta t$/s')
    a.set_title('Benchmark for $\dot{m}\sqrt{c_pT_0}/Ap_0$')
    a.grid(True)
    a.legend()
    f.tight_layout()
    plt.savefig('bench_mcpTo_APo.svg')


