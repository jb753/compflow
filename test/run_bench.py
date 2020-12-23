import numpy as np
import timeit
# from numba import njit
from compflow_fort import mash as func_fort
from compflow import Ma_from_Mash as func_pure

if __name__ == '__main__':

    N1 = 10
    N2 = 100
    Ma = np.linspace(0.3,0.99,N1)
    # Ma = np.array([0.6])
    ga = 1.4

    print(func_fort(Ma,ga))
    print(func_pure(Ma,ga))

    # for _ in range(N2):
    #     Po_P_fort(Ma,ga)

    for t in ['pure', 'fort']:
        fname = 'func_%s' % t
        print(fname)
        print(timeit.timeit(fname + '(Ma,ga)','from __main__ import %s,Ma,ga' % fname,
                number=N2))


