import numpy as np
import timeit
# from numba import njit
from compflow_fort import mcpto_apo as func_fort
from compflow import Ma_from_mcpTo_APo as func_pure

if __name__ == '__main__':

    N1 = 10
    N2 = 100
    # Ma = np.linspace(0.1,1.4,N1)
    Ma = np.array([1.1])
    ga = 1.4

    print(func_fort(Ma,ga,True))
    print(func_fort(Ma,ga,False))
    print(func_pure(Ma,ga,True))
    print(func_pure(Ma,ga,False))

    # for _ in range(N2):
    #     Po_P_fort(Ma,ga)

    for t in ['pure', 'fort']:
        fname = 'func_%s' % t
        print(fname)
        print(timeit.timeit(fname + '(Ma,ga)','from __main__ import %s,Ma,ga' % fname,
                number=N2))


