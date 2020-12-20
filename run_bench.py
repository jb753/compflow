import numpy as np
import timeit
# from numba import njit
from compflow_fort import po_p as func_fort
from compflow import V_cpTo_from_Ma as func_pure

# def Po_P_pure(Ma,ga):
# return (1. + 0.5*(ga - 1.)*Ma**2.)**(ga/(ga-1.0))

# @njit
# def Po_P_numba(Ma,ga):
#     return (1. + 0.5*(ga - 1.)*Ma**2.)**(ga/(ga-1.0))

if __name__ == '__main__':

    N1 = 100
    N2 = 10000
    Ma = np.linspace(0,1,N1)
    Ma = np.array([ 0.3 ])
    ga = 1.4

    # for _ in range(N2):
    #     Po_P_fort(Ma,ga)

    for t in ['pure', 'fort']:
        fname = 'func_%s' % t
        print(fname)
        print(timeit.timeit(fname + '(Ma,ga)','from __main__ import %s,Ma,ga' % fname,
                number=N2))


