import compflow
import numpy.random

N = 100000
ga = 1.4
Q1 = 0.01
Q2 = 1.28
Q = numpy.random.rand(N) * (Q2-Q1) + Q1
compflow.Ma_from_mcpTo_APo_vector(Q,ga)
#for i in range(N):
#    compflow.Ma_from_mcpTo_APo_vector(Q[[i]],ga)
    #compflow.Ma_from_mcpTo_APo_scalar(Q[i],ga)
