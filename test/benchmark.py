import compflow
import numpy as np
import numpy.random
#import matplotlib.pyplot as plt


# Input data
N = 1000
Ma_true = np.sort(numpy.random.rand(N))
ga = 1.4
Q = compflow.mcpTo_APo_from_Ma( Ma_true, ga)

# vectorised newton
#y = compflow.Ma_from_mcpTo_APo(Q,ga)
#err = np.max(np.abs(y - Ma_true))
#print(err)

# for loop
# table lookup
f = compflow.generate_lookup('mcpTo_APo',ga)
#y = f(Q)
y = np.empty_like(Q)
for k in range(100):
    for i in range(N):
        # vectorised newton, one at a time
        #compflow.Ma_from_mcpTo_APo_vector(Q[[i]],ga)
        # root_scalar
        #y[i] = compflow.Ma_from_mcpTo_APo_scalar(Q[i],ga)
        y[i] = f(Q[i])

err = np.max(np.abs(y - Ma_true))
print(err)

#
#print(err)
