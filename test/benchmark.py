import compflow
import numpy as np
import numpy.random
#import matplotlib.pyplot as plt


# Input data
N = 1000
Ma_true = np.sort(numpy.random.rand(N))
ga = 1.4
Q = compflow.mcpTo_APo_from_Ma( Ma_true, ga)

y = np.empty_like(Q)
for k in range(100):
    for i in range(N):
        y[i] = compflow.to_Ma('mcpTo_APo',Q[i],ga,use_lookup=True)

err = np.max(np.abs(y - Ma_true))
print(err)
