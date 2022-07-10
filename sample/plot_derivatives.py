#
# Plot the derivatives of compressible flow quantities as function of Mach
# calculated from explicit functions and a numerical approximation
#
import compflow as cf
import numpy as np
import matplotlib.pyplot as plt

# Input data
Ma = np.linspace(0., 3.)
ga = np.array(1.4)
varlist = ['To_T', 'Po_P', 'rhoo_rho', 'V_cpTo', 'mcpTo_APo', 'mcpTo_AP', 'A_Acrit', 'Mash', 'Posh_Po','F_mcpTo']

# Set up figure
fig, ax = plt.subplots()
fig.set_size_inches((7., 7.))

# Evaluate explicit derivative functions
for v in varlist:
    ax.plot(Ma, cf.derivative_from_Ma(v, Ma, ga), label=v.replace('_', ''))

# Numerical approximation
ax.set_prop_cycle(None)
for v in varlist:
    Y = cf.from_Ma(v, Ma, ga)
    dYdMa = np.gradient(Y, Ma)
    ax.plot(Ma, dYdMa, 'x')

# Decorate plot
ax.set_xlim(Ma[np.ix_([1, -1])])
ax.set_ylim((-10., 20.))
ax.legend()
ax.set_title(r'Compressible flow derivatives for perfect gas with $\gamma = {}$'.format(ga))
ax.set_xlabel(r'Mach Number, $\mathit{Ma}$')
ax.set_ylabel(r'Derivative')
fig.tight_layout()
plt.show()
