#
# Plot the values of the compressible flow quantities as function of Mach
# and inverted values of Mach from each quantity
#
import compflow as cf
import numpy as np
import matplotlib.pyplot as plt

# Input data
Ma = np.linspace(0., 3.)
ga = np.array(1.4)
varlist = ['To_T', 'Po_P', 'rhoo_rho', 'V_cpTo', 'mcpTo_APo', 'mcpTo_AP', 'A_Acrit', 'Mash', 'Posh_Po']

# Forwards calcuations
Y = {v: cf.from_Ma(v, Ma, ga) for v in varlist}

# Backwards calculations in sub and supersonic regimes
Xsub = {v: cf.to_Ma(v, Y[v][(Ma <= 1.) & np.isfinite(Y[v])], ga) for v in varlist}
Xsup = {v: cf.to_Ma(v, Y[v][(Ma > 1.) & np.isfinite(Y[v])], ga, supersonic=True) for v in varlist}

# Create figure
fig, ax = plt.subplots()
fig.set_size_inches((7., 7.))

# Plot forwards
for v in varlist:
    ax.plot(Ma, Y[v], label=v.replace('_', ''))

# Plot inverted subsonic
ax.set_prop_cycle(None)
for v in varlist:
    ax.plot(Xsub[v], Y[v][(Ma <= 1.) & np.isfinite(Y[v])], 'x')

# Plot inverted supersonic
ax.set_prop_cycle(None)
for v in varlist:
    ax.plot(Xsup[v], Y[v][(Ma > 1.) & np.isfinite(Y[v])], '+')

# Decorate plot
ax.legend()
ax.set_xlim(Ma[np.ix_([1, -1])])
ax.set_ylim((0., 3.))
ax.set_title(r'Compressible flow quantities for perfect gas with $\gamma = {}$'.format(ga))
ax.set_xlabel(r'Mach Number, $\mathit{Ma}$')
ax.set_ylabel(r'Quantity')
fig.tight_layout()

plt.show()
