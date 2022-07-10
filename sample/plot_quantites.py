#
# Plot the values of the compressible flow quantities as function of Mach
# and inverted values of Mach from each quantity
#
import compflow as cf
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams

# Input data
Ma = np.linspace(0., 4.,1000)
ga = 1.4
varlist = ['To_T', 'Po_P', 'rhoo_rho', 'V_cpTo', 'mcpTo_APo',
           'mcpTo_AP', 'A_Acrit', 'Mash', 'Posh_Po','F_mcpTo']
labels = [r'$T_0/T$', r'$p_0/p$', r'$\rho_0/\rho$', r'$V/\sqrt{c_pT_0}$',
          r'$\dot{m}\sqrt{c_pT_0}/Ap_0$',
          r'$\dot{m}\sqrt{c_pT_0}/Ap$', r'$A/A_*$',
          r'$\mathit{M\kern-.2ema}_\mathrm{sh}$', r'$p_{0,\mathrm{sh}}/p_0$',
          r'$F/\dot{m}\sqrt{c_pT_0}$']

# Forwards calcuations
Y = {v: cf.from_Ma(v, Ma, ga) for v in varlist}

# Remove non-physical values
for v in ['Mash','Posh_Po']:
    Y[v][Ma<1.] = np.nan

# # Backwards calculations in sub and supersonic regimes
# Xsub = {v: cf.to_Ma(v, Y[v][(Ma <= 1.) & np.isfinite(Y[v])], ga)
#         for v in varlist}
# Xsup = {v: cf.to_Ma(v, Y[v][(Ma > 1.) & np.isfinite(Y[v])], ga, supersonic=True)
#         for v in varlist}

# Set up plotting
rcParams['text.usetex'] = True
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = 'cm'
rcParams['axes.titlesize'] = 'medium'
rcParams['font.serif'] = 'cm'

# Create figure
fig, ax = plt.subplots()
fig.set_size_inches((7., 4.326))

# Plot forwards
for v, lv in zip(varlist, labels):
    ax.plot(Ma, Y[v], label=lv)

# # Plot inverted subsonic
# ax.set_prop_cycle(None)
# for v in varlist:
#     ax.plot(Xsub[v], Y[v][(Ma <= 1.) & np.isfinite(Y[v])], 'x')

# # Plot inverted supersonic
# ax.set_prop_cycle(None)
# for v in varlist:
#     ax.plot(Xsup[v], Y[v][(Ma > 1.) & np.isfinite(Y[v])], '+')

# Decorate plot
ax.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
ax.set_xlim(Ma[np.ix_([1, -1])])
ax.set_xticks([0.,1.,2.,3.,4.])
ax.set_yticks([0.,1.,2.,3.,4.])
ax.set_ylim((0., 3.))
ax.set_title(
    r'\noindent Compressible flow relations\\for a perfect gas with $\gamma = {}$'.format(ga))
ax.set_xlabel(r'Mach Number, $\mathit{M\kern-.2ema}$')
ax.set_ylabel(r'Non-dimensional Group')
fig.tight_layout(pad=0.2)

plt.savefig('sample.png',dpi=250)
