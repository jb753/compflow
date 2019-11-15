import compflow as cf
import numpy as np


# Load validation data from CUED Thermofluids Data Book
ga = 1.4
# Subsonic
var_sub = ['Ma', 'To_T', 'Po_P', 'rhoo_rho', 'V_cpTo', 'mcpTo_APo', 'mcpTo_AP', 'F_mcpTo', '4cfL_D', 'rhoVsq_2Po']
dat_sub = {var_sub[i]: np.loadtxt('compflow/test/cued_data_book_subsonic.csv').reshape((-1, 10))[1:, i]
           for i in range(len(var_sub))}
# Supersonic
var_sup = var_sub + ['Mash', 'Posh_Po', 'Psh_P', 'Posh_P', 'Tsh_T', 'nu', 'Ma2']
dat_sup = {var_sup[i]: np.loadtxt('compflow/test/cued_data_book_supersonic.csv').reshape((-1, 17))[1:, i]
           for i in range(len(var_sup))}

# Take reciprocals
for v in ['To_T', 'Po_P', 'rhoo_rho']:
    dat_sub[v] = 1. / dat_sub[v]
    dat_sup[v] = 1. / dat_sup[v]

# All variables which the module implements
var_test_sub = ['To_T', 'Po_P', 'rhoo_rho', 'V_cpTo', 'mcpTo_APo', 'mcpTo_AP']
var_test_sup = ['To_T', 'Po_P', 'rhoo_rho', 'V_cpTo', 'mcpTo_APo', 'mcpTo_AP', 'Mash', 'Posh_Po']

# Subsonic tests
for v in var_test_sub:
    #
    # Error w.r.t. table looking up Ma
    err_forward = np.abs(np.round(cf.from_Ma(v, dat_sub['Ma'], ga), 4) - dat_sub[v])
    assert np.max(err_forward) < 2.e-4

    # Error w.r.t. table finding Ma
    err_back = np.abs(cf.to_Ma(v, dat_sub[v], ga) - dat_sub['Ma'])
    assert np.max(err_back) < 5.e-3

# Supersonic tests
for v in var_test_sup:
    #
    # Error w.r.t. table looking up Ma
    err_forward = np.abs(cf.from_Ma(v, dat_sup['Ma'], ga) / dat_sup[v] - 1.)
    assert np.max(err_forward) < 0.005

    # Error w.r.t. table finding Ma
    err_back = np.abs(cf.to_Ma(v, dat_sup[v], ga, supersonic=True) - dat_sup['Ma'])
    assert np.sqrt(np.mean(err_back**2.) < 0.001)
    assert np.max(err_back) < 0.05

